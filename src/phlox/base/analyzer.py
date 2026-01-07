"""
Main analyzer class - orchestrates the entire analysis workflow.
"""

import logging
import time
import copy
from dataclasses import dataclass
from typing import Optional, List, Dict, Any
from pathlib import Path
from multiprocessing import Pool
from collections import Counter
import networkx as nx

from phlox.base.config import Config
from phlox.base.results import Results
from phlox.base.atom import Atom
from phlox.base.molecule import Molecule


logger = logging.getLogger(__name__)


@dataclass
class WorkerTask:
    """
    Work assignment for a single worker process.

    Attributes:
        worker_id: Unique worker identifier
        file_start: Index of first trajectory file
        file_end: Index of last trajectory file
        frame_start: Global frame number to start
        frame_end: Global frame number to end
        line_start: Line number to start reading (file format specific)
        line_end: Line number to end reading (-1 = end of file)
    """
    worker_id: int
    file_start: int
    file_end: int
    frame_start: int
    frame_end: int
    line_start: int = 0
    line_end: int = -1


class _Worker:
    """
    Worker process for analyzing trajectory chunks.

    Internal class used by Analyzer for parallel processing.
    """

    def __init__(self, task: WorkerTask, config: Config):
        """
        Initialize worker.

        Args:
            task: Work assignment
            config: Configuration object
        """
        self.task = task
        self.config = config
        self.worker_id = task.worker_id

        # Local storage
        self.local_graph = nx.MultiDiGraph()
        self.struct_dict = {}  # DicStruct: speciesID -> structure info
        self.molecule_info = {}  # DicMoleInfo: moleculeID -> atom indices
        self.time_series = []  # SpeciesCount: [Counter, ...]

        # Frame tracking
        self.current_frame = 0
        self.first_frame_molecules = None
        self.last_frame_molecules = None

        # Frame processor
        from phlox.base.frame_processor import FrameProcessor
        self.frame_processor = FrameProcessor(config, worker_id=task.worker_id)

    @staticmethod
    def process_chunk(task: WorkerTask, config: Config) -> Dict[str, Any]:
        """
        Static method for multiprocessing.Pool.starmap().

        Args:
            task: Work assignment
            config: Configuration

        Returns:
            Dictionary with worker results
        """
        # Add delay to prevent resource contention
        time.sleep(task.worker_id * 30.0)

        worker = _Worker(task, config)
        return worker.run()

    def run(self) -> Dict[str, Any]:
        """
        Execute worker task.

        Returns:
            Dictionary containing worker results
        """
        logger.info(f"Worker {self.worker_id} starting: "
                   f"frames {self.task.frame_start}-{self.task.frame_end}")

        start_time = time.time()

        # Process assigned trajectory chunk
        if self.config.traj_type in ["xyz", "extxyz"]:
            self._process_xyz_chunk()
        elif self.config.traj_type == "lmp":
            self._process_lammps_chunk()
        else:
            raise ValueError(f"Unsupported trajectory type: {self.config.traj_type}")

        # Clean local network
        # logger.info(f"Worker {self.worker_id} cleaning local network...")
        # self._clean_local_network()  # DISABLED: Skip cleaning to preserve all reaction events

        elapsed = time.time() - start_time
        logger.info(f"Worker {self.worker_id} completed in {elapsed:.2f}s: "
                   f"{self.current_frame} frames processed")

        return {
            'worker_id': self.worker_id,
            'graph': self.local_graph,
            'struct_dict': self.struct_dict,
            'molecule_info': self.molecule_info,
            'time_series': self.time_series,
            'first_frame_molecules': self.first_frame_molecules,
            'last_frame_molecules': self.last_frame_molecules,
            'frame_start': self.task.frame_start,
            'frame_end': self.task.frame_end,
        }

    def _process_xyz_chunk(self):
        """Process XYZ trajectory chunk."""
        from phlox.io.trajectory_reader import create_reader

        files_to_process = self.config.traj_files[
            self.task.file_start:self.task.file_end + 1
        ]

        logger.info(f"Worker {self.worker_id} processing {len(files_to_process)} files")

        # Calculate how many frames to skip before starting
        # This is needed when multiple workers process the same file or when
        # the worker starts in the middle of a file
        frames_before_worker = 0
        if self.task.file_start > 0:
            # Count frames in all files before this worker's file range
            for i in range(self.task.file_start):
                reader = create_reader(self.config.traj_files[i], self.config)
                frames_before_worker += reader.count_frames()

        # Calculate offset within the first file
        frames_to_skip = self.task.frame_start - frames_before_worker

        logger.info(f"Worker {self.worker_id} skipping {frames_to_skip} frames in first file")

        previous_molecules = None
        frame_count = 0

        for file_idx, traj_file in enumerate(files_to_process):
            with create_reader(traj_file, self.config) as reader:
                # Skip frames in the first file only
                if file_idx == 0 and frames_to_skip > 0:
                    reader.skip_frames(frames_to_skip)
                    logger.debug(f"Worker {self.worker_id}: skipped {frames_to_skip} frames")

                while True:
                    frame_data = reader.read_frame()
                    if frame_data is None:
                        break

                    atoms, pbc_box = frame_data

                    # Check if we should stop
                    global_frame = self.task.frame_start + frame_count
                    if global_frame > self.task.frame_end:
                        break

                    # Process frame using shared FrameProcessor
                    molecules, frame_species_ids = self.frame_processor.process_frame(
                        atoms, pbc_box, global_frame,
                        self.struct_dict, self.molecule_info
                    )

                    # Add frame counts to time series as Counter
                    self.time_series.append(Counter(frame_species_ids))

                    # Store first frame
                    if frame_count == 0:
                        self.first_frame_molecules = copy.deepcopy(molecules)

                    # Compare with previous frame
                    if previous_molecules is not None:
                        self._compare_frames(previous_molecules, molecules, global_frame)

                    previous_molecules = copy.deepcopy(molecules)
                    self.last_frame_molecules = molecules

                    frame_count += 1

                    # Update neighbor lists periodically (every 5 frames)
                    if frame_count % 100 == 0 and self.frame_processor.fragmenter is not None:
                        self.frame_processor.update_neighbor_lists(atoms, pbc_box)

                    if frame_count % 100 == 0:
                        logger.debug(f"Worker {self.worker_id}: processed {frame_count} frames")

        self.current_frame = frame_count

    def _process_lammps_chunk(self):
        """Process LAMMPS trajectory chunk."""
        from phlox.io.trajectory_reader import create_reader

        files_to_process = self.config.traj_files[
            self.task.file_start:self.task.file_end + 1
        ]

        logger.info(f"Worker {self.worker_id} processing {len(files_to_process)} files")

        # Calculate how many frames to skip before starting
        # This is needed when multiple workers process the same file or when
        # the worker starts in the middle of a file
        frames_before_worker = 0
        if self.task.file_start > 0:
            # Count frames in all files before this worker's file range
            for i in range(self.task.file_start):
                reader = create_reader(self.config.traj_files[i], self.config)
                frames_before_worker += reader.count_frames()

        # Calculate offset within the first file
        frames_to_skip = self.task.frame_start - frames_before_worker

        logger.info(f"Worker {self.worker_id} skipping {frames_to_skip} frames in first file")

        previous_molecules = None
        frame_count = 0

        for file_idx, traj_file in enumerate(files_to_process):
            with create_reader(traj_file, self.config) as reader:
                # Skip frames in the first file only
                if file_idx == 0 and frames_to_skip > 0:
                    reader.skip_frames(frames_to_skip)
                    logger.debug(f"Worker {self.worker_id}: skipped {frames_to_skip} frames")

                while True:
                    frame_data = reader.read_frame()
                    if frame_data is None:
                        break

                    atoms, pbc_box = frame_data

                    # Check if we should stop
                    global_frame = self.task.frame_start + frame_count
                    if global_frame > self.task.frame_end:
                        break

                    # Process frame using shared FrameProcessor
                    molecules, frame_species_ids = self.frame_processor.process_frame(
                        atoms, pbc_box, global_frame,
                        self.struct_dict, self.molecule_info
                    )

                    # Add frame counts to time series as Counter
                    self.time_series.append(Counter(frame_species_ids))

                    # Store first frame
                    if frame_count == 0:
                        self.first_frame_molecules = copy.deepcopy(molecules)

                    # Compare with previous frame
                    if previous_molecules is not None:
                        self._compare_frames(previous_molecules, molecules, global_frame)

                    previous_molecules = copy.deepcopy(molecules)
                    self.last_frame_molecules = molecules

                    frame_count += 1

                    # Update neighbor lists periodically (every 5 frames)
                    if frame_count % 100 == 0 and self.frame_processor.fragmenter is not None:
                        self.frame_processor.update_neighbor_lists(atoms, pbc_box)

                    if frame_count % 100 == 0:
                        logger.debug(f"Worker {self.worker_id}: processed {frame_count} frames")

        self.current_frame = frame_count

    def _compare_frames(
        self,
        molecules_prev: List[Molecule],
        molecules_curr: List[Molecule],
        frame_number: int
    ):
        """Compare consecutive frames to detect reactions."""
        from phlox.network.builder import NetworkBuilder

        logger.debug(f"Worker {self.worker_id}: comparing frames at {frame_number}")

        # Initialize network builder on first call
        if not hasattr(self, 'network_builder'):
            self.network_builder = NetworkBuilder(self.config)
            self.local_graph = self.network_builder.graph

        # Compare frames and build network
        self.network_builder.compare_frames(
            molecules_prev,
            molecules_curr,
            frame_number
        )

        # Update local graph reference
        self.local_graph = self.network_builder.graph

    def _clean_local_network(self):
        """Clean local reaction network."""
        from phlox.filter.oscillation import OscillationRemover
        from phlox.filter.contraction import NodeContractor

        logger.debug(f"Worker {self.worker_id}: cleaning network before "
                    f"({self.local_graph.number_of_nodes()} nodes, "
                    f"{self.local_graph.number_of_edges()} edges)")

        # Remove oscillations with relaxed threshold for workers
        oscillation_lag = max(1, self.config.stability_lag // 2)

        remover = OscillationRemover(stability_threshold=oscillation_lag)
        self.local_graph = remover.remove_oscillations(self.local_graph)

        # Contract nodes
        contractor = NodeContractor(stability_threshold=oscillation_lag)
        self.local_graph = contractor.contract_nodes(self.local_graph)

        logger.debug(f"Worker {self.worker_id}: cleaning network after "
                    f"({self.local_graph.number_of_nodes()} nodes, "
                    f"{self.local_graph.number_of_edges()} edges)")


class Analyzer:
    """
    Main analyzer for reactive MD trajectories.

    Orchestrates the full analysis pipeline:
    1. Read trajectory files
    2. Track molecular fragments
    3. Build reaction network
    4. Clean network
    5. Classify reactions
    6. Generate outputs
    """

    def __init__(self, config: Config):
        """
        Initialize analyzer.

        Args:
            config: Configuration object
        """
        self.config = config
        self.config.validate()

        # Create output directory
        self.config.output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Initialized analyzer with {config.n_cores} cores")
        logger.info(f"Trajectory type: {config.traj_type}")
        logger.info(f"Number of files: {len(config.traj_files)}")

    def run(self) -> Results:
        """
        Run complete analysis pipeline.

        Returns:
            Results object containing all analysis outputs

        Raises:
            RuntimeError: If analysis fails
        """
        logger.info("Starting Phlox analysis...")

        # Create results object
        results = Results(config=self.config)

        try:
            if self.config.n_cores == 1:
                logger.info("Running serial analysis...")
                self._run_serial(results)
            else:
                logger.info(f"Running parallel analysis with {self.config.n_cores} cores...")
                self._run_parallel(results)

            # Compute SMILES for all species in batch (ReaxANA approach)
            logger.info("Computing SMILES for all species...")
            self._compute_smiles(results)

            # Clean network
            # logger.info("Cleaning reaction network...")
            # self._clean_network(results)  # DISABLED: Skip cleaning to preserve all reaction events

            # Classify reactions
            logger.info("Classifying reactions...")
            self._classify_reactions(results)

            # Write outputs
            logger.info("Writing results...")
            self._write_outputs(results)

            logger.info("Analysis complete!")
            return results

        except Exception as e:
            logger.error(f"Analysis failed: {e}")
            raise RuntimeError(f"Analysis failed: {e}") from e

    def _run_serial(self, results: Results):
        """Run serial analysis (single core)."""
        from collections import Counter
        from phlox.base.frame_processor import FrameProcessor
        from phlox.network.builder import NetworkBuilder
        from phlox.io.trajectory_reader import XYZReader, LAMMPSReader

        frame_processor = FrameProcessor(self.config)
        network_builder = NetworkBuilder(self.config)

        # Read trajectory
        ReaderClass = XYZReader if self.config.traj_type in ["xyz", "extxyz"] else LAMMPSReader

        for traj_file in self.config.traj_files:
            with ReaderClass(str(traj_file), self.config) as reader:
                frame_num = 0
                molecules_prev = None

                while True:
                    try:
                        result = reader.read_frame()
                        if result is None:
                            break
                        atoms, pbc_box = result
                    except StopIteration:
                        break

                    frame_num += 1

                    # Process frame using shared FrameProcessor
                    molecules_curr, frame_species_ids = frame_processor.process_frame(
                        atoms, pbc_box, frame_num,
                        results.struct_dict, results.molecule_info
                    )

                    # Add frame counts to time series as Counter
                    results.time_series.append(Counter(frame_species_ids))

                    # Compare frames
                    if molecules_prev is not None:
                        network_builder.compare_frames(molecules_prev, molecules_curr, frame_num)

                    molecules_prev = molecules_curr

                    if frame_num % 100 == 0:
                        logger.info(f"Processed {frame_num} frames")

        results.graph = network_builder.graph
        results.total_frames = frame_num
        logger.info(f"Built network with {results.graph.number_of_nodes()} nodes and {results.graph.number_of_edges()} edges")

    def _run_parallel(self, results: Results):
        """
        Run parallel analysis (multi-core).

        Based on parallelPool() from AnalysFlow_parallel.py
        """
        from phlox.io.trajectory_reader import create_reader
        from phlox.network.builder import NetworkBuilder

        # Count frames
        logger.info("Counting frames in trajectory files...")
        total_frames = 0
        frames_per_file = []

        for traj_file in self.config.traj_files:
            logger.debug(f"Counting frames in {traj_file}")
            reader = create_reader(traj_file, self.config)
            n_frames = reader.count_frames()
            frames_per_file.append(n_frames)
            total_frames += n_frames
            logger.debug(f"  {n_frames} frames")

        logger.info(f"Total frames: {total_frames}")

        # Distribute work across workers
        work_assignments = self._distribute_work(total_frames, frames_per_file)
        logger.info(f"Work distributed across {len(work_assignments)} workers")

        # Execute workers in parallel
        logger.info(f"Starting parallel analysis with {self.config.n_cores} cores")
        logger.info("Spawning worker processes...")

        with Pool(processes=self.config.n_cores) as pool:
            worker_args = [(task, self.config) for task in work_assignments]
            worker_results = pool.starmap(_Worker.process_chunk, worker_args)

        logger.info("All workers completed")

        # Sort results by worker_id
        worker_results.sort(key=lambda x: x['worker_id'])

        # Merge results from all workers
        logger.info("Merging results from workers...")
        for worker_result in worker_results:
            worker_id = worker_result['worker_id']
            logger.debug(f"Merging results from worker {worker_id}")

            # Merge structure dictionary
            if 'struct_dict' in worker_result:
                results.struct_dict.update(worker_result['struct_dict'])

            # Merge molecule info dictionary
            if 'molecule_info' in worker_result:
                results.molecule_info.update(worker_result['molecule_info'])

            # Merge time series
            if 'time_series' in worker_result:
                results.time_series.extend(worker_result['time_series'])

            # Merge reaction graphs
            if 'graph' in worker_result:
                graph = worker_result['graph']
                results.graph.add_nodes_from(list(graph.nodes(data=True)))
                results.graph.add_edges_from(list(graph.edges(data=True)))

        logger.info(f"Merged {results.graph.number_of_nodes()} nodes, "
                   f"{results.graph.number_of_edges()} edges")

        # Bridge transitions between workers
        logger.info("Bridging transitions between workers...")
        builder = NetworkBuilder(self.config)
        builder.graph = results.graph

        for i in range(len(worker_results) - 1):
            worker_n = worker_results[i]
            worker_n1 = worker_results[i + 1]

            if 'last_frame_molecules' not in worker_n or 'first_frame_molecules' not in worker_n1:
                continue

            last_molecules = worker_n['last_frame_molecules']
            first_molecules = worker_n1['first_frame_molecules']

            if last_molecules is None or first_molecules is None:
                continue

            transition_time = worker_n['frame_end']
            logger.debug(f"Bridging workers {i} and {i+1} at frame {transition_time}")

            builder.compare_frames(last_molecules, first_molecules, transition_time)
            results.graph = builder.graph

        logger.info(f"Bridging complete: {results.graph.number_of_nodes()} nodes, "
                   f"{results.graph.number_of_edges()} edges")

        # Store results
        results.total_frames = total_frames

    def _distribute_work(self, total_frames: int, frames_per_file: List[int]) -> List[WorkerTask]:
        """
        Distribute trajectory frames across workers.

        Strategy:
        - First worker gets 20% more frames (warm-up)
        - Last worker gets 20% fewer frames
        - Middle workers get equal shares

        Args:
            total_frames: Total number of frames
            frames_per_file: List of frame counts per file

        Returns:
            List of WorkerTask assignments
        """
        n_cores = self.config.n_cores

        if n_cores == 1:
            # Single core - process everything
            task = WorkerTask(
                worker_id=0,
                file_start=0,
                file_end=len(self.config.traj_files) - 1,
                frame_start=0,
                frame_end=total_frames - 1,
                line_start=0,
                line_end=-1
            )
            return [task]

        # Calculate frames per core
        base_frames = total_frames // n_cores
        first_core_frames = int(base_frames * 1.2)
        last_core_frames = int(base_frames * 0.8)

        # Distribute frames
        frames_per_core = []
        for i in range(n_cores - 1):
            frames = first_core_frames + int(i * (last_core_frames - first_core_frames) / (n_cores - 1))
            frames_per_core.append(frames)

        # Last core gets remaining frames
        last_core_frames = total_frames - sum(frames_per_core)
        frames_per_core.append(last_core_frames)

        logger.info(f"Frames per core: {frames_per_core}")

        # Convert to file/frame assignments
        work_assignments = []
        current_frame = 0
        cumulative_frames = 0
        current_file = 0

        for worker_id in range(n_cores):
            start_frame = current_frame
            end_frame = current_frame + frames_per_core[worker_id] - 1

            # Find which files this worker needs
            file_start = current_file
            file_end = current_file

            frames_needed = frames_per_core[worker_id]
            frames_counted = 0

            temp_file = current_file
            temp_cumulative = cumulative_frames

            while frames_counted < frames_needed and temp_file < len(frames_per_file):
                available = frames_per_file[temp_file] - (temp_cumulative - sum(frames_per_file[:temp_file]))
                if frames_counted + available <= frames_needed:
                    frames_counted += available
                    file_end = temp_file
                    temp_file += 1
                    temp_cumulative += available
                else:
                    file_end = temp_file
                    frames_counted = frames_needed
                    temp_cumulative += (frames_needed - frames_counted)

            task = WorkerTask(
                worker_id=worker_id,
                file_start=file_start,
                file_end=file_end,
                frame_start=start_frame,
                frame_end=end_frame,
                line_start=0,
                line_end=-1
            )

            work_assignments.append(task)
            current_frame = end_frame + 1
            cumulative_frames = temp_cumulative
            current_file = temp_file

        logger.info(f"Created {len(work_assignments)} work assignments")
        return work_assignments

    def _compute_smiles(self, results: Results):
        """
        Compute SMILES for all species in batch.

        This is done AFTER trajectory processing to avoid redundant SMILES
        computation for the same species appearing in multiple frames.

        Based on printUnknowStruc() from ReaxANA/tool.py
        """
        from phlox.utils.batch_smiles import compute_species_smiles

        if not results.struct_dict:
            logger.warning("No species to compute SMILES for")
            return

        # Compute SMILES for all species (struct_dict keyed by speciesID)
        count = compute_species_smiles(results.struct_dict)
        logger.info(f"Computed SMILES for {count} species")

    def _clean_network(self, results: Results):
        """Clean reaction network (remove oscillations, contract nodes)."""
        from phlox.filter.oscillation import OscillationRemover
        from phlox.filter.contraction import NodeContractor

        if results.graph is None or results.graph.number_of_nodes() == 0:
            logger.warning("No network to clean")
            return

        # Remove oscillations
        osc_remover = OscillationRemover(stability_threshold=self.config.stability_lag)
        results.graph = osc_remover.remove_oscillations(results.graph)

        # Contract nodes
        contractor = NodeContractor(stability_threshold=self.config.stability_lag)
        results.graph = contractor.contract_nodes(results.graph)

        logger.info(f"Cleaned network: {results.graph.number_of_nodes()} nodes, {results.graph.number_of_edges()} edges")

    def _classify_reactions(self, results: Results):
        """Classify reactions into types."""
        from phlox.network.classifier import ReactionClassifier

        if results.graph is None or results.graph.number_of_nodes() == 0:
            logger.warning("No network to classify")
            return

        classifier = ReactionClassifier(stability_threshold=self.config.stability_lag)
        results.reactions = classifier.classify_all(results.graph)
        logger.info(f"Classified {len(results.reactions)} reaction types")

    def _write_outputs(self, results: Results):
        """Write all output files."""
        import networkx as nx

        if results.graph is None or results.graph.number_of_nodes() == 0:
            logger.warning("No network to write")
            return

        # Write network as DOT file
        output_file = self.config.output_dir / "reaction_network.dot"
        nx.drawing.nx_pydot.write_dot(results.graph, str(output_file))
        logger.info(f"Wrote network to {output_file}")

        # Write basic stats
        stats_file = self.config.output_dir / "stats.txt"
        with open(stats_file, 'w') as f:
            f.write(f"Total frames: {results.total_frames}\n")
            f.write(f"Network nodes: {results.graph.number_of_nodes()}\n")
            f.write(f"Network edges: {results.graph.number_of_edges()}\n")
            if results.reactions:
                f.write(f"Reaction types: {len(results.reactions)}\n")
                total_events = sum(len(v) for v in results.reactions.values())
                f.write(f"Total reaction events: {total_events}\n")
        logger.info(f"Wrote stats to {stats_file}")

        # Write species files
        self._write_species(results)

    def _write_species(self, results: Results):
        """Write species XYZ files and species list."""
        if not results.struct_dict:
            logger.warning("No species to write")
            return

        # Create species directory
        species_dir = self.config.output_dir / "species"
        species_dir.mkdir(parents=True, exist_ok=True)

        # Write species list file
        species_list_file = species_dir / "species_list.txt"

        logger.info(f"Writing {len(results.struct_dict)} species to {species_dir}")

        with open(species_list_file, 'w') as list_file:
            # struct_dict is keyed by speciesID (molecular structure hash)
            for species_id, species_data in results.struct_dict.items():
                # Get species data
                elements = species_data['elements']
                coordinates = species_data['coordinates']
                smiles = species_data.get('smiles', '')

                # Write XYZ file (named by speciesID)
                xyz_file = species_dir / f"{species_id}.xyz"
                with open(xyz_file, 'w') as f:
                    # Write number of atoms
                    f.write(f"{len(elements)}\n")
                    # Write comment line (speciesID and SMILES)
                    f.write(f"{species_id} {smiles}\n")
                    # Write atom lines
                    for elem, coord in zip(elements, coordinates):
                        f.write(f"{elem:2s} {coord[0]:12.6f} {coord[1]:12.6f} {coord[2]:12.6f}\n")

                # Write to species list
                list_file.write(f"species-{species_id}; {smiles}\n")

        logger.info(f"Wrote species XYZ files and species list to {species_dir}")
