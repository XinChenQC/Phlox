"""
Command-line interface for Phlox.
"""

import argparse
import sys
import logging
from pathlib import Path

from phlox import __version__
from phlox.io import read_config
from phlox.base.analyzer import Analyzer


def setup_logging(level: str = "INFO"):
    """Configure logging."""
    numeric_level = getattr(logging, level.upper(), None)
    if not isinstance(numeric_level, int):
        numeric_level = logging.INFO

    logging.basicConfig(
        level=numeric_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Phlox - Reaction mechanism discovery from reactive MD trajectories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze trajectory with configuration file
  phlox analyze input.conf

  # Specify custom output directory
  phlox analyze input.conf --output results/

  # Enable verbose logging
  phlox analyze input.conf --verbose

  # Show version
  phlox --version
        """
    )

    parser.add_argument(
        '--version',
        action='version',
        version=f'Phlox {__version__}'
    )

    subparsers = parser.add_subparsers(dest='command', help='Commands')

    # Analyze command
    analyze_parser = subparsers.add_parser(
        'analyze',
        help='Analyze reactive MD trajectory'
    )
    analyze_parser.add_argument(
        'config',
        type=str,
        help='Path to configuration file'
    )
    analyze_parser.add_argument(
        '--output', '-o',
        type=str,
        default=None,
        help='Output directory (overrides config file)'
    )
    analyze_parser.add_argument(
        '--cores', '-c',
        type=int,
        default=None,
        help='Number of CPU cores (overrides config file)'
    )
    analyze_parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )
    analyze_parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable debug output'
    )

    # Config command
    config_parser = subparsers.add_parser(
        'config',
        help='Generate example configuration file'
    )
    config_parser.add_argument(
        'output',
        type=str,
        nargs='?',
        default='input.conf',
        help='Output configuration file path (default: input.conf)'
    )

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    if args.command == 'analyze':
        # Setup logging
        if args.debug:
            log_level = "DEBUG"
        elif args.verbose:
            log_level = "INFO"
        else:
            log_level = "WARNING"

        setup_logging(log_level)
        logger = logging.getLogger(__name__)

        try:
            # Read configuration
            logger.info(f"Reading configuration from {args.config}")
            config = read_config(args.config)

            # Override with command-line arguments
            if args.output:
                config.output_dir = Path(args.output)
            if args.cores:
                config.n_cores = args.cores

            # Create analyzer
            logger.info("Initializing analyzer")
            analyzer = Analyzer(config)

            # Run analysis
            logger.info("Starting analysis...")
            results = analyzer.run()

            logger.info(f"Analysis complete!")
            logger.info(f"Found {len(results.struct_dict)} unique species")
            logger.info(f"Found {len(results.molecule_info)} molecular instances")
            logger.info(f"Detected {len(results.reactions)} reaction types")
            logger.info(f"Results saved to {config.output_dir}")

            return 0

        except Exception as e:
            logger.error(f"Analysis failed: {e}", exc_info=args.debug)
            return 1

    elif args.command == 'config':
        # Generate example configuration
        from phlox.io import write_config
        from phlox.base.config import Config

        # Create example config
        config = Config(
            n_cores=4,
            traj_type="lmp",
            traj_files=["trajectory.trj"],
            time_per_frame=(0.1, "fs"),
            pbc_molecule=True,
            pbc_box=((0.0, 0.0, 0.0), (100.0, 100.0, 100.0)),
            stability_time=(1.0, "ps"),
            lammps_atom_types={1: "C", 2: "H", 3: "O"},
        )

        write_config(config, args.output)
        print(f"Example configuration written to {args.output}")
        return 0


if __name__ == '__main__':
    sys.exit(main())
