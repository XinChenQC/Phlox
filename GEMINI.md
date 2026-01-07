# Gemini Code Assistant Context

This document provides context for the Phlox project to the Gemini code assistant.

## Project Overview

Phlox is a Python-based tool for discovering reaction mechanisms from reactive molecular dynamics (MD) trajectories. It analyzes trajectory files, identifies molecular fragments, tracks their evolution over time, and constructs a reaction network.

The project is structured as a Python package with a command-line interface (CLI) and a Python API. It supports parallel processing to handle large trajectories efficiently.

**Key Technologies:**

*   **Programming Language:** Python 3.8+
*   **Package Management:** `setuptools`, `pip`
*   **Dependencies:**
    *   `numpy`: For numerical operations.
    *   `scipy`: For scientific computing.
    *   `networkx`: For building and managing the reaction network graph.
    *   `rdkit`: For cheminformatics tasks like SMILES generation.
    *   `pydot`: For exporting the reaction network graph.
*   **Development Tools:**
    *   `pytest`: For testing.
    *   `black`, `isort`: For code formatting.
    *   `flake8`, `mypy`: For linting and static analysis.

**Architecture:**

The core logic is encapsulated in the `Analyzer` class (`src/phlox/base/analyzer.py`). This class manages the analysis workflow, including:

1.  **Configuration:** Reads analysis parameters from a configuration file.
2.  **Trajectory Reading:** Reads trajectory data from XYZ or LAMMPS files.
3.  **Parallel Processing:** Distributes the analysis workload across multiple CPU cores using the `multiprocessing` module. Worker processes (`_Worker` class) analyze chunks of the trajectory.
4.  **Frame Processing:** Identifies molecules and their structures in each frame.
5.  **Network Building:** Compares consecutive frames to detect reactions and builds a reaction network using `networkx`.
6.  **Reaction Classification:** Classifies reactions based on the changes in molecular structures.
7.  **Output Generation:** Writes the results, including the reaction network, species information, and statistics, to the output directory.

The command-line interface is defined in `src/phlox/cli.py` using the `argparse` module. It provides two main commands: `analyze` to run the analysis and `config` to generate an example configuration file.

## Building and Running

The project uses a `Makefile` to simplify common development tasks.

**Installation:**

*   To install the package in editable mode:
    ```bash
    make install
    ```
*   To install with development dependencies:
    ```bash
    make install-dev
    ```

**Running the Application:**

The application can be run as a command-line tool.

1.  **Generate a configuration file:**
    ```bash
    phlox config
    ```
2.  **Run the analysis:**
    ```bash
    phlox analyze input.conf
    ```

**Running Tests:**

*   To run the test suite:
    ```bash
    make test
    ```

**Linting and Formatting:**

*   To check for code style and type errors:
    ```bash
    make lint
    ```
*   To format the code:
    ```bash
    make format
    ```

## Development Conventions

*   **Code Style:** The project uses `black` for code formatting and `isort` for import sorting.
*   **Testing:** Tests are written using `pytest` and are located in the `tests/` directory.
*   **Type Hinting:** The project uses type hints, and `mypy` is used for static type checking.
*   **CLI:** The command-line interface is built with `argparse`.
*   **Configuration:** The application is configured through a custom configuration file format.
*   **Parallelism:** The application is designed to run in parallel using the `multiprocessing` module.
