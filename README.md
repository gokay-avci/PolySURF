# Crystal Surface Generator

A comprehensive Rust-based tool for generating crystal surface slabs from CIF files, with support for semantic analysis via MOFid and physical reconstruction.

## Project Tree

```
.
├── Cargo.toml                 # Project configuration and dependencies
├── src
│   ├── lib.rs                 # Library entry point (exports modules)
│   ├── main.rs                # CLI entry point
│   ├── analysis               # Topology analysis
│   │   ├── mod.rs
│   │   └── topology.rs        # Void crawling and safe offset detection
│   ├── chemistry              # Semantic analysis & MOFid integration
│   │   ├── mod.rs
│   │   └── tagging.rs         # Semantic tagging of atoms
│   ├── core                   # Core data structures
│   │   ├── mod.rs
│   │   ├── connectivity.rs    # Molecule finding algorithms
│   │   └── structure.rs       # Crystal, Atom, Lattice definitions
│   ├── io                     # Input/Output
│   │   ├── mod.rs
│   │   ├── parser.rs          # CIF parser
│   │   └── writer.rs          # CIF writer
│   ├── math                   # Mathematical utilities
│   │   ├── mod.rs
│   │   ├── integer_basis.rs
│   │   └── lll.rs             # LLL lattice reduction
│   └── synthesis              # Surface generation logic
│       ├── mod.rs
│       ├── builder.rs         # Slab geometry calculation
│       ├── ionic.rs           # Ionic reconstruction (dipole correction)
│       └── population.rs      # Populating the slab with atoms
├── 1_mofid_to_rust            # Sub-crate for MOFid analysis
├── A_sample_inputs            # Sample input CIF files
└── B_sample_outputs           # Sample output CIF files
```

## Architecture Overview

The project is organized into several modular components:

*   **Core (`src/core`):** Defines the fundamental data structures like `Crystal`, `Atom`, `Lattice`, and `Molecule`. It also handles connectivity analysis to identify discrete molecules within the crystal.
*   **IO (`src/io`):** Handles reading and writing of Crystallographic Information Files (CIF).
*   **Math (`src/math`):** Provides mathematical tools, including integer basis determination and LLL reduction, which are crucial for defining the surface plane.
*   **Analysis (`src/analysis`):** Contains logic for analyzing the crystal topology, such as finding "safe" cut offsets to avoid breaking molecules (Void Crawler).
*   **Synthesis (`src/synthesis`):** The heart of the generation process.
    *   `builder.rs`: Computes the transformation matrix and geometry for the requested (h k l) slab.
    *   `population.rs`: Fills the calculated slab geometry with atoms from the unit cell.
    *   `ionic.rs`: Handles surface reconstruction, such as dipole correction (Tasker III).
*   **Chemistry (`src/chemistry`):** Integrates with MOFid to understand the semantic structure of Metal-Organic Frameworks (nodes, linkers) and allows for chemically-aware surface termination (e.g., exposing nodes or linkers).

## Installation Guide

### Prerequisites

*   **Rust:** Ensure you have the Rust toolchain installed. You can install it from [rustup.rs](https://rustup.rs/).

### Steps

1.  Clone the repository:
    ```bash
    git clone <repository-url>
    cd <repository-directory>
    ```

2.  Build the project:
    ```bash
    cargo build --release
    ```

3.  The binary will be available at `target/release/crystal_surface_generator`.

## Usage

The tool is run via the command line.

### Basic Command

```bash
cargo run --release -- generate --input <INPUT_CIF> --output <OUTPUT_CIF> <H> <K> <L>
```

### Arguments

*   `generate`: The subcommand to generate a surface.
*   `--input, -i`: Path to the input CIF file.
*   `--output, -o`: Path for the output CIF file.
*   `<H> <K> <L>`: The Miller indices of the surface plane (integers).
*   `--thickness`: (Optional) Desired thickness of the slab in Angstroms (default: 15.0).
*   `--vacuum`: (Optional) Thickness of the vacuum layer in Angstroms (default: 15.0).
*   `--offset`: (Optional) Custom cut offset along the normal vector.
*   `--reconstruct`: (Optional) Enable dipole reconstruction (Tasker III).
*   `--with-mofid`: (Optional) Enable MOFid semantic decomposition.
*   `--expose-nodes`: (Optional) Prefer node-terminated surfaces.
*   `--expose-linkers`: (Optional) Prefer linker-terminated surfaces.

## Examples

### 1. Simple Surface Generation

Generate a (1 0 0) surface from `A_sample_inputs/1.cif` with default settings.

```bash
cargo run --release -- generate -i A_sample_inputs/1.cif -o output_100.cif 1 0 0
```

### 2. High-Index Surface with Custom Thickness

Generate a (2 1 4) surface with 20 Å thickness.

```bash
cargo run --release -- generate -i A_sample_inputs/1.cif -o output_214.cif 2 1 4 --thickness 20.0
```

### 3. MOFid Integration

Generate a surface using MOFid to expose nodes.

```bash
cargo run --release -- generate -i A_sample_inputs/1.cif -o output_nodes.cif 1 1 1 --with-mofid --expose-nodes
```

## Testing

The project includes a test suite to verify functionality.

Run the tests using:

```bash
cargo test
```

This will run unit tests and integration tests located in the `tests/` directory.
