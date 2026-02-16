use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::PathBuf;
use std::time::Instant;

use crystal_surface_generator::{
    parser, writer, generate_surface, SurfaceConfig, MoleculeFinder
};

#[derive(Parser)]
#[command(author, version, about = "Ultimate Crystal Surface Generator")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Generates a surface slab from a CIF file.
    Generate {
        #[arg(short, long)]
        input: PathBuf,

        #[arg(short, long)]
        output: PathBuf,
        
        h: i32, 
        k: i32, 
        l: i32,

        #[arg(long, default_value_t = 15.0)]
        thickness: f64,

        #[arg(long, default_value_t = 15.0)]
        vacuum: f64,

        #[arg(long)]
        offset: Option<f64>,

        /// Enable Tasker III dipole reconstruction (Physics).
        #[arg(long)]
        reconstruct: bool,

        /// Enable MOFid semantic decomposition (Chemistry).
        #[arg(long)]
        with_mofid: bool,

        /// Custom directory for MOFid intermediate files.
        #[arg(long)]
        mofid_work_dir: Option<PathBuf>,

        /// Prefer node-terminated surfaces (cap linkers)
        #[arg(long)]
        expose_nodes: bool,

        /// Prefer linker-terminated surfaces (cap nodes)
        #[arg(long)]
        expose_linkers: bool,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let start_time = Instant::now();

    match cli.command {
        Commands::Generate { 
            input, output, h, k, l, 
            thickness, vacuum, offset, reconstruct, 
            with_mofid, mofid_work_dir,
            expose_nodes, expose_linkers,
        } => {
            println!("--- Crystal Surface Generator ---");

            // Hard validation
            if expose_nodes && expose_linkers {
                anyhow::bail!("--expose-nodes and --expose-linkers cannot be used together.");
            }

            // 1. Parsing
            println!("Reading structure from {:?}...", input);
            let mut crystal = parser::from_cif(&input)?;
            println!("-> Loaded {} atoms.", crystal.atoms.len());

            // 2. Molecule Analysis
            println!("Analyzing connectivity...");
            let finder = MoleculeFinder::new(2.0);
            let molecules = finder.find_molecules(&crystal)?;
            if !molecules.is_empty() {
                println!("-> Detected {} discrete molecules.", molecules.len());
            }

            // 3. MOFid
            if with_mofid {
                println!("-> MOFid integration enabled.");
            }

            println!(
                "Surface mode: {}",
                if expose_nodes {
                    "node-terminated"
                } else if expose_linkers {
                    "linker-terminated"
                } else {
                    "neutral (no semantic capping)"
                }
            );

            // 4. Execution
            println!("Generating ({} {} {}) slab...", h, k, l);

            let config = SurfaceConfig {
                miller_indices: [h, k, l],
                thickness,
                vacuum,
                offset,
                reconstruct,

                // MOFid integration
                input_cif_path: Some(input),
                enable_mofid: with_mofid,
                mofid_output_root: mofid_work_dir,

                // NOTE:
                // Your current SurfaceConfig DOES NOT expose:
                //   expose_nodes / expose_linkers
                //
                // So we DO NOT set them here.
                // The CLI flags are parsed and validated,
                // and you can wire them into SurfaceConfig later
                // when you add the fields in lib.rs.
            };

            let (slab, report) = generate_surface(&mut crystal, &molecules, &config)?;

            println!("\nSuccess!");
            println!("{}", report);

            println!("Writing output to {:?}...", output);
            writer::to_cif(&slab, &output)?;

            println!(
                "Done in {:.2?}",
                start_time.elapsed()
            );
        }
    }

    Ok(())
}