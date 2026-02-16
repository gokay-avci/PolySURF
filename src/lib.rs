// ============================================================================
// MODULE DECLARATIONS
// ============================================================================
pub mod core;
pub mod io;
pub mod math;
pub mod analysis;
pub mod synthesis;
pub mod chemistry; // NEW: Export the chemistry module

// ============================================================================
// RE-EXPORTS (Public API)
// ============================================================================
pub use crate::core::structure::{Atom, Crystal, Lattice, Molecule, ComponentType};
pub use crate::core::connectivity::MoleculeFinder;
pub use crate::io::{parser, writer};

pub use crate::synthesis::builder::SlabBuilder;
pub use crate::synthesis::population::SlabPopulator;
pub use crate::synthesis::ionic::{IonicReconstructor, ReconstructionMode};
pub use crate::analysis::topology::VoidCrawler;
pub use crate::chemistry::tagging::SemanticTagger; // NEW

use anyhow::{Result, Context};
use std::path::PathBuf;
use nalgebra::Vector3;

// ============================================================================
// HIGH-LEVEL INTERFACE
// ============================================================================

/// Configuration for the surface generation pipeline.
#[derive(Debug, Clone)]
pub struct SurfaceConfig {
    pub miller_indices: [i32; 3],
    pub thickness: f64,
    pub vacuum: f64,
    pub offset: Option<f64>,
    pub reconstruct: bool,

    // --- NEW: Semantic Configuration ---
    /// Path to the original input CIF (required for MOFid execution).
    pub input_cif_path: Option<PathBuf>,
    /// Enable MOFid decomposition?
    pub enable_mofid: bool,
    /// Directory for MOFid intermediate files (default: "mofid_work").
    pub mofid_output_root: Option<PathBuf>,
}

/// The Master Pipeline function.
pub fn generate_surface(
    crystal: &mut Crystal, // CHANGED: Mutable to allow tagging
    molecules: &[Molecule], 
    config: &SurfaceConfig
) -> Result<(Crystal, String)> {
    
    let mut report_buffer = String::new();

    // 0. SEMANTIC PHASE (The Integration Logic)
    if config.enable_mofid {
        if let Some(cif_path) = &config.input_cif_path {
            let root = config.mofid_output_root.clone().unwrap_or_else(|| PathBuf::from("mofid_work"));
            
            // Call the Semantic Engine
            // This runs the external binaries and generates geometry files
            let artifacts = mofid_rust::analyze_cif(cif_path, &root)
                .context("MOFid Analysis Failed")?;
            
            // Map Semantics to Geometry
            let tag_msg = SemanticTagger::tag_structure(crystal, &artifacts)?;
            report_buffer.push_str(&format!("--- Chemistry ---\n{}\n", tag_msg));
        } else {
            report_buffer.push_str("Warning: --with-mofid requested but input path lost. Skipping semantics.\n");
        }
    }

    // 1. MATH PHASE
    let builder = SlabBuilder::new(
        config.miller_indices[0], 
        config.miller_indices[1], 
        config.miller_indices[2], 
        config.thickness, 
        config.vacuum
    );
    let geometry = builder.compute_geometry(crystal)?;

    // 2. TOPOLOGY PHASE
    let offset = if let Some(u) = config.offset {
        u
    } else {
        let hkl = Vector3::new(
            config.miller_indices[0] as f64, 
            config.miller_indices[1] as f64, 
            config.miller_indices[2] as f64
        );
        let normal: Vector3<f64> = crystal.lattice.reciprocal_matrix * hkl;
        
        let crawler = VoidCrawler::new(crystal, &normal);
        crawler.find_safe_offsets().first().map(|c| c.offset_z).unwrap_or(0.0)
    };

    // 3. SYNTHESIS PHASE
    // Note: SlabPopulator reads `component_type` from atoms. 
    // If tagged in Phase 0, it can now make smarter decisions (future upgrade).
    let mut slab_atoms = SlabPopulator::populate(crystal, &geometry, molecules, offset)?;

    // 4. PHYSICS PHASE
    let slab_lattice = crate::core::structure::Lattice::new(geometry.basis)
        .map_err(|e| anyhow::anyhow!(e))?;
        
    let mode = if config.reconstruct { ReconstructionMode::DipoleCorrection } else { ReconstructionMode::None };
    
    let phys_report = IonicReconstructor::stabilize(&mut slab_atoms, &slab_lattice, mode)?;

    // 5. REPORT GENERATION
    let actual_material_thickness = geometry.n_layers as f64 * geometry.d_hkl;
    
    let verbose_output = format!(
        "{}\n--- Surface Generation Report ---\n\
         • Plane:           ({} {} {})\n\
         • Interplanar Spacing: {:.4} Å\n\
         • Quantization:    Requested {:.2} Å → {} Full Layers\n\
         • Final Thickness: {:.4} Å (Material) + {:.2} Å (Vacuum)\n\
         • Cut Offset:      {:.4} Å\n\
         • Physics:         {}", 
        report_buffer,
        config.miller_indices[0], config.miller_indices[1], config.miller_indices[2],
        geometry.d_hkl,
        config.thickness, geometry.n_layers,
        actual_material_thickness, config.vacuum,
        offset,
        phys_report
    );

    Ok((
        Crystal { lattice: slab_lattice, atoms: slab_atoms },
        verbose_output
    ))
}