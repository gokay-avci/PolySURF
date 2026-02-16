use crate::core::structure::{Atom, Crystal, Molecule, ComponentType};
use crate::synthesis::builder::SlabGeometry;
use nalgebra::Vector3;
use anyhow::{Result, anyhow};

pub struct SlabPopulator;

impl SlabPopulator {
    pub fn populate(
        crystal: &Crystal, 
        geometry: &SlabGeometry, 
        molecules: &[Molecule], 
        offset_z: f64
    ) -> Result<Vec<Atom>> {
        
        // Store interim atoms along with their semantic component type. The third
        // entry preserves the ComponentType from the bulk structure or assigns
        // Unknown for molecules. This allows downstream capping logic to
        // distinguish metal nodes from organic linkers when reconstructing the
        // slab surfaces.
        let mut final_atoms = Vec::new();
        
        // 1. Calculate Miller Plane Projection
        // We want to work in "Layer Index Space" rather than Angstroms to avoid high-index precision errors.
        // A point r in bulk fractional coords (u,v,w) lies on plane P = h*u + k*v + l*w.
        // If P is an integer, it's on a lattice plane.
        // The spacing between P and P+1 is d_hkl.
        
        // Recover the Miller indices from the geometry normal? 
        // We implicitly have them via d_hkl and the normal, but let's use the explicit d_hkl logic.
        
        let slab_c = geometry.basis.column(2);
        let slab_normal = slab_c.normalize();
        
        // Convert offset_z (Angstroms) to Layer Index (Float)
        // offset_idx = 3.5 means we start cutting halfway between layer 3 and 4.
        let offset_idx = offset_z / geometry.d_hkl;
        let n_layers_float = geometry.n_layers as f64;
        
        // Epsilon is now in "Layer Units". 1e-3 layers is very safe.
        let epsilon = 1e-3;
        let min_idx = offset_idx - epsilon;
        let max_idx = offset_idx + n_layers_float - epsilon;

        // 2. Dynamic Bounding (High Index Safe)
        // Project bulk vectors onto normal to determine sweep range
        let proj_a = crystal.lattice.matrix.column(0).dot(&slab_normal).abs();
        let proj_b = crystal.lattice.matrix.column(1).dot(&slab_normal).abs();
        let proj_c = crystal.lattice.matrix.column(2).dot(&slab_normal).abs();
        
        // This is the height of one unit cell in Angstroms
        let cell_height_ang = proj_a.max(proj_b).max(proj_c);
        
        if cell_height_ang < 1e-9 { return Err(anyhow!("Degenerate unit cell (zero height).")); }

        // Convert total slab height to bulk cell repeats
        let total_slab_height = geometry.n_layers as f64 * geometry.d_hkl;
        
        // Multiplier 2.5 is sufficient because we check Com/Bounds later
        let repeats = (total_slab_height / cell_height_ang).ceil() as i32 + 3;

        let range_iter = (-repeats..=repeats).flat_map(|i| 
            (-repeats..=repeats).flat_map(move |j| 
                (-repeats..=repeats).map(move |k| Vector3::new(i as f64, j as f64, k as f64))
            )
        );

        if !molecules.is_empty() {
            // --- Molecular Mode ---
            for cell_shift_frac in range_iter {
                let cell_shift_cart = crystal.lattice.to_cartesian(&cell_shift_frac);
                
                for mol in molecules {
                    let shifted_com = mol.center_of_mass + cell_shift_cart;
                    // Project COM onto Normal
                    let z_ang = shifted_com.dot(&slab_normal);

                    // Convert to Layer Index Space
                    let layer_val = z_ang / geometry.d_hkl;

                    if layer_val >= min_idx && layer_val < max_idx {
                        for (element, rel_pos) in &mol.atoms {
                             let final_pos = rel_pos + cell_shift_cart;
                             // Molecule atoms lack semantic tagging; default to Unknown
                             final_atoms.push((element.clone(), final_pos, ComponentType::Unknown));
                        }
                    }
                }
            }
        } else {
            // --- Atomic Mode (High Precision) ---
            for cell_shift_frac in range_iter {
                let cell_shift_cart = crystal.lattice.to_cartesian(&cell_shift_frac);
                
                for atom in &crystal.atoms {
                    let pos_cart = crystal.lattice.to_cartesian(&atom.fractional_coords) + cell_shift_cart;
                    let z_ang = pos_cart.dot(&slab_normal);
                    
                    // Layer Index Check (The Fix for High Indices)
                    let layer_val = z_ang / geometry.d_hkl;

                    if layer_val >= min_idx && layer_val < max_idx {
                        // Preserve the semantic tag from the bulk atom
                        final_atoms.push((atom.element.clone(), pos_cart, atom.component_type));
                    }
                }
            }
        }

        if final_atoms.is_empty() {
            return Err(anyhow!("Generated slab is empty. Thickness might be too small for the selected plane."));
        }

        // 3. Post-Processing: Center the Slab
        // Find material bounds
        let mut min_z = f64::INFINITY;
        let mut max_z = f64::NEG_INFINITY;

        for (_, pos, _) in &final_atoms {
            let z = pos.dot(&slab_normal);
            if z < min_z { min_z = z; }
            if z > max_z { max_z = z; }
        }
        
        let current_material_thickness = max_z - min_z;
        let total_box_height = slab_c.norm(); 
        
        // Target: Center the material in the box
        let target_z_start = (total_box_height - current_material_thickness) / 2.0;
        let shift_val = target_z_start - min_z;
        
        let shift_vec = slab_normal * shift_val;
        let slab_basis_inv = geometry.basis.try_inverse().ok_or(anyhow!("Slab basis singular"))?;

        let result_atoms: Vec<Atom> = final_atoms.into_iter().map(|(el, pos, comp_type)| {
            let shifted_cart = pos + shift_vec;
            let fractional = slab_basis_inv * shifted_cart;
            Atom {
                element: el,
                fractional_coords: fractional,
                // Preserve the component type for semantic aware capping
                component_type: comp_type,
            }
        }).collect();

        Ok(result_atoms)
    }
}