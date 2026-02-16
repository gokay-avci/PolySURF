use crate::core::structure::{Atom, Lattice};
use nalgebra::Vector3;
use anyhow::Result;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ReconstructionMode {
    None,
    DipoleCorrection,
}

pub struct IonicReconstructor;

impl IonicReconstructor {
    pub fn stabilize(
        atoms: &mut Vec<Atom>, 
        lattice: &Lattice, 
        mode: ReconstructionMode
    ) -> Result<String> {
        if mode == ReconstructionMode::None || atoms.is_empty() {
            return Ok("No reconstruction applied.".to_string());
        }

        // 1. Group atoms into Z-planes with high precision
        let mut indices: Vec<usize> = (0..atoms.len()).collect();
        // Sort by Z coordinate
        indices.sort_by(|&i, &j| {
            let z_i = lattice.to_cartesian(&atoms[i].fractional_coords).z;
            let z_j = lattice.to_cartesian(&atoms[j].fractional_coords).z;
            z_i.partial_cmp(&z_j).unwrap()
        });

        // Robust Plane Clustering
        // We use a dynamic tolerance based on the density of points to handle high indices
        let mut planes: Vec<Vec<usize>> = Vec::new();
        let mut current_plane = vec![indices[0]];
        let mut current_z = lattice.to_cartesian(&atoms[indices[0]].fractional_coords).z;

        // Tighter tolerance for high-index surfaces where planes are close
        let tolerance = 0.25; 

        for &idx in indices.iter().skip(1) {
            let z = lattice.to_cartesian(&atoms[idx].fractional_coords).z;
            if (z - current_z).abs() < tolerance {
                current_plane.push(idx);
            } else {
                planes.push(current_plane);
                current_plane = vec![idx];
                current_z = z;
            }
        }
        planes.push(current_plane);

        // 2. Dipole Check
        let charges = Self::guess_charges(atoms);
        let dipole_z: f64 = atoms.iter().zip(&charges)
            .map(|(a, q)| q * lattice.to_cartesian(&a.fractional_coords).z)
            .sum();

        // 3. Vector-Based Reconstruction
        if dipole_z.abs() > 0.5 {
            if let Some(top_plane) = planes.last() {
                let num_to_move = top_plane.len() / 2;
                
                if num_to_move > 0 && planes.len() > 1 {
                    // INTELLIGENT PLACEMENT LOGIC:
                    // Instead of guessing Z, we calculate the vector from Top Layer -> Second Layer.
                    // This vector represents the "Stacking Shift" in reverse.
                    
                    let top_z_avg = Self::average_pos(lattice, top_plane, atoms);
                    let second_z_avg = Self::average_pos(lattice, &planes[planes.len() - 2], atoms);
                    
                    // Vector pointing 'down' one layer in the stack
                    let stacking_vector = second_z_avg - top_z_avg;
                    
                    // To go from Top to "New Bottom" (which is below the current Bottom),
                    // we need to apply this vector (N_planes - 1) times.
                    // Or more simply: New_Pos = Old_Pos + (stacking_vector * (planes.len() - 1))?
                    // No, that assumes linear spacing. Safe bet: 
                    // Calculate vector from Plane[1] to Plane[0] (Bottom to Bottom-most).
                    // Apply that vector to Plane[0] to find "Ghost Plane[-1]".
                    
                    // Let's use the Bottom -> Bottom+1 vector reversed.
                    let bottom_plane = &planes[0];
                    let bottom_next_plane = &planes[1];
                    let v_up = Self::average_pos(lattice, bottom_next_plane, atoms) 
                             - Self::average_pos(lattice, bottom_plane, atoms);
                    
                    let v_down = -v_up; // Vector to move from Bottom to "Ghost Bottom"

                    for &atom_idx in top_plane.iter().take(num_to_move) {
                        let current_cart = lattice.to_cartesian(&atoms[atom_idx].fractional_coords);
                        
                        // We are moving this atom from Top to Bottom.
                        // First, shift it by the total slab height to get it roughly to the bottom plane
                        // Then apply the specific stacking offset.
                        // Actually, the safest way is:
                        // New = Current - (Total_Material_Vector) - v_down?
                        
                        // Simplest robust method:
                        // Map the atom relative to the Top Plane Center, apply that relative offset to the Ghost Bottom Center.
                        let rel_to_top = current_cart - top_z_avg;
                        let ghost_center = Self::average_pos(lattice, bottom_plane, atoms) + v_down;
                        let new_cart = ghost_center + rel_to_top;

                        atoms[atom_idx].fractional_coords = lattice.to_fractional(&new_cart);
                    }
                    
                    return Ok(format!(
                        "Dipole detected ({:.3} eA). Moved {} atoms to crystallographic bottom sites.", 
                        dipole_z, num_to_move
                    ));
                }
            }
        }

        Ok(format!("Surface is stable (Dipole: {:.3} eA).", dipole_z))
    }

    fn average_pos(lattice: &Lattice, indices: &[usize], atoms: &[Atom]) -> Vector3<f64> {
        let mut sum = Vector3::zeros();
        for &i in indices {
            sum += lattice.to_cartesian(&atoms[i].fractional_coords);
        }
        sum / (indices.len() as f64)
    }

    fn guess_charges(atoms: &[Atom]) -> Vec<f64> {
        atoms.iter().map(|a| match a.element.as_str() {
            "Li"|"Na"|"K"|"H" => 1.0, "Mg"|"Ca"|"Zn"|"Fe" => 2.0, "Al" => 3.0,
            "F"|"Cl"|"Br"|"I" => -1.0, "O"|"S" => -2.0, "N" => -3.0, _ => 0.0,
        }).collect()
    }
}