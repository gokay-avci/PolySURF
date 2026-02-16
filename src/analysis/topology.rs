use crate::core::structure::{Atom, Crystal};
use nalgebra::Vector3;
use std::cmp::Ordering;

/// Represents a potential slice plane through the crystal.
#[derive(Debug, Clone, Copy)]
pub struct SafeCut {
    /// The Z-coordinate (along the surface normal) where the cut should occur.
    pub offset_z: f64,
    /// The size of the empty gap (in Å) at this location.
    pub gap_size: f64,
    /// A heuristic score (0.0 to 1.0) indicating how "clean" this cut is.
    pub quality_score: f64,
}

/// A rigorous 1D density analyzer for determining optimal slab terminations.
pub struct VoidCrawler {
    /// List of (center_z, radius) for every atom, projected onto the normal.
    projections: Vec<(f64, f64)>,
    /// The repeat distance along the normal direction.
    periodicity: f64,
}

impl VoidCrawler {
    /// Initializes the crawler by projecting the bulk crystal onto the surface normal.
    ///
    /// # Arguments
    /// * `crystal` - The bulk structure.
    /// * `surface_normal` - The Cartesian vector normal to the surface plane.
    pub fn new(crystal: &Crystal, surface_normal: &Vector3<f64>) -> Self {
        let normal_normalized = surface_normal.normalize();

        // 1. Calculate Periodicity along the normal.
        // We need to know how often the bulk repeats in this direction to handle wrapping.
        // We project the three lattice vectors onto the normal. 
        // The periodicity is determined by the specific stacking vector, but for a general 
        // 1D density map, the projection of the c-axis (or the largest projection) 
        // usually defines the "repeat block" for the crawler's domain.
        let p_a = crystal.lattice.matrix.column(0).dot(&normal_normalized).abs();
        let p_b = crystal.lattice.matrix.column(1).dot(&normal_normalized).abs();
        let p_c = crystal.lattice.matrix.column(2).dot(&normal_normalized).abs();
        
        // We take the max projection as the domain size to be safe.
        // In a perfect projection, this equals d_hkl * N_layers_in_unit_cell.
        let periodicity = p_a.max(p_b).max(p_c);

        // 2. Project all atoms
        let mut projections = Vec::with_capacity(crystal.atoms.len() * 3);
        
        for atom in &crystal.atoms {
            let cart = crystal.lattice.to_cartesian(&atom.fractional_coords);
            let z = cart.dot(&normal_normalized);
            let r = Self::get_vdw_radius(&atom.element);
            
            // Map z into [0, periodicity)
            let z_mod = z.rem_euclid(periodicity);
            
            // To handle periodic boundary conditions (gaps wrapping around edges),
            // we replicate the spheres at -P, 0, and +P.
            projections.push((z_mod, r));
            projections.push((z_mod - periodicity, r));
            projections.push((z_mod + periodicity, r));
        }

        // Sort by Z center for the sweep-line algorithm
        projections.sort_by(|a, b| a.0.total_cmp(&b.0));

        Self {
            projections,
            periodicity,
        }
    }

    /// Finds the largest gaps in the atomic density.
    /// Returns a list of safe offsets sorted by gap size (best/largest gap first).
    pub fn find_safe_offsets(&self) -> Vec<SafeCut> {
        if self.projections.is_empty() {
            return vec![SafeCut { offset_z: 0.0, gap_size: 10.0, quality_score: 1.0 }];
        }

        // 1. Merge Intervals (The "Sweep Line" Algorithm)
        // Convert spheres into 1D exclusion zones [z-r, z+r]
        let mut intervals: Vec<(f64, f64)> = self.projections.iter()
            .map(|(z, r)| (z - r, z + r))
            .collect();

        // Sort by start point
        intervals.sort_by(|a, b| a.0.total_cmp(&b.0));

        let mut merged = Vec::new();
        if let Some(first) = intervals.first() {
            let mut current_start = first.0;
            let mut current_end = first.1;

            for next in intervals.iter().skip(1) {
                if next.0 < current_end {
                    // Overlap detected: extend the current block
                    current_end = current_end.max(next.1);
                } else {
                    // Gap detected: push current block and start a new one
                    merged.push((current_start, current_end));
                    current_start = next.0;
                    current_end = next.1;
                }
            }
            merged.push((current_start, current_end));
        }

        // 2. Find Gaps between merged intervals
        let mut cuts = Vec::new();
        
        for window in merged.windows(2) {
            let occupied_end = window[0].1;
            let occupied_start = window[1].0;
            
            let gap_size = occupied_start - occupied_end;
            
            // Filter: Positive gaps only
            if gap_size > 0.01 {
                let mid_point = occupied_end + gap_size / 2.0;
                
                // Only return cuts that fall effectively within the primary unit cell [0, P]
                if mid_point >= 0.0 && mid_point <= self.periodicity {
                     // Score logic: 3.0 Å is considered a "perfect" Van der Waals gap.
                     let score = (gap_size / 3.0).min(1.0);
                     
                     cuts.push(SafeCut {
                        offset_z: mid_point,
                        gap_size,
                        quality_score: score,
                    });
                }
            }
        }

        // Sort by gap size descending (Best first)
        cuts.sort_by(|a, b| b.gap_size.total_cmp(&a.gap_size));
        
        cuts
    }

    /// Returns the Van der Waals radius for a given element.
    /// Data Source: Alvarez, S. (2013). Dalton Trans., 42, 8617-8636.
    fn get_vdw_radius(element: &str) -> f64 {
        match element {
            // Period 1
            "H" => 1.20, "He" => 1.40,
            // Period 2
            "Li" => 1.82, "Be" => 1.53, "B" => 1.92, "C" => 1.70, 
            "N" => 1.55, "O" => 1.52, "F" => 1.47, "Ne" => 1.54,
            // Period 3
            "Na" => 2.27, "Mg" => 1.73, "Al" => 1.84, "Si" => 2.10, 
            "P" => 1.80, "S" => 1.80, "Cl" => 1.75, "Ar" => 1.88,
            // Period 4
            "K" => 2.75, "Ca" => 2.31, "Sc" => 2.11, "Ti" => 2.00, "V" => 2.00, "Cr" => 2.00,
            "Mn" => 2.00, "Fe" => 2.00, "Co" => 2.00, "Ni" => 1.63, "Cu" => 1.40, "Zn" => 1.39,
            "Ga" => 1.87, "Ge" => 2.11, "As" => 1.85, "Se" => 1.90, "Br" => 1.85, "Kr" => 2.02,
            // Period 5
            "Rb" => 3.03, "Sr" => 2.49, "Pd" => 1.63, "Ag" => 1.72, "Cd" => 1.58,
            "In" => 1.93, "Sn" => 2.17, "Sb" => 2.06, "Te" => 2.06, "I" => 1.98, "Xe" => 2.16,
            // Period 6
            "Cs" => 3.43, "Ba" => 2.68, "Pt" => 1.75, "Au" => 1.66, "Hg" => 1.55,
            "Tl" => 1.96, "Pb" => 2.02, "Bi" => 2.07, "Po" => 1.97, "At" => 2.02, "Rn" => 2.20,
            // Default fallback
            _ => 1.80, 
        }
    }
}