use crate::core::structure::{Atom, Crystal, ComponentType, Lattice};
use mofid_rust::types::MofArtifacts;
use std::fs;
use std::path::Path;
use glob::glob;
use nalgebra::Vector3;
use anyhow::{Result, Context, anyhow};
use std::collections::HashMap;
use crate::io::parser;

/// Engine for mapping semantic identities (Node/Linker) onto geometric atoms.
pub struct SemanticTagger;

impl SemanticTagger {
    /// Reads MOFid outputs (XYZ fragments) and updates the Crystal's atom tags.
    ///
    /// # Arguments
    /// * `crystal` - The mutable bulk crystal structure.
    /// * `artifacts` - The file manifest returned by `mofid_rust::analyze_cif`.
    ///
    /// # Returns
    /// A summary report string.
    pub fn tag_structure(crystal: &mut Crystal, artifacts: &MofArtifacts) -> Result<String> {
        let mut node_count = 0;
        let mut linker_count = 0;

        // 1. Build Spatial Index for Fast Lookup
        // We create a map from (grid_i, grid_j, grid_k) -> Vec<atom_index>
        // This reduces matching from O(N*M) to O(M).
        let spatial_index = SpatialIndex::new(crystal, 2.0); // 2.0 Angstrom cell size

        // 2. Load and Tag Nodes (SBUs)
        //
        // Legacy outputs produce individual XYZ files under `Nodes/`.
        // Newer AllNode outputs produce a `nodes.cif` inside the directory.
        // Some AllNode outputs may contain additional `.cif` files (e.g. branch_points.cif)
        // which are not of interest and may be empty. We restrict the fallback to
        // explicitly named files and gracefully handle parse failures.
        let node_xyz_pattern = artifacts.nodes_dir.join("*.xyz");
        let node_xyz_paths: Vec<_> = glob(node_xyz_pattern.to_str().unwrap())
            .context("Invalid glob pattern for nodes")?
            .collect();
        if !node_xyz_paths.is_empty() {
            for entry in node_xyz_paths {
                if let Ok(path) = entry {
                    let coords = match Self::parse_fragment_xyz(&path) {
                        Ok(c) => c,
                        Err(_) => continue,
                    };
                    if Self::apply_tag(crystal, &spatial_index, &coords, ComponentType::MetalNode) {
                        node_count += 1;
                    }
                }
            }
        } else {
            // Fallback: parse nodes.cif if present
            let cif_path = artifacts.nodes_dir.join("nodes.cif");
            if cif_path.exists() {
                if let Ok(coords) = Self::parse_fragment_cif(&cif_path) {
                    if Self::apply_tag(crystal, &spatial_index, &coords, ComponentType::MetalNode) {
                        node_count += 1;
                    }
                }
            }
        }

        // 3. Load and Tag Linkers (Edges)
        // As with nodes, we first look for XYZ fragments. If none are present,
        // we look for explicitly named CIFs. Recent MOFid versions may produce
        // an "edges.cif" rather than "linkers.cif"; we support both names.
        let linker_xyz_pattern = artifacts.linkers_dir.join("*.xyz");
        let linker_xyz_paths: Vec<_> = glob(linker_xyz_pattern.to_str().unwrap())
            .context("Invalid glob pattern for linkers")?
            .collect();
        if !linker_xyz_paths.is_empty() {
            for entry in linker_xyz_paths {
                if let Ok(path) = entry {
                    let coords = match Self::parse_fragment_xyz(&path) {
                        Ok(c) => c,
                        Err(_) => continue,
                    };
                    if Self::apply_tag(crystal, &spatial_index, &coords, ComponentType::OrganicLinker) {
                        linker_count += 1;
                    }
                }
            }
        } else {
            // Fallback: try edges.cif, then linkers.cif
            let candidate_names = ["edges.cif", "linkers.cif"];
            for name in &candidate_names {
                let cif_path = artifacts.linkers_dir.join(name);
                if cif_path.exists() {
                    if let Ok(coords) = Self::parse_fragment_cif(&cif_path) {
                        if Self::apply_tag(crystal, &spatial_index, &coords, ComponentType::OrganicLinker) {
                            linker_count += 1;
                        }
                    }
                    // Use the first existing candidate; do not continue to next
                    break;
                }
            }
        }

        // 4. Statistics Check
        let total = crystal.atoms.len();
        let tagged = crystal.atoms.iter().filter(|a| a.component_type != ComponentType::Unknown).count();
        
        Ok(format!(
            "Semantic Tagging Complete.\n\
             • Metal Nodes Found: {}\n\
             • Linkers Found:     {}\n\
             • Coverage:          {}/{} atoms ({:.1}%) tagged.", 
            node_count, linker_count, tagged, total, (tagged as f64 / total as f64) * 100.0
        ))
    }

    /// Helper to assign a tag to bulk atoms that match fragment coordinates.
    fn apply_tag(
        crystal: &mut Crystal, 
        index: &SpatialIndex, 
        fragment_coords: &[Vector3<f64>], 
        tag: ComponentType
    ) -> bool {
        let tolerance_sq = 0.5_f64.powi(2); // 0.5 Angstrom squared tolerance (loose to handle rounding)
        let mut matched_any = false;

        for frag_pos_cart in fragment_coords {
            // Convert Fragment (Cartesian) -> Fractional for PBC check
            let frag_pos_frac = crystal.lattice.to_fractional(frag_pos_cart);

            // Query the spatial index for nearby bulk atoms
            // We need to check the exact grid cell and neighbors because of edge cases
            let candidates = index.query(frag_pos_cart);

            for &atom_idx in &candidates {
                let atom = &mut crystal.atoms[atom_idx];

                // PRIORITY RULE: Never overwrite a MetalNode with a Linker.
                // Metals are the anchors.
                if atom.component_type == ComponentType::MetalNode && tag == ComponentType::OrganicLinker {
                    continue;
                }

                // Distance Check with PBC
                let dist_vec = crystal.lattice.get_shortest_distance_vector(
                    &atom.fractional_coords,
                    &frag_pos_frac
                );

                if dist_vec.norm_squared() < tolerance_sq {
                    atom.component_type = tag;
                    matched_any = true;
                    // We keep searching in case of overlaps, but practically we found our atom.
                }
            }
        }
        matched_any
    }

    /// Minimal, robust XYZ parser.
    fn parse_fragment_xyz(path: &Path) -> Result<Vec<Vector3<f64>>> {
        let content = fs::read_to_string(path).context("Failed to read XYZ fragment")?;
        let mut coords = Vec::with_capacity(content.lines().count());
        
        // Skip header (Atom Count + Comment)
        // Standard XYZ has 2 header lines.
        for line in content.lines().skip(2) {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 4 {
                // Format: Element X Y Z
                // We only care about geometry here
                let x: f64 = parts[1].parse().unwrap_or(0.0);
                let y: f64 = parts[2].parse().unwrap_or(0.0);
                let z: f64 = parts[3].parse().unwrap_or(0.0);
                coords.push(Vector3::new(x, y, z));
            }
        }
        Ok(coords)
    }

    /// Parse a CIF fragment file and return Cartesian coordinates.
    ///
    /// Uses the internal CIF parser to load the structure and converts
    /// fractional coordinates to Cartesian coordinates using the fragment's
    /// lattice. Only positional data are extracted; chemical identity is ignored.
    fn parse_fragment_cif(path: &Path) -> Result<Vec<Vector3<f64>>> {
        // Use the crate's CIF parser to build a temporary Crystal
        let frag_crystal = parser::from_cif(path)
            .with_context(|| format!("Failed to parse CIF fragment at {:?}", path))?;
        let mut coords = Vec::new();
        for atom in &frag_crystal.atoms {
            let cart = frag_crystal.lattice.to_cartesian(&atom.fractional_coords);
            coords.push(cart);
        }
        Ok(coords)
    }
}

// ============================================================================
// SPATIAL ACCELERATION
// ============================================================================

/// A simple spatial hash grid to accelerate atom lookups.
struct SpatialIndex {
    cell_size: f64,
    grid: HashMap<(i32, i32, i32), Vec<usize>>,
}

impl SpatialIndex {
    fn new(crystal: &Crystal, cell_size: f64) -> Self {
        let mut grid = HashMap::new();
        for (i, atom) in crystal.atoms.iter().enumerate() {
            let cart = crystal.lattice.to_cartesian(&atom.fractional_coords);
            let key = Self::get_key(&cart, cell_size);
            grid.entry(key).or_insert_with(Vec::new).push(i);
        }
        Self { cell_size, grid }
    }

    fn get_key(pos: &Vector3<f64>, cell_size: f64) -> (i32, i32, i32) {
        (
            (pos.x / cell_size).floor() as i32,
            (pos.y / cell_size).floor() as i32,
            (pos.z / cell_size).floor() as i32,
        )
    }

    /// Returns candidate atom indices near the query position.
    fn query(&self, pos: &Vector3<f64>) -> Vec<usize> {
        let (cx, cy, cz) = Self::get_key(pos, self.cell_size);
        let mut candidates = Vec::new();

        // Check 3x3x3 neighborhood to account for boundary conditions
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let key = (cx + dx, cy + dy, cz + dz);
                    if let Some(indices) = self.grid.get(&key) {
                        candidates.extend_from_slice(indices);
                    }
                }
            }
        }
        candidates
    }
}