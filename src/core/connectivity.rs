use crate::core::structure::{Crystal, Molecule};
use petgraph::graph::{NodeIndex, UnGraph};
use petgraph::visit::Bfs;
use nalgebra::Vector3;
use std::collections::{HashMap, HashSet, VecDeque};
use anyhow::Result; // Defensive error handling

// ============================================================================
// GRAPH REPRESENTATION
// ============================================================================

/// Represents the connectivity of atoms within the crystal.
/// Used to identify distinct molecules or clusters.
pub struct GraphRepresentation {
    /// Undirected graph where nodes are atom indices and edges represent bonds.
    pub graph: UnGraph<usize, ()>,
}

impl GraphRepresentation {
    /// Builds the connectivity graph using the Minimum Image Convention (MIC).
    /// Atoms closer than `cutoff` are considered bonded.
    ///
    /// # Complexity
    /// O(N^2) currently. For systems > 10,000 atoms, a Cell List or KD-Tree should be used.
    pub fn from_crystal(crystal: &Crystal, cutoff: f64) -> Self {
        let num_atoms = crystal.atoms.len();
        let mut graph = UnGraph::<usize, ()>::with_capacity(num_atoms, num_atoms * 3);
        
        // Add all nodes first to maintain index mapping
        let node_indices: Vec<NodeIndex> = (0..num_atoms)
            .map(|i| graph.add_node(i))
            .collect();
            
        let cutoff_sq = cutoff.powi(2);

        // Iterate pairwise to find bonds
        for i in 0..num_atoms {
            for j in (i + 1)..num_atoms {
                // Calculate distance considering Periodic Boundary Conditions
                let dist_vec = crystal.lattice.get_shortest_distance_vector(
                    &crystal.atoms[i].fractional_coords,
                    &crystal.atoms[j].fractional_coords,
                );

                if dist_vec.norm_squared() < cutoff_sq {
                    graph.add_edge(node_indices[i], node_indices[j], ());
                }
            }
        }
        Self { graph }
    }

    /// Finds all connected components (subgraphs) in the graph.
    /// Returns a list of vectors, where each vector contains the atom indices of a molecule.
    pub fn find_connected_components(&self) -> Vec<Vec<usize>> {
        let mut visited = vec![false; self.graph.node_count()];
        let mut all_components = Vec::new();

        for i in 0..self.graph.node_count() {
            if !visited[i] {
                let mut component = Vec::new();
                let start_node = NodeIndex::new(i);
                
                let mut bfs = Bfs::new(&self.graph, start_node);
                while let Some(nx) = bfs.next(&self.graph) {
                    let atom_index = self.graph[nx]; // The payload of the node
                    if !visited[atom_index] {
                        visited[atom_index] = true;
                        component.push(atom_index);
                    }
                }
                
                // Sort for determinism
                component.sort_unstable();
                if !component.is_empty() {
                    all_components.push(component);
                }
            }
        }
        all_components
    }
}

// ============================================================================
// MOLECULE FINDER
// ============================================================================

/// Engine for detecting molecules in a periodic crystal.
pub struct MoleculeFinder {
    bond_cutoff: f64,
}

impl MoleculeFinder {
    pub fn new(cutoff: f64) -> Self {
        Self { bond_cutoff: cutoff }
    }

    /// Primary entry point: Finds molecules and returns them as robust `Molecule` objects.
    pub fn find_molecules(&self, crystal: &Crystal) -> Result<Vec<Molecule>> {
        let (molecules, _) = self.find_molecules_with_indices(crystal)?;
        Ok(molecules)
    }

    /// Advanced entry point: Returns molecules AND a set of atom indices that were assigned.
    /// Useful for debugging (finding "orphan" atoms).
    pub fn find_molecules_with_indices(
        &self,
        crystal: &Crystal,
    ) -> Result<(Vec<Molecule>, HashSet<usize>)> {
        if crystal.atoms.is_empty() {
            return Ok((Vec::new(), HashSet::new()));
        }

        // 1. Build Graph
        let crystal_graph = GraphRepresentation::from_crystal(crystal, self.bond_cutoff);
        
        // 2. Find Components (Indices only)
        let components = crystal_graph.find_connected_components();

        let mut molecules = Vec::with_capacity(components.len());
        let mut assigned_indices = HashSet::new();
        let original_atoms = &crystal.atoms;

        for indices in components {
            if indices.is_empty() { continue; }
            
            // Mark indices as processed
            indices.iter().for_each(|&i| { assigned_indices.insert(i); });

            // 3. Reassemble the Molecule (Unwrap PBC)
            // We cannot just use the raw coordinates because they might be wrapped.
            // We start at one atom and "walk" the graph, adding the relative vectors.
            
            let mut reassembled_atoms: HashMap<usize, Vector3<f64>> = HashMap::with_capacity(indices.len());
            let mut queue = VecDeque::new();
            
            let start_index = indices[0];
            // Fix the first atom at its Cartesian position
            reassembled_atoms.insert(
                start_index,
                crystal.lattice.to_cartesian(&original_atoms[start_index].fractional_coords),
            );
            queue.push_back(start_index);

            // Standard BFS to unwrap coordinates
            while let Some(current_idx) = queue.pop_front() {
                let current_pos = reassembled_atoms[&current_idx]; // Should always exist
                
                for neighbor in crystal_graph.graph.neighbors(NodeIndex::new(current_idx)) {
                    let neighbor_idx = *crystal_graph.graph.node_weight(neighbor).unwrap();
                    
                    // Only process if part of this component (graph ensures this) 
                    // and not yet visited in this reconstruction pass
                    if !reassembled_atoms.contains_key(&neighbor_idx) {
                        // Crucial Step: Get the shortest vector across PBC
                        let shortest_vec = crystal.lattice.get_shortest_distance_vector(
                            &original_atoms[current_idx].fractional_coords,
                            &original_atoms[neighbor_idx].fractional_coords,
                        );
                        
                        // Construct neighbor position relative to current
                        let neighbor_pos = current_pos + shortest_vec;
                        
                        reassembled_atoms.insert(neighbor_idx, neighbor_pos);
                        queue.push_back(neighbor_idx);
                    }
                }
            }

            // 4. Calculate Center of Mass (COM)
            let mut com = Vector3::zeros();
            for pos in reassembled_atoms.values() {
                com += pos;
            }
            com /= reassembled_atoms.len() as f64;

            // 5. Shift Logic (Visualization Polish)
            // The molecule might be reconstructed far outside the unit cell.
            // We shift the COM to be inside the [0, 0, 0] unit cell for neatness.
            let frac_com = crystal.lattice.to_fractional(&com);
            let shift_frac = -frac_com.map(|x| x.floor()); // Integer shift to bring [0..1]
            let shift_cart = crystal.lattice.to_cartesian(&shift_frac);

            let mut final_com = Vector3::zeros();
            let mut atoms_cart = Vec::with_capacity(indices.len());

            // Apply final shift
            for &idx in &indices {
                // If the graph was disconnected, reassembled_atoms might miss an index 
                // (shouldn't happen with correct BFS).
                if let Some(&pos) = reassembled_atoms.get(&idx) {
                    let final_pos = pos + shift_cart;
                    final_com += final_pos;
                    atoms_cart.push((original_atoms[idx].element.clone(), final_pos));
                }
            }
            final_com /= atoms_cart.len() as f64;

            molecules.push(Molecule {
                atoms: atoms_cart,
                center_of_mass: final_com,
            });
        }

        Ok((molecules, assigned_indices))
    }
}