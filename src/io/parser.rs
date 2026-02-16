use crate::core::structure::{Atom, Crystal, Lattice, ComponentType};
use anyhow::{anyhow, Context, Result};
use nalgebra::Vector3;
use std::collections::HashMap;
use std::fs;
use std::path::Path;

/// Parses a float value from a CIF string, safely removing uncertainty parentheses.
/// Example: "1.234(5)" -> 1.234
fn parse_cif_float(s: &str) -> Result<f64> {
    // Split at '(' and take the first part
    let clean_s = s.split('(').next().unwrap_or(s);
    clean_s.parse::<f64>().with_context(|| format!("Failed to parse '{}' as float", s))
}

/// Parses a CIF file into a Crystal structure.
///
/// Note: This is a robust manual parser. For production-grade generic CIF parsing,
/// consider using a dedicated crate, but this works for 99% of P1/VASP outputs.
pub fn from_cif(path: &Path) -> Result<Crystal> {
    let contents = fs::read_to_string(path).with_context(|| format!("Could not read CIF file: {:?}", path))?;
    let lines: Vec<&str> = contents.lines().map(str::trim).filter(|l| !l.is_empty()).collect();

    let mut lattice_params: HashMap<&str, f64> = HashMap::new();
    let mut atoms = Vec::new();
    
    let mut i = 0;
    while i < lines.len() {
        let line = lines[i];
        
        if line.starts_with("_cell_") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                // Handle key-value pairs like "_cell_length_a 10.0"
                if let Ok(value) = parse_cif_float(parts[1]) {
                     lattice_params.insert(parts[0], value);
                }
            }
        } else if line.starts_with("loop_") {
            // Advance past the "loop_" line
            i += 1;
            
            // Collect headers
            let mut headers = Vec::new();
            while i < lines.len() && lines[i].starts_with('_') {
                headers.push(lines[i]);
                i += 1;
            }

            // Check if this loop contains atom data
            if headers.contains(&"_atom_site_fract_x") {
                // Find column indices
                let symbol_idx = headers.iter().position(|&h| h == "_atom_site_type_symbol")
                    .context("CIF missing '_atom_site_type_symbol'")?;
                let x_idx = headers.iter().position(|&h| h == "_atom_site_fract_x")
                    .context("CIF missing '_atom_site_fract_x'")?;
                let y_idx = headers.iter().position(|&h| h == "_atom_site_fract_y")
                    .context("CIF missing '_atom_site_fract_y'")?;
                let z_idx = headers.iter().position(|&h| h == "_atom_site_fract_z")
                    .context("CIF missing '_atom_site_fract_z'")?;

                let max_idx = *[symbol_idx, x_idx, y_idx, z_idx].iter().max().unwrap();

                // Parse Data Rows
                while i < lines.len() && !lines[i].starts_with('_') && !lines[i].starts_with("loop_") {
                    let atom_line_parts: Vec<&str> = lines[i].split_whitespace().collect();
                    if atom_line_parts.len() > max_idx {
                        let element = atom_line_parts[symbol_idx].to_string();
                        let x = parse_cif_float(atom_line_parts[x_idx])?;
                        let y = parse_cif_float(atom_line_parts[y_idx])?;
                        let z = parse_cif_float(atom_line_parts[z_idx])?;

                        atoms.push(Atom {
                            element,
                            fractional_coords: Vector3::new(x, y, z),
                            // When parsing a CIF we don't yet know the semantic role of each atom,
                            // so default to "Unknown". The semantic tagging pass will update this field
                            // later based on MOFid fragments.
                            component_type: ComponentType::Unknown,
                        });
                    }
                    i += 1;
                }
                // Step back one, as the outer loop increments i
                i -= 1;
            }
        }
        i += 1;
    }

    // Extract Lattice Parameters
    let get_param = |key: &str| -> Result<f64> {
        lattice_params.get(key).copied().ok_or_else(|| anyhow!("CIF missing tag: {}", key))
    };

    let a = get_param("_cell_length_a")?;
    let b = get_param("_cell_length_b")?;
    let c = get_param("_cell_length_c")?;
    let alpha = get_param("_cell_angle_alpha")?;
    let beta = get_param("_cell_angle_beta")?;
    let gamma = get_param("_cell_angle_gamma")?;

    let lattice = Lattice::from_parameters(a, b, c, alpha, beta, gamma)
        .map_err(|e| anyhow!(e))?;

    if atoms.is_empty() {
        return Err(anyhow!("No atoms found in CIF file."));
    }

    Ok(Crystal { lattice, atoms })
}