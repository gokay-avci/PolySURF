use crate::types::{AnnotatedFragment, FragmentAtom};
use anyhow::{Result, Context};
use nalgebra::Vector3;
use std::fs;
use std::path::Path;
use glob::glob;

/// Scans the output directory of the SBU tool and parses the .xyz/.mol files
/// to retrieve the 3D geometry of nodes and linkers.
pub fn parse_sbu_output(output_dir: &Path) -> Result<(Vec<AnnotatedFragment>, Vec<AnnotatedFragment>)> {
    let mut nodes = Vec::new();
    let mut linkers = Vec::new();

    // 1. Parse Nodes
    let node_pattern = output_dir.join("Nodes").join("*.xyz"); // or .mol
    for entry in glob(node_pattern.to_str().unwrap())? {
        match entry {
            Ok(path) => {
                let atoms = parse_xyz(&path)?;
                // Read the SMILES if available (usually in filename or separate map)
                // For now, we use filename as placeholder ID
                let name = path.file_stem().unwrap().to_string_lossy().to_string();
                nodes.push(AnnotatedFragment { smiles: name, atoms });
            },
            Err(e) => eprintln!("Error reading node file: {:?}", e),
        }
    }

    // 2. Parse Linkers
    let linker_pattern = output_dir.join("Linkers").join("*.xyz");
    for entry in glob(linker_pattern.to_str().unwrap())? {
        match entry {
            Ok(path) => {
                let atoms = parse_xyz(&path)?;
                let name = path.file_stem().unwrap().to_string_lossy().to_string();
                linkers.push(AnnotatedFragment { smiles: name, atoms });
            },
            Err(e) => eprintln!("Error reading linker file: {:?}", e),
        }
    }

    Ok((nodes, linkers))
}

fn parse_xyz(path: &Path) -> Result<Vec<FragmentAtom>> {
    let content = fs::read_to_string(path).context("Failed to read XYZ")?;
    let mut atoms = Vec::new();
    
    // Simple XYZ parser
    // Skip first 2 lines (count and comment)
    for line in content.lines().skip(2) {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 4 {
            let element = parts[0].to_string();
            let x: f64 = parts[1].parse()?;
            let y: f64 = parts[2].parse()?;
            let z: f64 = parts[3].parse()?;
            
            atoms.push(FragmentAtom {
                element,
                position: Vector3::new(x, y, z),
            });
        }
    }
    Ok(atoms)
}