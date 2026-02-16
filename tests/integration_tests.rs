use crystal_surface_generator::{parser, generate_surface, SurfaceConfig, MoleculeFinder};
use std::path::PathBuf;

#[test]
fn test_generate_surface_all_samples() {
    let inputs = vec![
        "A_sample_inputs/1.cif",
        "A_sample_inputs/2.cif",
        "A_sample_inputs/3.cif",
        "A_sample_inputs/4.cif",
        "A_sample_inputs/5.cif",
        "A_sample_inputs/manual.cif",
        "A_sample_inputs/mgo.cif",
    ];

    let root_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

    for input_filename in inputs {
        let input_path = root_dir.join(input_filename);
        println!("Testing: {:?}", input_path);

        // Ensure file exists
        assert!(input_path.exists(), "Test file not found: {:?}", input_path);

        // 1. Parsing
        let mut crystal = parser::from_cif(&input_path).expect("Failed to parse CIF");

        // 2. Molecule Analysis
        // Using a bond tolerance of 2.0 Angstroms
        let finder = MoleculeFinder::new(2.0);
        let molecules = finder.find_molecules(&crystal).expect("Failed to find molecules");

        // 3. Configuration
        // Testing a standard (1 0 0) surface generation
        let config = SurfaceConfig {
            miller_indices: [1, 0, 0],
            thickness: 15.0,
            vacuum: 15.0,
            offset: None,
            reconstruct: false,
            input_cif_path: Some(input_path.clone()),
            enable_mofid: false,
            mofid_output_root: None,
        };

        // 4. Execution
        let result = generate_surface(&mut crystal, &molecules, &config);

        if let Err(e) = &result {
            println!("Error generating surface for {:?}: {}", input_path, e);
        }
        assert!(result.is_ok(), "Failed to generate surface for {:?}", input_path);

        let (slab, report) = result.unwrap();

        // Basic assertions on output
        assert!(slab.atoms.len() > 0, "Generated slab has no atoms");
        assert!(!report.is_empty(), "Report is empty");

        println!("Successfully generated surface for {:?}. Atoms: {}", input_filename, slab.atoms.len());
    }
}
