use nalgebra::{Matrix3, Vector3};
use std::fmt;

// ============================================================================
// TRAITS
// ============================================================================

pub trait CifRepresentable {
    fn lattice(&self) -> &Lattice;
    fn atoms(&self) -> &Vec<Atom>;
}

// ============================================================================
// ENUMS & SEMANTICS
// ============================================================================

/// Defines the semantic role of an atom within a Framework (MOF/COF/Zeolite).
/// Used by the logic engine to determine where to slice.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ComponentType {
    /// Standard bulk atom (default).
    Unknown,
    /// Part of a Metal Node / Secondary Building Unit (SBU).
    /// *Action:* Avoid cutting if possible.
    MetalNode,
    /// Part of an Organic Linker / Strut.
    /// *Action:* Cut here, then saturate broken bonds.
    OrganicLinker,
    /// Solvent or Guest molecule in the pore.
    /// *Action:* Usually removed or ignored.
    Solvent,
    /// Gas molecule (e.g., CO2) adsorbed.
    Adsorbate,
}

// ============================================================================
// DATA STRUCTURES
// ============================================================================

#[derive(Debug, Clone)]
pub struct Lattice {
    pub matrix: Matrix3<f64>,
    pub reciprocal_matrix: Matrix3<f64>,
}

impl Lattice {
    pub fn new(matrix: Matrix3<f64>) -> Result<Self, &'static str> {
        if matrix.determinant().abs() < 1e-6 {
            return Err("Lattice has zero or near-zero volume.");
        }
        let reciprocal_matrix = matrix
            .try_inverse()
            .ok_or("Lattice is not invertible.")?
            .transpose();
        Ok(Self {
            matrix,
            reciprocal_matrix,
        })
    }

    pub fn from_parameters(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> Result<Self, &'static str> {
        let alpha_r = alpha.to_radians();
        let beta_r = beta.to_radians();
        let gamma_r = gamma.to_radians();

        let term = 1.0 - alpha_r.cos().powi(2) - beta_r.cos().powi(2) - gamma_r.cos().powi(2)
            + 2.0 * alpha_r.cos() * beta_r.cos() * gamma_r.cos();
        
        if term <= 0.0 { return Err("Invalid lattice angles."); }
        
        let v_factor = term.sqrt();
        let matrix = Matrix3::new(
            a, b * gamma_r.cos(), c * beta_r.cos(),
            0.0, b * gamma_r.sin(), c * (alpha_r.cos() - beta_r.cos() * gamma_r.cos()) / gamma_r.sin(),
            0.0, 0.0, a * b * c * v_factor / (a * b * gamma_r.sin()),
        );
        Self::new(matrix)
    }

    pub fn to_cartesian(&self, frac: &Vector3<f64>) -> Vector3<f64> { self.matrix * frac }
    pub fn to_fractional(&self, cart: &Vector3<f64>) -> Vector3<f64> { self.reciprocal_matrix.transpose() * cart }
    
    pub fn to_parameters(&self) -> (f64, f64, f64, f64, f64, f64) {
        let a = self.matrix.column(0).norm();
        let b = self.matrix.column(1).norm();
        let c = self.matrix.column(2).norm();
        let alpha = (self.matrix.column(1).dot(&self.matrix.column(2)) / (b * c)).acos().to_degrees();
        let beta = (self.matrix.column(0).dot(&self.matrix.column(2)) / (a * c)).acos().to_degrees();
        let gamma = (self.matrix.column(0).dot(&self.matrix.column(1)) / (a * b)).acos().to_degrees();
        (a, b, c, alpha, beta, gamma)
    }

    pub fn get_shortest_distance_vector(&self, f1: &Vector3<f64>, f2: &Vector3<f64>) -> Vector3<f64> {
        let mut d = f2 - f1;
        d.x -= d.x.round();
        d.y -= d.y.round();
        d.z -= d.z.round();
        self.to_cartesian(&d)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub element: String,
    pub fractional_coords: Vector3<f64>,
    /// Semantic tag identifying the atom's role in the framework.
    pub component_type: ComponentType,
}

#[derive(Debug, Clone)]
pub struct Crystal {
    pub lattice: Lattice,
    pub atoms: Vec<Atom>,
}

impl CifRepresentable for Crystal {
    fn lattice(&self) -> &Lattice { &self.lattice }
    fn atoms(&self) -> &Vec<Atom> { &self.atoms }
}

#[derive(Debug, Clone)]
pub struct Molecule {
    pub atoms: Vec<(String, Vector3<f64>)>,
    pub center_of_mass: Vector3<f64>,
}