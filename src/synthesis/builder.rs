use crate::core::structure::{Crystal, Lattice};
use crate::math::{integer_basis, lll};
use nalgebra::{Matrix3, Vector3};
use anyhow::{Result, anyhow};

// --- STRICT TYPE DEFINITIONS FOR CLARITY ---
type Cartesian3 = Vector3<f64>;

#[derive(Debug, Clone)]
pub struct SlabGeometry {
    pub basis: Matrix3<f64>,
    pub d_hkl: f64,
    pub n_layers: usize,
    pub vacuum_thickness: f64,
}

pub struct SlabBuilder {
    miller_indices: [i32; 3],
    target_thickness: f64,
    vacuum: f64,
}

impl SlabBuilder {
    pub fn new(h: i32, k: i32, l: i32, thickness: f64, vacuum: f64) -> Self {
        Self {
            miller_indices: [h, k, l],
            target_thickness: thickness,
            vacuum,
        }
    }

    pub fn compute_geometry(&self, crystal: &Crystal) -> Result<SlabGeometry> {
        let (h, k, l) = (self.miller_indices[0], self.miller_indices[1], self.miller_indices[2]);

        // 1. INTEGER PHASE
        let (u_raw, v_raw) = integer_basis::find_primitive_in_plane_basis(h, k, l)?;

        // 2. INTEGER REDUCTION
        let (u_int, v_int) = lll::reduce_2d_integer(u_raw, v_raw);

        // 3. CARTESIAN CONVERSION
        let u_cart: Cartesian3 = crystal.lattice.to_cartesian(&Vector3::new(u_int.x as f64, u_int.y as f64, u_int.z as f64));
        let v_cart: Cartesian3 = crystal.lattice.to_cartesian(&Vector3::new(v_int.x as f64, v_int.y as f64, v_int.z as f64));

        // 4. ASPECT RATIO CHECK
        let len_u = u_cart.norm();
        let len_v = v_cart.norm();
        let ratio = if len_u > len_v { len_u / len_v } else { len_v / len_u };
        
        if ratio > 5.0 {
            eprintln!("Warning: High Aspect Ratio ({:.1}) detected for surface ({}{}{}).", ratio, h,k,l);
        }

        // 5. STACKING
        let reciprocal_n = crystal.lattice.reciprocal_matrix * Vector3::new(h as f64, k as f64, l as f64);
        let g_norm = reciprocal_n.norm();
        if g_norm < 1e-9 { return Err(anyhow!("Invalid Miller indices.")); }
        let d_hkl = 1.0 / g_norm;

        let n_layers = (self.target_thickness / d_hkl).round().max(1.0) as usize;
        let slab_height = n_layers as f64 * d_hkl;

        let normal = reciprocal_n.normalize();
        let c_slab: Cartesian3 = normal * (slab_height + self.vacuum);

        let mut basis = Matrix3::from_columns(&[u_cart, v_cart, c_slab]);

        // 6. 3D REDUCTION (Optional but good)
        // We define a temp basis with a huge Z to prevent mixing C into A/B
        let temp_basis = Matrix3::from_columns(&[
            u_cart,
            v_cart,
            normal * 10000.0 
        ]);
        
        let reduced = lll::lll_reduce(temp_basis);
        basis.set_column(0, &reduced.column(0));
        basis.set_column(1, &reduced.column(1));

        Ok(SlabGeometry {
            basis,
            d_hkl,
            n_layers,
            vacuum_thickness: self.vacuum,
        })
    }
}