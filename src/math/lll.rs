use nalgebra::{Matrix3, Vector3};

/// Performs LLL Lattice Reduction on a 3D basis (Floating Point).
pub fn lll_reduce(basis: Matrix3<f64>) -> Matrix3<f64> {
    let delta = 0.75; 
    let mut b = basis;
    let n = 3;
    let mut k = 1;

    while k < n {
        let mut b_star = b; 
        for i in 0..k+1 {
            for j in 0..i {
                let mu = b.column(i).dot(&b_star.column(j)) / b_star.column(j).dot(&b_star.column(j));
                let proj = b_star.column(j) * mu;
                let mut col_i = b_star.column(i).into_owned();
                col_i -= proj;
                b_star.set_column(i, &col_i);
            }
        }

        for j in (0..k).rev() {
            let mu = b.column(k).dot(&b_star.column(j)) / b_star.column(j).dot(&b_star.column(j));
            if mu.abs() > 0.5 {
                let coefficient = mu.round();
                let sub = b.column(j) * coefficient;
                let mut col_k = b.column(k).into_owned();
                col_k -= sub;
                b.set_column(k, &col_k);
            }
        }

        let mu_k_km1 = b.column(k).dot(&b_star.column(k - 1)) / b_star.column(k - 1).dot(&b_star.column(k - 1));
        
        // We calculate norms for the condition check
        let b_star_km1_norm = b_star.column(k - 1).norm_squared();
        
        let condition_met = b_star.column(k).norm_squared() >= (delta - mu_k_km1.powi(2)) * b_star_km1_norm;

        if !condition_met {
            b.swap_columns(k, k - 1);
            k = 1.max(k - 1); 
        } else {
            k += 1;
        }
    }
    b
}

/// Performs Lagrange-Gauss reduction on two 3D integer vectors.
/// FIXED: Uses .dot() instead of .norm_squared() for integer types.
pub fn reduce_2d_integer(mut u: Vector3<i32>, mut v: Vector3<i32>) -> (Vector3<i32>, Vector3<i32>) {
    // Fix: u.norm_squared() -> u.dot(&u)
    if u.dot(&u) > v.dot(&v) {
        std::mem::swap(&mut u, &mut v);
    }

    loop {
        let dot = u.dot(&v);
        // Fix: u.norm_squared() -> u.dot(&u)
        let norm_sq = u.dot(&u); 
        
        if norm_sq == 0 { return (u, v); } 

        let mu = (dot as f64 / norm_sq as f64).round() as i32;

        if mu == 0 {
            return (u, v);
        }

        let v_new = v - u * mu;

        // Fix: v.norm_squared() -> v.dot(&v)
        if v_new.dot(&v_new) >= v.dot(&v) {
            return (u, v);
        }

        v = v_new;

        // Fix: u.norm_squared() -> u.dot(&u)
        if u.dot(&u) > v.dot(&v) {
            std::mem::swap(&mut u, &mut v);
        }
    }
}