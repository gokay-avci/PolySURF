use nalgebra::Vector3;
use anyhow::{Result, anyhow};

/// Calculates the Greatest Common Divisor (Euclidean algorithm).
fn gcd(a: i32, b: i32) -> i32 {
    let mut a = a.abs();
    let mut b = b.abs();
    while b != 0 {
        let temp = b;
        b = a % b;
        a = temp;
    }
    a
}

/// Finds two integer vectors (u, v) that span the plane orthogonal to (h, k, l).
/// Guaranteed to be a primitive basis (smallest possible integer area).
pub fn find_primitive_in_plane_basis(h: i32, k: i32, l: i32) -> Result<(Vector3<i32>, Vector3<i32>)> {
    if h == 0 && k == 0 && l == 0 {
        return Err(anyhow!("Miller indices cannot be (0,0,0)."));
    }

    // 1. Construct the normal vector n
    let n = Vector3::new(h, k, l);

    // 2. Find a "trial" vector t that is NOT collinear with n.
    // We try basis vectors x, y, z. One must be non-collinear.
    // If n is along Z, we pick Y.
    let mut t = Vector3::new(0, 0, 1);
    
    // FIX: Use strict equality check against zero vector for integers
    // (norm_squared is not implemented for i32 types in nalgebra)
    if n.cross(&t) == Vector3::zeros() {
        t = Vector3::new(0, 1, 0);
    }

    // 3. First basis vector u = (n x t) / gcd
    // This gives an integer vector orthogonal to n.
    let u_raw = n.cross(&t);
    let common_u = gcd(gcd(u_raw.x, u_raw.y), u_raw.z);
    let u = u_raw / common_u;

    // 4. Second basis vector v = (n x u) / gcd
    // This guarantees orthogonality between n and v, and ensures linear independence from u.
    let v_raw = n.cross(&u);
    let common_v = gcd(gcd(v_raw.x, v_raw.y), v_raw.z);
    let v = v_raw / common_v;

    Ok((u, v))
}

pub fn find_stacking_vector(h: i32, k: i32, l: i32) -> Vector3<i32> {
    let n = Vector3::new(h, k, l);
    let target = gcd(gcd(h, k), l);
    
    for x in -10..=10 {
        for y in -10..=10 {
            for z in -10..=10 {
                let w = Vector3::new(x, y, z);
                if n.dot(&w) == target {
                    return w;
                }
            }
        }
    }
    
    if h != 0 { return Vector3::new(1,0,0); }
    if k != 0 { return Vector3::new(0,1,0); }
    Vector3::new(0,0,1)
}