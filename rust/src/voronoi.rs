use qhull_sys as q;
use std::collections::BTreeSet;
use std::os::raw::{c_char, c_int, c_void};

struct Work {
    qh: *mut q::qhT,
    dim: c_int,
    n: c_int,
    points: *mut f64,
    cmd: *mut c_char,
    exitcode: c_int,
}

unsafe extern "C" fn run_qhull(data: *mut c_void) {
    let w = &mut *(data as *mut Work);
    w.exitcode = q::qh_new_qhull(
        w.qh,
        w.dim,
        w.n,
        w.points,
        0u32,
        w.cmd,
        std::ptr::null_mut(),
        q::qhull_sys__stderr(),
    );
}

/// Return cluster indices sharing a Voronoi facet with point 0.
///
/// Uses qhull's Delaunay mode with SciPy's flags (`d Qbb Qc Qz Q12`), matching
/// pymatgen's VoronoiNN(tol=0) neighbor set for OMS clusters.
pub fn center_neighbors(points: &[[f64; 3]]) -> Result<Vec<usize>, String> {
    if points.len() < 2 {
        return Ok(Vec::new());
    }
    unsafe {
        let mut qh: q::qhT = std::mem::zeroed();
        q::qh_zero(&mut qh, q::qhull_sys__stderr());

        let mut coords: Vec<f64> = points.iter().flat_map(|p| p.iter().copied()).collect();
        let mut cmd: Vec<u8> = b"qhull d Qbb Qc Qz Q12\0".to_vec();
        let mut work = Work {
            qh: &mut qh,
            dim: 3,
            n: points.len() as c_int,
            points: coords.as_mut_ptr(),
            cmd: cmd.as_mut_ptr() as *mut c_char,
            exitcode: 0,
        };

        let tryrc =
            q::qhull_sys__try_on_qh(&mut qh, Some(run_qhull), &mut work as *mut _ as *mut c_void);
        if tryrc != 0 || work.exitcode != 0 {
            q::qh_freeqhull(&mut qh, 1u32);
            return Err(format!(
                "qhull failed: tryrc={tryrc} exitcode={}",
                work.exitcode
            ));
        }

        let mut neighbors = BTreeSet::new();
        let mut facet = qh.facet_list;
        while !facet.is_null() && facet != qh.facet_tail {
            if (*facet).upperdelaunay() == 0 {
                let vertices = (*facet).vertices;
                if !vertices.is_null() {
                    let base = (*vertices).e.as_ptr();
                    let mut ids = Vec::new();
                    let mut i = 0isize;
                    while i < 4096 {
                        let ptr = (*base.offset(i)).p;
                        if ptr.is_null() {
                            break;
                        }
                        let vertex = ptr as *mut q::vertexT;
                        let id = q::qh_pointid(&mut qh, (*vertex).point);
                        if id >= 0 && (id as usize) < points.len() {
                            ids.push(id as usize);
                        }
                        i += 1;
                    }
                    if ids.contains(&0) {
                        for id in ids {
                            if id != 0 {
                                neighbors.insert(id);
                            }
                        }
                    }
                }
            }
            facet = (*facet).next;
        }
        q::qh_freeqhull(&mut qh, 1u32);
        Ok(neighbors.into_iter().collect())
    }
}

#[cfg(test)]
mod tests {
    use super::center_neighbors;

    #[test]
    fn tetrahedron_center_neighbors() {
        let points = vec![
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
        ];
        assert_eq!(center_neighbors(&points).unwrap(), vec![1, 2, 3, 4]);
    }
}
