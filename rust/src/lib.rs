use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

pub mod contacts;
pub mod eqeq;
pub mod graph;
pub mod oms;
pub mod pbc;
pub mod voronoi;

fn image_distance_to_object(
    py: Python<'_>,
    image_distance: pbc::ImageDistance,
) -> PyResult<PyObject> {
    let dict = PyDict::new_bound(py);
    dict.set_item("delta_frac", image_distance.delta_frac.to_vec())?;
    dict.set_item("image", image_distance.image.to_vec())?;
    dict.set_item("distance", image_distance.distance)?;
    Ok(dict.into_py(py))
}

fn short_contact_to_object(py: Python<'_>, contact: &contacts::ShortContact) -> PyResult<PyObject> {
    let dict = PyDict::new_bound(py);
    dict.set_item("i", contact.i)?;
    dict.set_item("j", contact.j)?;
    dict.set_item("image_j", contact.image_j.to_vec())?;
    dict.set_item("distance", contact.distance)?;
    dict.set_item("cutoff", contact.cutoff)?;
    Ok(dict.into_py(py))
}

fn neighbor_candidate_to_object(
    py: Python<'_>,
    candidate: &contacts::NeighborCandidate,
) -> PyResult<PyObject> {
    let dict = PyDict::new_bound(py);
    dict.set_item("i", candidate.i)?;
    dict.set_item("j", candidate.j)?;
    dict.set_item("image_j", candidate.image_j.to_vec())?;
    dict.set_item("distance", candidate.distance)?;
    Ok(dict.into_py(py))
}

#[pyfunction]
fn minimum_image_distance(
    py: Python<'_>,
    frac_i: [f64; 3],
    frac_j: [f64; 3],
    lattice: [[f64; 3]; 3],
) -> PyResult<PyObject> {
    image_distance_to_object(py, pbc::minimum_image_distance(frac_i, frac_j, lattice))
}

#[pyfunction]
fn find_short_contacts(
    py: Python<'_>,
    frac_coords: Vec<[f64; 3]>,
    atomic_numbers: Vec<u8>,
    lattice: [[f64; 3]; 3],
    cutoff_matrix: Vec<Vec<f64>>,
    scale: f64,
) -> PyResult<PyObject> {
    let contacts = contacts::find_short_contacts(
        &frac_coords,
        &atomic_numbers,
        lattice,
        &cutoff_matrix,
        scale,
    );
    let items = contacts
        .iter()
        .map(|contact| short_contact_to_object(py, contact))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyList::new_bound(py, items).into_py(py))
}

#[pyfunction]
fn find_neighbor_candidates(
    py: Python<'_>,
    frac_coords: Vec<[f64; 3]>,
    lattice: [[f64; 3]; 3],
    cutoff: f64,
) -> PyResult<PyObject> {
    let candidates = contacts::find_neighbor_candidates(&frac_coords, lattice, cutoff);
    let items = candidates
        .iter()
        .map(|candidate| neighbor_candidate_to_object(py, candidate))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyList::new_bound(py, items).into_py(py))
}

#[pyfunction]
fn connected_components(n_atoms: usize, edges: Vec<(usize, usize)>) -> PyResult<Vec<Vec<usize>>> {
    for &(a, b) in &edges {
        if a >= n_atoms || b >= n_atoms {
            return Err(PyValueError::new_err("edge endpoint out of bounds"));
        }
    }
    Ok(graph::connected_components(n_atoms, &edges))
}

#[pyfunction]
fn node_degrees(n_atoms: usize, edges: Vec<(usize, usize)>) -> PyResult<Vec<usize>> {
    for &(a, b) in &edges {
        if a >= n_atoms || b >= n_atoms {
            return Err(PyValueError::new_err("edge endpoint out of bounds"));
        }
    }
    Ok(graph::node_degrees(n_atoms, &edges))
}

#[pyfunction]
#[pyo3(signature = (n_atoms, edges, length_bound, max_cycles=0))]
fn bounded_simple_cycles_undirected(
    n_atoms: usize,
    edges: Vec<(usize, usize)>,
    length_bound: usize,
    max_cycles: usize,
) -> PyResult<Vec<Vec<usize>>> {
    for &(a, b) in &edges {
        if a >= n_atoms || b >= n_atoms {
            return Err(PyValueError::new_err("edge endpoint out of bounds"));
        }
    }
    graph::bounded_simple_cycles_undirected(n_atoms, &edges, length_bound, max_cycles).ok_or_else(
        || PyValueError::new_err(format!("cycle count exceeded max_cycles={max_cycles}")),
    )
}

#[pyfunction]
fn oms_is_open(cn: usize, center: [f64; 3], neighbor_coords: Vec<[f64; 3]>) -> PyResult<bool> {
    Ok(oms::oms_is_open(cn, center, &neighbor_coords))
}

#[pyfunction]
fn voronoi_center_neighbors(points: Vec<[f64; 3]>) -> PyResult<Vec<usize>> {
    voronoi::center_neighbors(&points).map_err(PyValueError::new_err)
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn eqeq_charges(
    frac_coords: Vec<[f64; 3]>,
    cell_lengths: [f64; 3],
    cell_angles: [f64; 3],
    electronegativity: Vec<f64>,
    hardness: Vec<f64>,
    total_charge: f64,
    lambda: f64,
    eta: f64,
    real_cells: i32,
    recip_cells: i32,
    charge_precision: i32,
) -> PyResult<Vec<f64>> {
    let input = eqeq::EqeqInput {
        frac_coords: &frac_coords,
        cell_lengths,
        cell_angles,
        electronegativity: &electronegativity,
        hardness: &hardness,
        total_charge,
        lambda,
        eta,
        real_cells,
        recip_cells,
        charge_precision,
    };
    eqeq::equilibrate(&input).map_err(PyValueError::new_err)
}

#[pymodule]
fn _rust(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(minimum_image_distance, m)?)?;
    m.add_function(wrap_pyfunction!(find_short_contacts, m)?)?;
    m.add_function(wrap_pyfunction!(find_neighbor_candidates, m)?)?;
    m.add_function(wrap_pyfunction!(connected_components, m)?)?;
    m.add_function(wrap_pyfunction!(node_degrees, m)?)?;
    m.add_function(wrap_pyfunction!(bounded_simple_cycles_undirected, m)?)?;
    m.add_function(wrap_pyfunction!(oms_is_open, m)?)?;
    m.add_function(wrap_pyfunction!(voronoi_center_neighbors, m)?)?;
    m.add_function(wrap_pyfunction!(eqeq_charges, m)?)?;
    Ok(())
}
