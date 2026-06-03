use crate::pbc::minimum_image_distance;

#[derive(Debug, Clone, PartialEq)]
pub struct ShortContact {
    pub i: usize,
    pub j: usize,
    pub image_j: [i32; 3],
    pub distance: f64,
    pub cutoff: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct NeighborCandidate {
    pub i: usize,
    pub j: usize,
    pub image_j: [i32; 3],
    pub distance: f64,
}

pub fn find_neighbor_candidates(
    frac_coords: &[[f64; 3]],
    lattice: [[f64; 3]; 3],
    cutoff: f64,
) -> Vec<NeighborCandidate> {
    let mut candidates = Vec::new();
    for i in 0..frac_coords.len().saturating_sub(1) {
        for j in (i + 1)..frac_coords.len() {
            let image_distance = minimum_image_distance(frac_coords[i], frac_coords[j], lattice);
            if image_distance.distance < cutoff {
                candidates.push(NeighborCandidate {
                    i,
                    j,
                    image_j: image_distance.image,
                    distance: image_distance.distance,
                });
            }
        }
    }

    candidates.sort_by(|a, b| {
        (a.i, a.j, a.image_j)
            .cmp(&(b.i, b.j, b.image_j))
            .then_with(|| a.distance.total_cmp(&b.distance))
    });
    candidates
}

pub fn find_short_contacts(
    frac_coords: &[[f64; 3]],
    atomic_numbers: &[u8],
    lattice: [[f64; 3]; 3],
    cutoff_matrix: &[Vec<f64>],
    scale: f64,
) -> Vec<ShortContact> {
    assert_eq!(frac_coords.len(), atomic_numbers.len());

    let mut contacts = Vec::new();
    for i in 0..frac_coords.len().saturating_sub(1) {
        let zi = atomic_numbers[i] as usize;
        for j in (i + 1)..frac_coords.len() {
            let zj = atomic_numbers[j] as usize;
            let cutoff = cutoff_matrix[zi][zj] * scale;
            let image_distance = minimum_image_distance(frac_coords[i], frac_coords[j], lattice);
            if image_distance.distance < cutoff {
                contacts.push(ShortContact {
                    i,
                    j,
                    image_j: image_distance.image,
                    distance: image_distance.distance,
                    cutoff,
                });
            }
        }
    }

    contacts.sort_by(|a, b| {
        (a.i, a.j, a.image_j)
            .cmp(&(b.i, b.j, b.image_j))
            .then_with(|| a.distance.total_cmp(&b.distance))
    });
    contacts
}

#[cfg(test)]
mod tests {
    use super::{find_neighbor_candidates, find_short_contacts};

    #[test]
    fn finds_single_contact() {
        let coords = [[0.0, 0.0, 0.0], [0.05, 0.0, 0.0], [0.5, 0.5, 0.5]];
        let atomic_numbers = [6, 6, 8];
        let lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
        let mut cutoffs = vec![vec![0.0; 9]; 9];
        cutoffs[6][6] = 0.76;
        cutoffs[6][8] = 0.66;
        cutoffs[8][6] = 0.66;
        let contacts = find_short_contacts(&coords, &atomic_numbers, lattice, &cutoffs, 1.0);
        assert_eq!(contacts.len(), 1);
        assert_eq!(contacts[0].i, 0);
        assert_eq!(contacts[0].j, 1);
        assert_eq!(contacts[0].image_j, [0, 0, 0]);
    }

    #[test]
    fn includes_periodic_image() {
        let coords = [[0.98, 0.0, 0.0], [1.02, 0.0, 0.0]];
        let atomic_numbers = [1, 1];
        let lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
        let mut cutoffs = vec![vec![0.0; 2]; 2];
        cutoffs[1][1] = 0.5;
        let contacts = find_short_contacts(&coords, &atomic_numbers, lattice, &cutoffs, 1.0);
        assert_eq!(contacts.len(), 1);
        assert_eq!(contacts[0].image_j, [0, 0, 0]);
        assert!((contacts[0].distance - 0.4).abs() < 1e-12);
    }

    #[test]
    fn neighbor_candidates_use_scalar_cutoff() {
        let coords = [[0.0, 0.0, 0.0], [0.04, 0.0, 0.0], [0.5, 0.5, 0.5]];
        let lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
        let candidates = find_neighbor_candidates(&coords, lattice, 0.5);
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].i, 0);
        assert_eq!(candidates[0].j, 1);
        assert_eq!(candidates[0].image_j, [0, 0, 0]);
        assert!((candidates[0].distance - 0.4).abs() < 1e-12);
    }
}
