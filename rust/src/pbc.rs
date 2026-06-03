#[derive(Debug, Clone, PartialEq)]
pub struct ImageDistance {
    pub delta_frac: [f64; 3],
    pub image: [i32; 3],
    pub distance: f64,
}

fn round_ties_even(value: f64) -> f64 {
    let rounded = value.round();
    let diff = (value - rounded).abs();
    if diff == 0.5 {
        let lower = value.floor();
        let upper = value.ceil();
        if (lower as i64).rem_euclid(2) == 0 {
            lower
        } else {
            upper
        }
    } else {
        rounded
    }
}

pub fn minimum_image_distance(
    frac_i: [f64; 3],
    frac_j: [f64; 3],
    lattice: [[f64; 3]; 3],
) -> ImageDistance {
    let mut base_image = [0; 3];
    let mut d = [0.0; 3];
    for k in 0..3 {
        d[k] = frac_j[k] - frac_i[k];
        base_image[k] = -round_ties_even(d[k]) as i32;
    }

    let mut best = image_distance_for_image(d, base_image, lattice);
    for x in -1..=1 {
        for y in -1..=1 {
            for z in -1..=1 {
                if x == 0 && y == 0 && z == 0 {
                    continue;
                }
                let image = [base_image[0] + x, base_image[1] + y, base_image[2] + z];
                let candidate = image_distance_for_image(d, image, lattice);
                if candidate.distance < best.distance {
                    best = candidate;
                }
            }
        }
    }

    best
}

fn image_distance_for_image(
    d: [f64; 3],
    image: [i32; 3],
    lattice: [[f64; 3]; 3],
) -> ImageDistance {
    let delta_frac = [
        d[0] + image[0] as f64,
        d[1] + image[1] as f64,
        d[2] + image[2] as f64,
    ];
    let mut delta_cart = [0.0; 3];
    for col in 0..3 {
        delta_cart[col] = delta_frac[0] * lattice[0][col]
            + delta_frac[1] * lattice[1][col]
            + delta_frac[2] * lattice[2][col];
    }
    let distance = (delta_cart[0] * delta_cart[0]
        + delta_cart[1] * delta_cart[1]
        + delta_cart[2] * delta_cart[2])
        .sqrt();

    ImageDistance {
        delta_frac,
        image,
        distance,
    }
}

#[cfg(test)]
mod tests {
    use super::minimum_image_distance;

    fn assert_close(left: f64, right: f64, tolerance: f64) {
        assert!((left - right).abs() <= tolerance, "{left} != {right}");
    }

    #[test]
    fn identical_positions() {
        let lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
        let result = minimum_image_distance([0.1, 0.2, 0.3], [0.1, 0.2, 0.3], lattice);
        assert_eq!(result.delta_frac, [0.0, 0.0, 0.0]);
        assert_eq!(result.image, [0, 0, 0]);
        assert_close(result.distance, 0.0, 1e-12);
    }

    #[test]
    fn cubic_no_wrapping() {
        let lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
        let result = minimum_image_distance([0.1, 0.1, 0.1], [0.2, 0.3, 0.4], lattice);
        assert_eq!(result.image, [0, 0, 0]);
        assert_close(result.distance, (1.0_f64 + 4.0 + 9.0).sqrt(), 1e-12);
    }

    #[test]
    fn cubic_x_boundary_wrapping() {
        let lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
        let result = minimum_image_distance([0.95, 0.0, 0.0], [0.05, 0.0, 0.0], lattice);
        assert_eq!(result.image, [1, 0, 0]);
        assert_close(result.delta_frac[0], 0.1, 1e-12);
        assert_close(result.distance, 1.0, 1e-12);
    }

    #[test]
    fn y_z_wrapping() {
        let lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
        let result = minimum_image_distance([0.0, 0.9, 0.95], [0.0, 0.1, 0.05], lattice);
        assert_eq!(result.image, [0, 1, 1]);
        assert_close(result.distance, (2.0_f64 * 2.0 + 1.0 * 1.0).sqrt(), 1e-12);
    }

    #[test]
    fn non_orthogonal_lattice() {
        let lattice = [[4.0, 0.0, 0.0], [1.0, 3.0, 0.0], [0.5, 0.25, 5.0]];
        let result = minimum_image_distance([0.0, 0.0, 0.0], [0.25, 0.25, 0.0], lattice);
        assert_eq!(result.image, [0, 0, 0]);
        assert_close(result.distance, (1.25_f64 * 1.25 + 0.75 * 0.75).sqrt(), 1e-12);
    }

    #[test]
    fn skew_lattice_avoids_component_wrapping_pitfall() {
        let lattice = [[1.0, 0.0, 0.0], [0.99, 0.1, 0.0], [0.0, 0.0, 10.0]];
        let result = minimum_image_distance([0.0, 0.0, 0.0], [0.49, 0.49, 0.0], lattice);
        assert_eq!(result.image, [0, -1, 0]);
        assert_close(result.delta_frac[0], 0.49, 1e-12);
        assert_close(result.delta_frac[1], -0.51, 1e-12);
        assert!(result.distance < 0.06);
    }

    #[test]
    fn negative_fractional_coordinates() {
        let lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
        let result = minimum_image_distance([-0.1, 0.0, 0.0], [0.1, 0.0, 0.0], lattice);
        assert_eq!(result.image, [0, 0, 0]);
        assert_close(result.distance, 2.0, 1e-12);
    }

    #[test]
    fn coordinates_greater_than_one() {
        let lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
        let result = minimum_image_distance([1.1, 0.0, 0.0], [1.9, 0.0, 0.0], lattice);
        assert_eq!(result.image, [-1, 0, 0]);
        assert_close(result.delta_frac[0], -0.2, 1e-12);
        assert_close(result.distance, 2.0, 1e-12);
    }
}
