use std::f64::consts::PI;

const EPS: f64 = 1e-12;

#[derive(Clone, Copy)]
struct Param {
    ta: f64,
    igw_ta: f64,
    min_spp: f64,
    igw_spp: f64,
    igw_ep: f64,
    w_spp: f64,
    fac_aa: f64,
    exp_cos_aa: f64,
}

fn param(name: &str) -> Param {
    match name {
        "sq_plan" => Param {
            min_spp: 2.792526803190927,
            igw_spp: 15.0,
            igw_ep: 18.0,
            w_spp: 1.0,
            fac_aa: 1.0,
            exp_cos_aa: 2.0,
            ..Default::default()
        },
        "see_saw_rect" => Param {
            min_spp: 2.356194490192345,
            igw_spp: 11.5,
            igw_ep: 27.0,
            fac_aa: 2.0,
            exp_cos_aa: 2.0,
            ..Default::default()
        },
        "tet" => Param {
            ta: 0.6081734479693927,
            igw_ta: 15.0,
            fac_aa: 1.5,
            exp_cos_aa: 2.0,
            ..Default::default()
        },
        "tri_pyr" => Param {
            igw_ep: 15.5,
            fac_aa: 1.5,
            exp_cos_aa: 2.0,
            ..Default::default()
        },
        "pent_plan" => Param {
            ta: 0.6,
            igw_ta: 18.0,
            ..Default::default()
        },
        "sq_pyr" => Param {
            igw_ep: 14.9,
            fac_aa: 2.0,
            exp_cos_aa: 2.0,
            ..Default::default()
        },
        "tri_bipyr" => Param {
            min_spp: 2.356194490192345,
            igw_spp: 12.0,
            igw_ep: 16.6,
            fac_aa: 1.5,
            exp_cos_aa: 2.0,
            w_spp: 1.0,
            ..Default::default()
        },
        "pent_pyr" => Param {
            igw_ep: 13.8,
            fac_aa: 2.5,
            exp_cos_aa: 2.0,
            ..Default::default()
        },
        "oct" => Param {
            min_spp: 2.792526803190927,
            igw_spp: 15.0,
            igw_ep: 18.0,
            w_spp: 3.0,
            fac_aa: 2.0,
            exp_cos_aa: 2.0,
            ..Default::default()
        },
        "hex_pyr" => Param {
            igw_ep: 12.5,
            fac_aa: 3.0,
            exp_cos_aa: 2.0,
            ..Default::default()
        },
        "pent_bipyr" => Param {
            min_spp: 2.356194490192345,
            igw_spp: 12.5,
            igw_ep: 14.75,
            fac_aa: 2.5,
            exp_cos_aa: 2.0,
            w_spp: 1.0,
            ..Default::default()
        },
        "hex_bipyr" => Param {
            min_spp: 2.356194490192345,
            igw_spp: 14.1,
            igw_ep: 13.6,
            fac_aa: 3.0,
            exp_cos_aa: 2.0,
            w_spp: 1.0,
            ..Default::default()
        },
        _ => Default::default(),
    }
}

impl Default for Param {
    fn default() -> Self {
        Self {
            ta: 0.0,
            igw_ta: 0.0,
            min_spp: 0.0,
            igw_spp: 0.0,
            igw_ep: 0.0,
            w_spp: 1.0,
            fac_aa: 1.0,
            exp_cos_aa: 2.0,
        }
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
fn unit(a: [f64; 3]) -> Option<[f64; 3]> {
    let n = norm(a);
    (n > EPS).then(|| scale(a, 1.0 / n))
}
fn clamp(x: f64) -> f64 {
    x.max(-1.0).min(1.0)
}
fn gs(v: [f64; 3], z: [f64; 3]) -> [f64; 3] {
    sub(v, scale(z, dot(v, z)))
}

fn gaussian(x: f64) -> f64 {
    (-0.5 * x * x).exp()
}

fn geom_op(name: &str, rij_norm: &[[f64; 3]]) -> Option<f64> {
    let n = rij_norm.len();
    if n <= 1 {
        return None;
    }
    let p = param(name);
    let mut qsp = vec![vec![Vec::<f64>::new(); n]; n];
    let mut norms = vec![vec![Vec::<f64>::new(); n]; n];
    let ipi = 1.0 / PI;

    for j in 0..n {
        let zaxis = rij_norm[j];
        let mut kc = 0;
        for k in 0..n {
            if j == k {
                continue;
            }
            qsp[j][kc].push(0.0);
            norms[j][kc].push(0.0);
            let thetak = clamp(dot(zaxis, rij_norm[k])).acos();
            let xaxis0 = gs(rij_norm[k], zaxis);
            let xnorm = norm(xaxis0);
            let flag_xaxis = xnorm < EPS;
            let xaxis = if flag_xaxis {
                [0.0; 3]
            } else {
                scale(xaxis0, 1.0 / xnorm)
            };

            match name {
                "tet" => { /* needs m contribution too */ }
                "tri_pyr" | "sq_pyr" | "pent_pyr" | "hex_pyr" => {
                    qsp[j][kc][0] += gaussian(p.igw_ep * (thetak * ipi - 0.5));
                    norms[j][kc][0] += 1.0;
                }
                "sq_plan" | "oct" => {
                    if thetak >= p.min_spp {
                        qsp[j][kc][0] += p.w_spp * gaussian(p.igw_spp * (thetak * ipi - 1.0));
                        norms[j][kc][0] += p.w_spp;
                    }
                }
                "see_saw_rect" | "tri_bipyr" | "pent_bipyr" | "hex_bipyr" => {
                    if thetak < p.min_spp {
                        qsp[j][kc][0] += gaussian(p.igw_ep * (thetak * ipi - 0.5));
                        norms[j][kc][0] += 1.0;
                    }
                }
                "pent_plan" => { /* below */ }
                _ => {}
            }

            let mut gaussthetak = 0.0;
            if name == "tet" {
                gaussthetak = gaussian(p.igw_ta * (thetak * ipi - p.ta));
            } else if name == "pent_plan" {
                let tmp = if thetak <= p.ta * PI { 0.4 } else { 0.8 };
                gaussthetak = gaussian(p.igw_ta * (thetak * ipi - tmp));
            }

            for m in 0..n {
                if m == j || m == k || flag_xaxis {
                    continue;
                }
                let thetam = clamp(dot(zaxis, rij_norm[m])).acos();
                let xtwo0 = gs(rij_norm[m], zaxis);
                let xtwonorm = norm(xtwo0);
                if xtwonorm < EPS {
                    continue;
                }
                let xtwo = scale(xtwo0, 1.0 / xtwonorm);
                let phi = clamp(dot(xtwo, xaxis)).acos();
                let cos_term = (p.fac_aa * phi).cos().powf(p.exp_cos_aa);
                match name {
                    "tet" => {
                        qsp[j][kc][0] +=
                            gaussthetak * gaussian(p.igw_ta * (thetam * ipi - p.ta)) * cos_term;
                        norms[j][kc][0] += 1.0;
                    }
                    "pent_plan" => {
                        let tmp = if thetam <= p.ta * PI { 0.4 } else { 0.8 };
                        let c = phi.cos();
                        qsp[j][kc][0] +=
                            gaussthetak * gaussian(p.igw_ta * (thetam * ipi - tmp)) * c * c;
                        norms[j][kc][0] += 1.0;
                    }
                    "tri_pyr" | "sq_pyr" | "pent_pyr" | "hex_pyr" => {
                        qsp[j][kc][0] += cos_term * gaussian(p.igw_ep * (thetam * ipi - 0.5));
                        norms[j][kc][0] += 1.0;
                    }
                    "sq_plan" | "oct" => {
                        if thetak < p.min_spp && thetam < p.min_spp {
                            qsp[j][kc][0] += cos_term * gaussian(p.igw_ep * (thetam * ipi - 0.5));
                            norms[j][kc][0] += 1.0;
                        }
                    }
                    "tri_bipyr" | "pent_bipyr" | "hex_bipyr" => {
                        if thetam >= p.min_spp {
                            qsp[j][kc][0] += gaussian(p.igw_spp * (thetam * ipi - 1.0));
                            norms[j][kc][0] += 1.0;
                        }
                        if thetam < p.min_spp && thetak < p.min_spp {
                            qsp[j][kc][0] += cos_term * gaussian(p.igw_ep * (thetam * ipi - 0.5));
                            norms[j][kc][0] += 1.0;
                        }
                    }
                    "see_saw_rect" => {
                        if thetam < p.min_spp && thetak < p.min_spp && phi < 0.75 * PI {
                            qsp[j][kc][0] += cos_term * gaussian(p.igw_ep * (thetam * ipi - 0.5));
                            norms[j][kc][0] += 1.0;
                        }
                    }
                    _ => {}
                }
            }
            kc += 1;
        }
    }

    match name {
        "tet" | "sq_plan" | "oct" | "pent_plan" => {
            let mut sum = 0.0;
            let mut den = 0.0;
            for j in 0..n {
                for k in 0..qsp[j].len() {
                    sum += qsp[j][k].iter().sum::<f64>();
                    den += norms[j][k].iter().sum::<f64>();
                }
            }
            (den > EPS).then_some(sum / den)
        }
        _ => {
            let mut best: Option<f64> = None;
            for j in 0..n {
                for k in 0..qsp[j].len() {
                    let den: f64 = norms[j][k].iter().sum();
                    let val = if den > EPS {
                        qsp[j][k].iter().sum::<f64>() / den
                    } else {
                        0.0
                    };
                    best = Some(best.map_or(val, |b| b.max(val)));
                }
            }
            best
        }
    }
}

fn sq_op(center: [f64; 3], neigh: &[[f64; 3]], rij_norm: &[[f64; 3]]) -> Option<f64> {
    let n = neigh.len();
    if n < 3 {
        return None;
    }
    let mut aijs = Vec::new();
    for i in 0..rij_norm.len() {
        for j in (i + 1)..rij_norm.len() {
            aijs.push(clamp(dot(rij_norm[i], rij_norm[j])).acos());
        }
    }
    aijs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let centroid = neigh.iter().fold([0.0; 3], |acc, &v| {
        [acc[0] + v[0], acc[1] + v[1], acc[2] + v[2]]
    });
    let centroid = scale(centroid, 1.0 / n as f64);
    let h = norm(sub(centroid, center));
    let mut dists = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            dists.push(norm(sub(neigh[j], neigh[i])));
        }
    }
    if dists.is_empty() {
        return Some(0.0);
    }
    let b = dists.iter().copied().fold(f64::INFINITY, f64::min);
    let dhalf = dists.iter().copied().fold(0.0, f64::max) / 2.0;
    let a = 2.0 * (b / (2.0 * (h * h + dhalf * dhalf).sqrt())).asin();
    let mut op = 1.0;
    for angle in aijs.iter().take(n.min(4)) {
        op *= gaussian((angle - a) * 30.0);
    }
    Some(op)
}

pub fn oms_is_open(cn: usize, center: [f64; 3], neighbor_coords: &[[f64; 3]]) -> bool {
    if cn <= 3 {
        return true;
    }
    if cn > 8 {
        return false;
    }
    let rij_norm: Vec<[f64; 3]> = neighbor_coords
        .iter()
        .filter_map(|&v| unit(sub(v, center)))
        .collect();
    let ops: Vec<(Option<f64>, f64, bool)> = match cn {
        4 => vec![
            (geom_op("sq_plan", &rij_norm), 0.2, true),
            (sq_op(center, neighbor_coords, &rij_norm), 0.1, true),
            (geom_op("see_saw_rect", &rij_norm), 0.1, true),
            (geom_op("tet", &rij_norm), 0.5, false),
            (geom_op("tri_pyr", &rij_norm), 0.5, true),
        ],
        5 => vec![
            (geom_op("pent_plan", &rij_norm), 1.0, true),
            (geom_op("sq_pyr", &rij_norm), 0.5, true),
            (geom_op("tri_bipyr", &rij_norm), 0.5, false),
        ],
        6 => vec![
            (geom_op("pent_pyr", &rij_norm), 0.3, true),
            (geom_op("oct", &rij_norm), 0.7, false),
        ],
        7 => vec![
            (geom_op("hex_pyr", &rij_norm), 0.7, true),
            (geom_op("pent_bipyr", &rij_norm), 0.3, false),
        ],
        8 => vec![(geom_op("hex_bipyr", &rij_norm), 1.0, false)],
        _ => return false,
    };
    let mut open = 0.0;
    let mut total = 0.0;
    for (op, weight, is_open) in ops {
        if let Some(v) = op {
            let c = v * weight;
            total += c;
            if is_open {
                open += c;
            }
        }
    }
    total > EPS && open / total > 0.5
}

#[cfg(test)]
mod tests {
    use super::oms_is_open;

    #[test]
    fn low_cn_is_open_high_cn_is_closed() {
        assert!(oms_is_open(3, [0.0; 3], &[]));
        assert!(!oms_is_open(9, [0.0; 3], &[]));
    }
}
