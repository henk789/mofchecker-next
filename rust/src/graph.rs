pub fn connected_components(n_atoms: usize, edges: &[(usize, usize)]) -> Vec<Vec<usize>> {
    let mut adjacency = vec![Vec::new(); n_atoms];
    for &(a, b) in edges {
        assert!(a < n_atoms, "edge endpoint out of bounds");
        assert!(b < n_atoms, "edge endpoint out of bounds");
        if a == b {
            continue;
        }
        adjacency[a].push(b);
        adjacency[b].push(a);
    }
    for neighbors in &mut adjacency {
        neighbors.sort_unstable();
        neighbors.dedup();
    }

    let mut seen = vec![false; n_atoms];
    let mut components = Vec::new();
    for start in 0..n_atoms {
        if seen[start] {
            continue;
        }
        let mut stack = vec![start];
        seen[start] = true;
        let mut component = Vec::new();
        while let Some(node) = stack.pop() {
            component.push(node);
            for &neighbor in adjacency[node].iter().rev() {
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
        component.sort_unstable();
        components.push(component);
    }
    components.sort_by_key(|component| component[0]);
    components
}

fn adjacency(n_atoms: usize, edges: &[(usize, usize)]) -> Vec<Vec<usize>> {
    let mut adjacency = vec![Vec::new(); n_atoms];
    for &(a, b) in edges {
        assert!(a < n_atoms, "edge endpoint out of bounds");
        assert!(b < n_atoms, "edge endpoint out of bounds");
        if a == b {
            continue;
        }
        adjacency[a].push(b);
        adjacency[b].push(a);
    }
    for neighbors in &mut adjacency {
        neighbors.sort_unstable();
        neighbors.dedup();
    }
    adjacency
}

pub fn node_degrees(n_atoms: usize, edges: &[(usize, usize)]) -> Vec<usize> {
    adjacency(n_atoms, edges)
        .into_iter()
        .map(|neighbors| neighbors.len())
        .collect()
}

/// Enumerate all simple cycles up to `length_bound`. Returns `None` if the
/// number of cycles would exceed `max_cycles` (0 = unlimited). The count cap
/// bounds memory analytically: dense/degenerate graphs (e.g. structures with
/// overlapping atoms) have exponentially many short cycles, which otherwise
/// grows this Vec without bound and OOM-kills the process. A real MOF linker
/// set has at most a few thousand nonmetal cycles, so exceeding the cap means
/// the structure is pathological and ring/charge checks are not meaningful.
pub fn bounded_simple_cycles_undirected(
    n_atoms: usize,
    edges: &[(usize, usize)],
    length_bound: usize,
    max_cycles: usize,
) -> Option<Vec<Vec<usize>>> {
    if length_bound < 3 {
        return Some(Vec::new());
    }
    let adjacency = adjacency(n_atoms, edges);
    let mut cycles = Vec::new();
    let mut path = Vec::with_capacity(length_bound);
    let mut seen = vec![false; n_atoms];

    for start in 0..n_atoms {
        if adjacency[start].len() < 2 {
            continue;
        }
        path.clear();
        path.push(start);
        seen[start] = true;
        for &next in &adjacency[start] {
            if next <= start {
                continue;
            }
            seen[next] = true;
            path.push(next);
            let aborted =
                dfs_cycles(start, &adjacency, length_bound, max_cycles, &mut seen, &mut path, &mut cycles);
            path.pop();
            seen[next] = false;
            if aborted {
                return None;
            }
        }
        seen[start] = false;
    }
    Some(cycles)
}

/// Returns true if the cycle cap was hit and enumeration should abort.
fn dfs_cycles(
    start: usize,
    adjacency: &[Vec<usize>],
    length_bound: usize,
    max_cycles: usize,
    seen: &mut [bool],
    path: &mut Vec<usize>,
    cycles: &mut Vec<Vec<usize>>,
) -> bool {
    let current = *path.last().unwrap();
    for &next in &adjacency[current] {
        if next == start {
            if path.len() >= 3 && path[1] < *path.last().unwrap() {
                cycles.push(path.clone());
                if max_cycles != 0 && cycles.len() > max_cycles {
                    return true;
                }
            }
        } else if next > start && !seen[next] && path.len() < length_bound {
            seen[next] = true;
            path.push(next);
            let aborted = dfs_cycles(start, adjacency, length_bound, max_cycles, seen, path, cycles);
            path.pop();
            seen[next] = false;
            if aborted {
                return true;
            }
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::{bounded_simple_cycles_undirected, connected_components, node_degrees};

    #[test]
    fn handles_empty_graph() {
        let components = connected_components(0, &[]);
        assert_eq!(components, Vec::<Vec<usize>>::new());
    }

    #[test]
    fn returns_isolated_atoms() {
        let components = connected_components(3, &[]);
        assert_eq!(components, vec![vec![0], vec![1], vec![2]]);
    }

    #[test]
    fn finds_multiple_components_stably() {
        let components = connected_components(6, &[(4, 5), (1, 2), (2, 3)]);
        assert_eq!(components, vec![vec![0], vec![1, 2, 3], vec![4, 5]]);
    }

    #[test]
    fn ignores_self_edges_and_duplicate_edges() {
        let components = connected_components(3, &[(0, 0), (0, 1), (1, 0), (1, 2)]);
        assert_eq!(components, vec![vec![0, 1, 2]]);
    }

    #[test]
    fn computes_unique_node_degrees() {
        let degrees = node_degrees(4, &[(0, 1), (1, 0), (1, 2), (2, 2)]);
        assert_eq!(degrees, vec![1, 2, 1, 0]);
    }

    #[test]
    fn bounded_cycles_are_unique() {
        let cycles =
            bounded_simple_cycles_undirected(4, &[(0, 1), (1, 2), (2, 0), (0, 3), (3, 2)], 4, 0).unwrap();
        assert_eq!(cycles, vec![vec![0, 1, 2], vec![0, 1, 2, 3], vec![0, 2, 3]]);
    }

    #[test]
    fn cycle_cap_truncates() {
        // K4 has 7 simple cycles; cap of 3 must abort.
        let k4 = &[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
        assert!(bounded_simple_cycles_undirected(4, k4, 4, 3).is_none());
        assert!(bounded_simple_cycles_undirected(4, k4, 4, 0).is_some());
    }
}
