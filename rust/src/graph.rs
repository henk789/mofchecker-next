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

pub fn bounded_simple_cycles_undirected(
    n_atoms: usize,
    edges: &[(usize, usize)],
    length_bound: usize,
) -> Vec<Vec<usize>> {
    if length_bound < 3 {
        return Vec::new();
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
            dfs_cycles(start, &adjacency, length_bound, &mut seen, &mut path, &mut cycles);
            path.pop();
            seen[next] = false;
        }
        seen[start] = false;
    }
    cycles
}

fn dfs_cycles(
    start: usize,
    adjacency: &[Vec<usize>],
    length_bound: usize,
    seen: &mut [bool],
    path: &mut Vec<usize>,
    cycles: &mut Vec<Vec<usize>>,
) {
    let current = *path.last().unwrap();
    for &next in &adjacency[current] {
        if next == start {
            if path.len() >= 3 && path[1] < *path.last().unwrap() {
                cycles.push(path.clone());
            }
        } else if next > start && !seen[next] && path.len() < length_bound {
            seen[next] = true;
            path.push(next);
            dfs_cycles(start, adjacency, length_bound, seen, path, cycles);
            path.pop();
            seen[next] = false;
        }
    }
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
        let cycles = bounded_simple_cycles_undirected(4, &[(0, 1), (1, 2), (2, 0), (0, 3), (3, 2)], 4);
        assert_eq!(cycles, vec![vec![0, 1, 2], vec![0, 1, 2, 3], vec![0, 2, 3]]);
    }
}
