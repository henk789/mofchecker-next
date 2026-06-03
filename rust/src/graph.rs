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

pub fn node_degrees(n_atoms: usize, edges: &[(usize, usize)]) -> Vec<usize> {
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
    adjacency
        .into_iter()
        .map(|mut neighbors| {
            neighbors.sort_unstable();
            neighbors.dedup();
            neighbors.len()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::{connected_components, node_degrees};

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
}
