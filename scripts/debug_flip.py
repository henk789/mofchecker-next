"""For a flip case, determine GROUND TRUTH: is the component mine flags actually a
finite molecule (offset-consistent: no cycle with nonzero net lattice shift) or
periodic (my false positive)? Compares ref vs mine vs the rigorous test."""
import sys, warnings, glob
warnings.filterwarnings("ignore")
from pymatgen.core import Structure
from mofchecker_next.checks.graph import build_structure_graph
from mofchecker_next.checks._subgraph_rx import finite_component_indices
from structuregraph_helpers.subgraph import get_subgraphs_as_molecules
import rustworkx as rx


def offset_consistent(sg, comp_nodes):
    """True if the subgraph on comp_nodes admits a consistent integer image
    offset per node (=> genuinely FINITE). False => periodic."""
    comp = set(comp_nodes)
    adj = {i: [] for i in comp}
    for u, v, d in sg.graph.edges(data=True):
        if u in comp and v in comp:
            j = d.get("to_jimage", (0, 0, 0))
            adj[u].append((v, (j[0], j[1], j[2])))
            adj[v].append((u, (-j[0], -j[1], -j[2])))
    start = next(iter(comp))
    off = {start: (0, 0, 0)}
    stack = [start]
    while stack:
        x = stack.pop()
        ox = off[x]
        for nb, d in adj[x]:
            imp = (ox[0] + d[0], ox[1] + d[1], ox[2] + d[2])
            if nb in off:
                if off[nb] != imp:
                    return False
            else:
                off[nb] = imp
                stack.append(nb)
    return True


def main(path):
    s = Structure.from_file(path)
    sg = build_structure_graph(s, "vesta")
    ref = get_subgraphs_as_molecules(sg, return_unique=False)[2]
    mine = finite_component_indices(sg)
    print(f"\n{path}\n  n_sites={len(s)}  ref_molecules={len(ref)}  mine_molecules={len(mine)}")

    # undirected connected components of the ORIGINAL graph
    n = len(s)
    g = rx.PyGraph(multigraph=True); g.add_nodes_from(range(n))
    for u, v, d in sg.graph.edges(data=True):
        g.add_edge(int(u), int(v), None)
    comps = [sorted(c) for c in rx.connected_components(g)]
    print(f"  original-graph components: {len(comps)} (sizes {sorted(len(c) for c in comps)})")
    for mol in mine:
        ms = set(mol)
        # which original component is it
        parent = next((c for c in comps if ms <= set(c)), None)
        fin = offset_consistent(sg, parent) if parent else None
        print(f"  mine flags molecule idx={mol}")
        print(f"     -> its full original component size={len(parent) if parent else '?'}, "
              f"offset_consistent(FINITE)={fin}")
        in_ref = any(set(r) == ms for r in ref)
        print(f"     -> present in ref? {in_ref}")


if __name__ == "__main__":
    for p in sys.argv[1:]:
        for f in sorted(glob.glob(p)):
            main(f)
