"""rustworkx-backed replacement for the `indices` output of
structuregraph_helpers.get_subgraphs_as_molecules(..., return_unique=False).

The reference builds a 3x3x3 supercell of the StructureGraph
(pymatgen StructureGraph.__mul__ -- networkx union/relabel + KD-tree, ~2.3 s/
struct) and keeps connected components that do NOT cross the *supercell*
boundary, then filters to molecules with an atom in the origin cell.

We reproduce that EXACTLY but with integer cell arithmetic + rustworkx
(no networkx union, no KD-tree):

- 27 cells c in {0,1,2}^3 (pymatgen's lattice_points_in_supercell order); the
  supercell node for (cell c, original site i) is `cell_order(c)*n + i`, whose
  pymatgen "idx" attribute is i.
- An original directed edge (u, v, J) replicated in cell c targets cell t = c+J.
  If t is inside the block it is an INTERNAL edge (supercell to_jimage 0); else
  it is a BOUNDARY-crossing edge (nonzero supercell to_jimage), which marks both
  its endpoints (node(c,u) and node(t mod 3, v)) as boundary-touching.
- A connected component (over internal edges) is a molecule iff it contains no
  boundary-touching node -- identical to pymatgen's "no edge with to_jimage !=
  (0,0,0)" test (a boundary edge excludes the components of both endpoints).
- Keep molecules with >=1 atom in the origin cell (original-lattice frac <= 1,
  i.e. cell c + frac_i <= 1 elementwise), deduped by sorted index set.

Validated to identical molecule-index-sets vs the reference; see
scripts/validate_subgraph_rx.py.
"""
from __future__ import annotations

import numpy as np
import rustworkx as rx
from pymatgen.util.coord import lattice_points_in_supercell


def _cells():
    """27 integer cell vectors in pymatgen's lattice_points_in_supercell order."""
    frac = lattice_points_in_supercell(3 * np.eye(3, dtype=int))  # supercell-frac
    cell_int = np.rint(frac * 3).astype(int)  # -> original-lattice integer units
    return cell_int


def finite_component_indices(structure_graph):
    """One sorted index list per molecule, matching
    get_subgraphs_as_molecules(return_unique=False)[2] as a set of index-sets."""
    n = len(structure_graph)
    frac = structure_graph.structure.frac_coords  # (n, 3)

    edges = []
    for u, v, d in structure_graph.graph.edges(data=True):
        j = d.get("to_jimage", (0, 0, 0))
        edges.append((int(u), int(v), np.array((int(j[0]), int(j[1]), int(j[2])))))

    cells = _cells()
    order = {tuple(c): ci for ci, c in enumerate(cells)}  # cell vec -> 0..26
    n_super = n * len(cells)

    def node(ci, i):
        return ci * n + i

    g = rx.PyGraph(multigraph=True)
    g.add_nodes_from(range(n_super))
    boundary = set()  # nodes incident to a boundary-crossing edge

    for ci, c in enumerate(cells):
        for u, v, J in edges:
            t = c + J
            tt = tuple(int(x) for x in t)
            if tt in order:  # target inside block -> internal edge
                g.add_edge(node(ci, u), node(order[tt], v), False)
            else:  # crosses supercell boundary: wrap target into block, keep as
                # connectivity (so multi-cell copies merge like pymatgen) but
                # flag both endpoints so the whole component is excluded.
                tw = tuple(int(x % 3) for x in t)
                a, b = node(ci, u), node(order[tw], v)
                g.add_edge(a, b, True)
                boundary.add(a)
                boundary.add(b)

    out = []
    seen = set()
    for comp in rx.connected_components(g):
        if comp & boundary:
            continue
        # in-origin-cell filter: any atom with (cell + frac_i) <= 1 elementwise
        in_cell = any(
            np.all(cells[nd // n] + frac[nd % n] <= 1.0) for nd in comp
        )
        if not in_cell:
            continue
        idx = sorted({nd % n for nd in comp})
        key = tuple(idx)
        if key not in seen:
            seen.add(key)
            out.append(idx)
    return out
