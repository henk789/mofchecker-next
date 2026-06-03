# MOFChecker 2.0 Reference

The official MOFChecker 2.0 repository is used as a behavioral reference and parity oracle only.

Do not edit files under `external/mofchecker_2_ref/`.

## Source

- Remote URL: `https://github.com/Au-4/mofchecker_2.0`
- Local path: `external/mofchecker_2_ref/`
- Branch: `main`
- Commit: `a8076cd1dc854960eabe6b56fe08942c6f0d312c`
- Added as: plain clone, not git submodule, because this workspace was not a git repository when the reference was added.

Reproduce with:

```bash
mkdir -p external
git clone https://github.com/Au-4/mofchecker_2.0 external/mofchecker_2_ref
cd external/mofchecker_2_ref
git checkout a8076cd1dc854960eabe6b56fe08942c6f0d312c
```

## License

- License file: `external/mofchecker_2_ref/LICENSE`
- License name: Academic Non-Commercial Share-Alike License, ANCSA 1.0
- Copyright: 2025 Berend Smit
- Authors named in attribution clause: Xin Jin, Kevin Jablonka, Berend Smit

License implications for this project:

- Direct code copying is restricted by non-commercial and share-alike terms.
- Do not copy large code blocks into this project.
- Do not copy fixtures/examples unless license compatibility is explicitly accepted.
- Use the repository as a behavioral oracle, source of edge cases, and source-location guide.
- Prefer independently written tests and implementations based on observed behavior.

## Install Notes

Reference README recommends Python `<=3.9` and installation with:

```bash
pip install git+https://github.com/Au-4/mofchecker_2.0
```

Optional porosity support requires Zeo++:

```bash
conda install -c conda-forge zeopp-lsmo
```

`setup.py` declares Python `>=3.8` and dependencies:

- `pyeqeq`
- `click`
- `networkx>=2.5`
- `backports.cached-property`
- `ase`
- `pyyaml`
- `structuregraph_helpers`
- `element_coder`
- `typing_extensions`
- `libconeangle`
- `psutil`

The source imports additional packages including `numpy`, `scipy`, and `pymatgen`.

uv reference environment:

```bash
uv venv --python 3.9 .venv-ref
uv pip install --python .venv-ref/bin/python -e external/mofchecker_2_ref numpy scipy pymatgen pandas
```

Smoke tests:

```bash
.venv-ref/bin/python -c "from mofchecker import MOFChecker; import mofchecker; print(mofchecker.__version__)"
.venv-ref/bin/mofchecker --help
.venv-ref/bin/mofchecker --no-primitive -d name -d formula -d has_atomic_overlaps -d has_metal external/mofchecker_2_ref/test_cases/str_check/cifs/FONQIJ_clean.cif
```

Observed smoke-test output for `FONQIJ_clean.cif` descriptor subset:

```json
[
  {
    "name": "FONQIJ_clean",
    "formula": "Na4 Mn4 H64 C56 N8 O16",
    "has_atomic_overlaps": false,
    "has_metal": true
  }
]
```

Initial system environment check before creating `.venv-ref`:

- Python: `3.9.21`
- Present among checked dependencies: `click`, `yaml`
- Missing among checked dependencies: `pymatgen`, `ase`, `scipy`, `numpy`, `networkx`, `pyeqeq`, `structuregraph_helpers`, `element_coder`, `libconeangle`, `psutil`
- Therefore the reference checker is not currently runnable without installing dependencies.

## Relevant Source Files

- Main public API and descriptor orchestration: `src/mofchecker/__init__.py`
- CLI: `src/mofchecker/cli.py`
- Check base classes and flagged index conventions: `src/mofchecker/checks/check_base.py`
- Atomic overlap check: `src/mofchecker/checks/local_structure/overlapping_atoms.py`
- Local structure geometry helpers: `src/mofchecker/checks/local_structure/geometry.py`
- General geometry helpers: `src/mofchecker/checks/utils/geometry.py`
- Neighbor/index helpers: `src/mofchecker/checks/utils/get_indices.py`
- Radii/constants: `src/mofchecker/checks/data/definitions.py`
- Radii fallback behavior: `src/mofchecker/checks/data/__init__.py`
- Global graph/connectedness checks: `src/mofchecker/checks/global_structure/graphcheck.py`
- Floating solvent / lone molecule check: `src/mofchecker/checks/floating_solvent.py`
- Examples/test cases: `test_cases/`

## Initial Geometry/PBC Observations

- MOFChecker 2.0 uses `pymatgen.core.Structure` and `IStructure` internally.
- CIF loading is via `pymatgen.io.cif.CifParser` in `MOFChecker.from_cif`.
- The atomic overlap check uses `Structure.distance_matrix`, which is pymatgen's periodic distance matrix.
- Other local geometry checks use `Structure.get_distance(...)` and `structure.lattice.get_distance_and_image(...)` in specific routines.
- The overlap criterion is `dist < tolerance * min(covalent_radius_i, covalent_radius_j)` with default `tolerance = 1.0`.
- Overlap output is a list of flagged atom indices, not pair records.
- The overlap check does not expose periodic image information.
- Pymatgen neighbor finding is used directly in several checks through `Structure.get_neighbors(...)` and structure graphs via `structuregraph_helpers`.
- `has_3d_connected_graph` uses pymatgen's `get_dimensionality_larsen` on the Python `StructureGraph`.
- `has_lone_molecule` / floating solvent uses `structuregraph_helpers.subgraph.get_subgraphs_as_molecules` on the Python `StructureGraph` and returns flagged atom indices.
- The reference floating-solvent helper copies the structure graph, expands it to a 3x3x3 supercell, converts the graph to undirected connected components, and treats components with any nonzero `to_jimage` edge as boundary-crossing/non-molecular. It then returns index lists for non-periodic components.

Parity acceptance for the new atomic-overlap diagnostic slice:

- Match MOFChecker 2.0 on `has_atomic_overlaps`.
- Match MOFChecker 2.0 on the set of involved atom indices where `get_overlapping_indices()` is available.
- Treat contact pairs, distances, and `image_j` vectors as richer diagnostics added by `mofchecker-next`, not as fields that can be directly compared to MOFChecker 2.0.
- Preserve the `mofchecker-next` image convention: `r_ij` points from atom `i` to the selected periodic image of atom `j`, with `delta_frac = frac_j + image_j - frac_i`.
- Current tests compare `mofchecker-next` MIC distances against pymatgen `Lattice.get_distance_and_image`. Representative cases also match pymatgen's returned image vector, but our public convention remains defined by `r_ij`/`image_j` rather than by pymatgen internals.
- The Rust connected-components kernel accepts explicit Python-built edges only; it does not perform bond perception, dimensionality analysis, or CIF parsing.
- The current `mofchecker-next` non-periodic component helper accepts explicit edges plus explicit edge image vectors. It is diagnostic scaffolding and not yet a full parity port of `structuregraph_helpers.subgraph.get_subgraphs_as_molecules`.
- `structuregraph-helpers==0.0.9` package metadata reports MIT license. It is used as a production Python dependency for graph parity. Selected pieces may be evaluated for Rust ports later for speed, but only after parity tests lock behavior.

Overcoordination observations:

- `OverCoordinatedCarbonCheck` flags C atoms with graph coordination number `> 4`, unless any neighbor is a metal, and unless any neighbor is boron.
- `OverCoordinatedNitrogenCheck` flags N atoms with graph coordination number `> 4`, unless any neighbor is a metal.
- `OverCoordinatedHydrogenCheck` does not use graph degree; it calls `structure.get_neighbors(site, vdw_radius("H"))` and flags H atoms with more than one neighbor.
- The current `mofchecker-next` degree helper accepts explicit edges, target atomic number, max degree, and caller-supplied neighbor exclusions.
- The current H overcoordination helper mirrors the reference rule with pymatgen neighbor counting and a caller-supplied H VDW radius.

## Public API Notes

Python:

```python
from mofchecker import MOFChecker

mofchecker = MOFChecker.from_cif("structure.cif")
has_overlaps = mofchecker.has_atomic_overlaps
indices = mofchecker.get_overlapping_indices()
descriptors = mofchecker.get_mof_descriptors()
```

CLI:

```bash
mofchecker structure1.cif structure2.cif
mofchecker -d has_metal -d has_atomic_overlaps *.cif
```

Expected input:

- Standard CIF files.
- README recommends reading and rewriting CIFs with ASE before running.

Expected output:

- CLI prints a JSON list of descriptor dictionaries.
- Python `get_mof_descriptors()` returns an `OrderedDict` keyed by descriptor names.
- Individual checks expose booleans and, for index checks, flagged atom indices.

## Parity Plan Sketch

Start with a small wrapper script that runs one reference CIF and stores normalized JSON output. Do not implement broad parity yet.

Proposed script path:

- `scripts/run_reference_mofchecker.py`

Proposed behavior:

- Accept one CIF path and optional descriptor names.
- Import `MOFChecker` from the cloned reference path or installed package.
- For overlap parity, use `MOFChecker.from_cif(cif, primitive=False, symprec=None, angle_tolerance=None)` so the comparison targets the input structure and avoids reference failures from symmetry preprocessing on artificial fixtures.
- Emit normalized JSON with sorted keys and simple JSON types only.
- Include reference commit hash and descriptor list in the output metadata.
- With `--include-graph`, emit reference-built `StructureGraph` edges and edge image vectors for local graph-kernel parity scaffolding. These graph arrays are generated locally from the reference checkout and are not copied into committed fixtures.

Do not copy fixture CIFs from the reference repository into our tests unless license compatibility is approved.
