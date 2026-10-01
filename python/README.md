# Python port

A Python (NumPy/SciPy) port of this repository's workflow:

| MATLAB | Python | Notes |
|---|---|---|
| `hindcast_run.m` | `hindcast_run.py` | same settings; command-line options for period/output |
| `initial_grid.m` | `initial_grid.py` | pandas + pyproj (EPSG:4167 → EPSG:2193) |
| `create_nc.m` | `create_nc.py` | contourpy + SciPy Delaunay interpolation; writes to `python/output/` |
| `project_output.m` | `project_output.py` | vectorised `get_polydistance`; matplotlib figure |
| `plot_waves.m` | `plot_waves.py` | reads the MATLAB v7.3 `wave_data.mat` or the Python one |
| `ShorelineS_functions/` | `shorelines/` | **the code path used by this project** (see below) |

```bash
pip install -r requirements.txt
git lfs pull                                   # data files (from the repository root)
python hindcast_run.py --end 2000-03-01        # short run; default is 2000-01-01 .. 2024-12-30
python create_nc.py                            # rebuild the wave input (to python/output/)
python -m pytest tests                         # tests against GNU Octave results
```

## Scope of the ShorelineS port

ShorelineS is ~18,000 lines of MATLAB with many optional processes. The port
covers the processes that `hindcast_run.m` switches on:

- open (non-cyclic) coastline sections separated by NaN, regridding
  (`griddingmethod` 2), coastline and foreshore orientation;
- wave time series at stations along a depth contour (`S.wvcfile`, NetCDF),
  `alongshoremapping` interpolation, time-series recycling;
- refraction/shoaling to `S.dnearshore`, iterative breaking height,
  critical wave angle (`get_Sphimax`), Kamphuis (`KAMP`) transport;
- shadowing by the coastline itself, smoothing at sharp angles, upwind
  correction for high-angle waves, `periodic`/`closed`/`neumann` boundaries;
- coastline change, overwash of narrow spits, removal of self-intersecting
  loops, merging of sections (when they do not cross), NaN cleanup;
- output storage with the same fields as `output.mat` (`O`).

Anything else (structures/groynes, revetments, dunes, mud, nourishments,
diffraction, tides, other transport formulas, wave climates, variable time
step, output grids) raises `NotImplementedError` rather than being ignored.

The port follows the MATLAB code statement by statement, including its
quirks, so that results can be compared exactly. Examples:
`get_Sphimax` fits its parabola with the transport of the *first* grid point
(`a = A\B` with a 3 × nq right-hand side, then `a(1)`, `a(2)`);
`get_foreshore_orientation` fixes the foreshore orientation at the first
time step; `spitdsf` uses the default `S.d`. MATLAB/Octave numerics
(`mod`, `sind`/`cosd`, `interp1`, colon ranges, `unique`, sequential `sum`)
are reproduced in `shorelines/mfuncs.py`.

## Verification against GNU Octave

The reference is the repository's MATLAB code run in GNU Octave 8.4.

- **Model, 2/6/30 days** (17/49/241 time steps): coastline coordinates
  (`O.x`, `O.y`, …) are **bit-identical** to Octave's. Transport and waves
  agree to floating-point round-off (relative differences ≤ 1e-14). See
  `tests/test_shorelines.py::test_model_2days_matches_octave`.
- **`get_intersections`**: bit-identical to Octave on 400 random cases
  (NaN separators, vertical segments, self-intersections).
- **`initial_grid`**: the initial coastline matches the one MATLAB stored in
  `30_sept_int_w/output.mat` to 4e-9 m (all 609 points).
- **`create_nc`**: station positions match MATLAB's `input_25_contour.nc` to
  6e-8 m; Hs and Tp to ≤ 1.4e-10; wave direction to 2.5e-7°.

NumPy's SIMD `exp`/`pow`/`atan2`/… can differ from glibc's libm (used by
Octave) in the last bit. With `SHORELINES_EXACT_LIBM=1` the port calls libm
instead (slower). This is only needed for bit-level comparisons.

RESULTS_PLACEHOLDER

## Things found while porting

- `initial_grid.m` reads `Range B2:KH517` with `readtable`, which takes row 2
  (the first CoastSat image, 1999-09-29) as the variable names. The removed
  transects are therefore held at the **second** image's position. The port
  reproduces this, because that is what the stored MATLAB run used.
- MATLAB's `contourc` places points only on grid-cell edges; contourpy's
  corner masking can add points on cell diagonals next to NaN cells. These
  are removed so that the 94 stations match MATLAB's.
