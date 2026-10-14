# Python port

A Python (NumPy/SciPy) port of this repository's workflow:

| MATLAB | Python | Notes |
|---|---|---|
| `hindcast_run.m` | `hindcast_run.py` | command-line options for period/output/wave file; `hindcast_settings()` returns the settings |
| `run_ensemble.m` | `run_ensemble.py` | several parameter sets in parallel, one process per run |
| `initial_grid.m` | `initial_grid.py` | pandas + pyproj (EPSG:4167 → EPSG:2193) |
| `create_nc.m` | `create_nc.py` | contourpy + SciPy Delaunay interpolation; writes to `python/output/` |
| `project_output.m` | `project_output.py` | vectorised `get_polydistance`; matplotlib figure |
| `plot_waves.m` | `plot_waves.py` | reads the MATLAB v7.3 `wave_data.mat` or the Python one |
| – | `make_figures.py` | the README figures (`docs/figures/`) from `results/` |
| – | `compare_model_coastsat.py` | the statistics behind the README section on why the model and CoastSat change differently |
| `ShorelineS_functions/` | `shorelines/` | **the code path used by this project** (see below) |

```bash
pip install -r requirements.txt
git lfs pull                                   # data files (from the repository root)
python hindcast_run.py --end 2000-03-01        # short run; default is 1999-08-01 .. 2024-12-30
python hindcast_run.py --wave-file ../sinthe_waves.nc --end 2004-06-30   # synthetic waves
python run_ensemble.py --set qscal=0.5,1,2 --end 2004-06-30 -j 3          # 3 runs in parallel
python create_nc.py                            # rebuild the wave input (to python/output/)
python -m pytest tests                         # tests against GNU Octave results
```

The reference data in `tests/data/` (`.mat`) are stored with Git LFS, like the
other data files.

## Running several parameter sets in parallel

One run cannot use more than one CPU (the time loop is sequential), so
`run_ensemble.py` runs parameter sets side by side, one process each
(`concurrent.futures.ProcessPoolExecutor`, one BLAS thread per process).

```bash
# grid: every combination of the --set values (3 x 2 = 6 runs), 4 processes at a time
python run_ensemble.py --set qscal=0.5,1,2 --set gamma=0.5,0.72 --end 2004-06-30 -j 4

# a list of parameter sets: JSON list of objects, or CSV with one set per row
python run_ensemble.py --runs runs.csv --out-dir output/my_ensemble --plot
```

```csv
qscal,d,boundaryconditionstart,boundaryconditionend,wave_file
0.5,8,closed,closed,../sinthe_waves.nc
1,10,periodic,periodic,../input_25_contour.nc
```

- **Keys.** A parameter set can override any ShorelineS setting the port
  supports (MATLAB field names, case-insensitive: `qscal`, `d`, `d50`,
  `gamma`, `ds0`, `ddeep`, `dnearshore`, `boundaryconditionstart`, …), plus
  `start`, `end`, `wave_file` and `interval` (storage interval in days).
  `--runs` and `--set` combine: every row of the file is run with every `--set`
  combination.
- **Checks first.** Every set is checked before anything starts: unknown keys,
  options the port does not support, missing wave files and start dates before
  1999-08-01 stop the command with an error.
- **Shared data.** `initial_grid` runs once. Each run takes its initial
  coastline from the row of its own start date.
- **Failures.** A run that fails does not stop the others. Its traceback goes
  to `error.txt` in its folder, and its row in `summary.csv` is marked
  `failed`.
- **Output** (in `--out-dir`, default `output/ensemble`):
  - `run_<i>_<parameters>/` holds `output.mat` and `params.json`, plus
    `project_output.png` with `--plot`;
  - `summary.csv` has one row per run: the parameters, `status`, `steps`, the
    run time, and model-vs-CoastSat statistics from `project_output`
    (`rmse_m`, `bias_m`, `corr` over all compared points and times, and
    `model_mean_abs_change_end_m`);
  - `ensemble.json` records the parameter sets.

From Python: `run_ensemble.run_ensemble([{"qscal": 0.5}, {"qscal": 2}], "out", jobs=2)`
returns the summary rows.

In MATLAB, `run_ensemble.m` does the same with `parfor`, which needs the
Parallel Computing Toolbox (start a pool with `parpool(n)`); without it the
loop runs serially:

```matlab
S0 = ...;   % the settings of hindcast_run.m, with S0.xmc / S0.ymc from initial_grid
params = struct('qscal', {0.5, 1, 2});
run_ensemble(params, S0, 'ensemble');
```

## Scope of the ShorelineS port

ShorelineS is ~18,000 lines of MATLAB with many optional processes. The port
covers the processes that `hindcast_run.m` switches on, plus closed
(cyclic) coastline sections:

- coastline sections separated by NaN, open or closed (see below),
  regridding (`griddingmethod` 2), coastline and foreshore orientation;
- wave time series at stations along a depth contour (`S.wvcfile`, NetCDF),
  `alongshoremapping` interpolation, time-series recycling;
- refraction/shoaling to `S.dnearshore`, iterative breaking height,
  critical wave angle (`get_Sphimax`), Kamphuis (`KAMP`) transport;
- shadowing by the coastline itself, smoothing at sharp angles, upwind
  correction for high-angle waves, `periodic`/`closed`/`neumann` boundaries;
- coastline change, overwash of narrow spits, removal of self-intersecting
  loops, merging of sections (also sections that cross each other, e.g. two
  beaches growing into each other: `shorelines/merge.py`), NaN cleanup;
- output storage with the same fields as `output.mat` (`O`).

### Closed (cyclic) sections: islands and lakes

A section is closed when its last point is within 0.2 · `ds0` of its first
point (ShorelineS' `get_cyclic`). Give it as a polygon that repeats its first
point at the end:

```python
x_mc = np.r_[x_beach, np.nan, x_island]   # x_island[-1] == x_island[0]
y_mc = np.r_[y_beach, np.nan, y_island]
```

- **Orientation matters.** The land must be on the right-hand side, as for
  open coasts (sea on the left when walking along the coastline). An island
  is therefore **clockwise**, and a lake or lagoon (water inside) is
  counter-clockwise. As in ShorelineS, the smoothing at sharp angles and the
  high-angle upwind correction are applied only to clockwise sections.
- **Differences from open sections:**
  - grid spacing, transport points and orientation wrap around the closure;
  - the foreshore orientation is smoothed around the loop and mapped with
    inverse-distance weights;
  - no boundary condition is applied;
  - the first and last point are moved together, in `coastline_change`,
    `merge_coastlines` and overwash.
- **Ordering.** When a closed section comes before an open one,
  `merge_coastlines_mc` swaps them after the first time step. The MATLAB code
  does the same.
- **Crossings.** Sections that cross each other (an island growing into the
  coast, two beaches growing into each other) are merged as in ShorelineS'
  `merge_coastlines_mc`.

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
- **Closed sections, 6 days** (49 steps, daily output), in four layouts:
  - the two beaches plus an island;
  - the island listed first;
  - the island alone;
  - the beaches plus a thin 35–70 m wide barrier island, so that overwash is
    active.

  Coastline coordinates are bit-identical to Octave's at every stored time.
  Waves and transport agree to ≤ 1e-14 relative, except at 13 lee-side
  points of the thin island, where they agree to 3e-11. There, waves
  arrive from the land side, the angle is clipped to 89.999°, and dividing
  by cos(89.999°) ≈ 2e-5 magnifies rounding differences. A separate test
  exercises the first/last-point merge of closed sections during overwash,
  and `get_clockpoly` is checked on 300 random polygons.
- **Merging of crossing sections** (`merge_coastlines_mc`): bit-identical to
  Octave on 400 random pairs of crossing sections:
  - two open beaches;
  - a beach and an island, in either order;
  - two islands;
  - wavy lines crossing many times.

  305 of these cases merge. See
  `tests/test_shorelines.py::test_merge_crossing_sections_matches_octave`.
- **`get_intersections`**: bit-identical to Octave on 400 random cases
  (NaN separators, vertical segments, self-intersections).
- **Current `hindcast_run` settings, 2 days** (start 1999-08-01, closed
  boundaries, synthetic waves `sinthe_waves.nc`, `d50 = 1e-4`, `gamma = 0.5`):
  coastline coordinates are bit-identical to Octave's. Waves and transport
  agree to ≤ 1e-13, except at 3 points (3e-11 relative), where synthetic waves
  arrive from the land side (see above). See
  `tests/test_shorelines.py::test_model_synthetic_waves_closed_matches_octave`.
- **`initial_grid`** (with the earlier `A_site2.mat` PCA data): the initial coastline matches the one MATLAB stored in
  `30_sept_int_w/output.mat` to 4e-9 m (all 609 points).
- **`create_nc`**: station positions match MATLAB's `input_25_contour.nc` to
  6e-8 m; Hs and Tp to ≤ 1.4e-10; wave direction to 2.5e-7°.

NumPy's SIMD `exp`/`pow`/`atan2`/… can differ from glibc's libm (used by
Octave) in the last bit. With `SHORELINES_EXACT_LIBM=1` the port calls libm
instead (slower). This is only needed for bit-level comparisons.

**Against the real MATLAB run.** `30_sept_int_w/output.mat` is a run of the
same configuration made in MATLAB (output up to 2004-03-10). Over those
4 years and 2 months the Python coastline uses the same grid at every stored
time, and positions differ by at most 1e-7 m (median 9e-10 m in cross-shore
distance). There is no drift.

## Runtime

Same 4-core machine, each run alone, GNU Octave 8.4 vs. Python 3.11 /
NumPy 2.4:

| | Octave, original MATLAB code | Octave, optimised MATLAB code | Python |
|---|---:|---:|---:|
| Model, 2 days (17 steps) | 171.0 s | 52.4 s | 0.9 s |
| Model, 6 days (49 steps) | 455.3 s | 104.3 s | 2.1 s |
| Model, 30 days (241 steps) | – | 369.8 s | 9.0 s |
| Model, per time step | 8.9 s | 1.4–1.6 s | 0.035 s |
| **Model, full 2000–2024 (73,041 steps)** | ≈ 180 h (est.) | ≈ 30 h (est.) | **43.9 min (measured)** |
| `project_output`, 9,131 daily coastlines | ≈ 2.1 h (est.) | 26.5 s | 25.8 s |
| `initial_grid` | – | – | 4.8 s |
| `create_nc` | – | – | 6 s |

Runtime does not depend on the number of CPUs: the time loop is sequential,
and the per-step arrays are too small for multithreaded BLAS (measured: same
time on 1 and 3 cores).

**Why Python is faster:** most of the difference comes from interpreter
overhead. The MATLAB model calls small functions hundreds of thousands of
times per simulated day (`sind`, `interp1`, `unique`, and `get_intersections`
for every coastline point in `find_shadows_mc` / `find_overwash_mc`), and
Octave, which has no JIT, pays 10–100 µs per call. The port batches these
(e.g. all shadow rays at once: 276 ms per call in Octave vs. 3 ms in Python).
Where both versions are vectorised (`project_output`), they run at the same
speed. MATLAB's JIT should make the MATLAB code considerably faster than in
Octave; MATLAB was not available to measure it.

## Things found while porting

- `initial_grid.m` reads `Range B2:KH517` with `readtable`, which takes row 2
  (the first CoastSat image, 1999-09-29) as the variable names. The removed
  transects are therefore held at the **second** image's position. The port
  reproduces this, because that is what the stored MATLAB run used.
- MATLAB's `contourc` places points only on grid-cell edges; contourpy's
  corner masking can add points on cell diagonals next to NaN cells. These
  are removed so that the 94 stations match MATLAB's.
