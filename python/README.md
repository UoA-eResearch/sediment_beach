# Python port

A Python (NumPy/SciPy) port of this repository's workflow:

| MATLAB | Python | Notes |
|---|---|---|
| `hindcast_run.m` | `hindcast_run.py` | command-line options for period/output/wave file, `--set KEY=VALUE` for model parameters; records the values used (`settings.json`, `S` in `output.mat`) |
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
python hindcast_run.py --set qscal=0.5 --set gamma=0.72                   # other model parameters
python hindcast_run.py --help                  # all options and parameters, with their defaults
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

- **Keys.** A parameter set can set any of the parameters below (MATLAB field
  names, case-insensitive), or any other ShorelineS setting the port supports.
  Parameters it does not set keep their defaults. `--runs` and `--set`
  combine: every row of the file is run with every `--set` combination.
  `python run_ensemble.py --help` (or `hindcast_run.py --help`) prints this
  table with the defaults computed from the code.
- **Checks first.** Every set is checked before anything starts: unknown keys,
  options the port does not support, missing wave files and start dates before
  1999-08-01 stop the command with an error.
- **Shared data.** `initial_grid` runs once. Each run takes its initial
  coastline from the row of its own start date.
- **Failures.** A run that fails does not stop the others. Its traceback goes
  to `error.txt` in its folder, and its row in `summary.csv` is marked
  `failed`.
- **Output** (in `--out-dir`, default `output/ensemble`). The values used are
  recorded in three places:
  - `run_<i>_<parameters>/output.mat` holds the model output `O`, and `S`
    with every ShorelineS setting used (defaults filled in), as in the
    MATLAB `output.mat`;
  - `run_<i>_<parameters>/params.json` holds the parameter set as given, the
    value of every parameter in the table, and every ShorelineS setting
    (plus `project_output.png` with `--plot`, and `error.txt` if the run
    failed);
  - `summary.csv` has one row per run: `run`, `status`, then the parameters
    that differ between runs, the statistics, and the value of every other
    parameter in the table. The statistics are `steps`, the run times, and
    the model-vs-CoastSat statistics from `project_output`: `rmse_m`,
    `bias_m` and `corr` over all compared points and times, and
    `model_mean_abs_change_end_m`.

  `ensemble.json` records the command's defaults and parameter sets.
  `hindcast_run.py` writes `output.mat` (with `S`) and `settings.json` to
  its output folder in the same way.

| Parameter | Default | Unit | Meaning |
|---|---|---|---|
| `start` | 1999-08-01 | date | start of the simulation; the initial coastline is the CoastSat/PCA coastline of this date |
| `end` | 2024-12-30 | date | end of the simulation |
| `wave_file` | input_25_contour_av.nc |  | wave input (NetCDF stations), e.g. ../sinthe_waves.nc for the synthetic waves |
| `interval` | 30 | days | storage interval of the output |
| `qscal` | 1 | - | calibration factor of the transport |
| `d50` | 0.0001 | m | median grain diameter |
| `gamma` | 0.72 | - | breaking coefficient Hs/h |
| `porosity` | 0.4 | - | porosity of the sediment |
| `tanbeta` | 0.03 | - | mean bed slope |
| `rhos` | 2650 | kg/m3 | density of sand |
| `rhow` | 1025 | kg/m3 | density of water |
| `g` | 9.81 | m/s2 | gravitational acceleration |
| `suppresshighangle` | 0 | 0/1 | limit the transport angle to the critical high angle |
| `twopoints` | 1 | 0/1 | high-angle correction spread over two downdrift cells |
| `maxangle` | 60 | deg | maximum coastline re-orientation between cells (smoothing of sharp angles) |
| `ds0` | 75 | m | alongshore grid size |
| `d` | 10 | m | active profile height |
| `boundaryconditionstart` | closed |  | start of open sections: closed, periodic, neumann, fixed (= neumann), angleconstant or gradient (= angleconstant); closed:Q fixes the transport at Q m3/yr, angleconstant:A the coastline angle at A degrees |
| `boundaryconditionend` | closed |  | end of open sections: as boundaryconditionstart |
| `ddeep` | 25 | m | water depth of the wave input |
| `dnearshore` | 25 | m | depth of the nearshore point (refraction/shoaling down to here) |
| `interpolationmethod` | alongshoremapping |  | wave stations to coast: alongshoremapping or weighted_distance |
| `smoothrefrac` | 0 | 0-1 | smoothing of the coastline angle used only for the refraction from the nearshore point to breaking (not for the transport) |
| `dt` | 1/(365·8) (3 h) | yr | fixed time step (1/(365*8) = 3 hours) |
| `spitwidth` | 50 | m | critical spit width below which overwash occurs |
| `owscale` | 0.1 | - | fraction of the overwash deficit moved per time step |
| `owtimescale` | 0 | yr | overwash time scale (used instead of owscale when > 0) |
| `bheight` | 2 | m | berm height used for overwash |
| `spitdsf` | 8 | m | shoreface depth used for the spit width |
| `spitdbb` | 4 | m | back-barrier depth used for the spit width |

Fixed in this configuration: `trform = KAMP`, `tc = 0`, `griddingmethod = 2`.

From Python: `run_ensemble.run_ensemble([{"qscal": 0.5}, {"qscal": 2}], "out", jobs=2)`
returns the summary rows.

**Measured** (4-CPU machine): 4 runs of 1999-08-01 to 2004-06-30 (14,361 steps
each; `qscal` 0.5 and 1 × synthetic and hindcast waves) with `-j 4` took
722 s in total. One of these runs alone takes 693 s, so the four cost about
the time of one: 3.3× faster than running them one after another. The
summary is in [`results/ensemble_1999-2004_summary.csv`](../results/ensemble_1999-2004_summary.csv).

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
  correction for high-angle waves;
- all ShorelineS boundary conditions of open sections: `closed` (or
  `{'Closed', Q}`, a fixed transport), `periodic`, `neumann`, `fixed` (the same
  as `neumann` in ShorelineS), `angleconstant`/`gradient` (or
  `{'Angleconstant', angle}`); on the command line `closed:Q`/`angleconstant:A`;
- `smoothrefrac`: smoothed coastline angle for the refraction from the
  nearshore point to breaking;
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
- **`smoothrefrac` and the boundary conditions, 2 days each**:
  - `smoothrefrac` 0.5 and 1;
  - `Angleconstant` at both ends;
  - `Fixed` at both ends;
  - `{'Closed', 5000}` with `{'Gradient', 200}` and `smoothrefrac` 0.3.

  In all five, the coastline is bit-identical to Octave's, and each option
  moves the coastline 0.7–4.9 m away from the plain closed run. See
  `tests/test_shorelines.py::test_smoothrefrac_and_boundary_conditions_match_octave`.
- **`get_intersections`**: bit-identical to Octave on 400 random cases
  (NaN separators, vertical segments, self-intersections).
- **Current `hindcast_run` settings, 2 days** (start 1999-08-01, closed
  boundaries, synthetic waves `sinthe_waves.nc`, `d50 = 1e-4`, `gamma = 0.5`):
  coastline coordinates are bit-identical to Octave's. Waves and transport
  agree to ≤ 1e-13, except at 3 points (3e-11 relative), where synthetic waves
  arrive from the land side (see above). See
  `tests/test_shorelines.py::test_model_synthetic_waves_closed_matches_octave`.
- **`initial_grid`** (with the earlier `A_site2.mat` PCA data): the initial coastline matches the one MATLAB stored in
  `30_sept_int_w/output.mat` to 4e-9 m (all 609 points). (That version kept
  the transects left out of the PCA; the current `initial_grid` drops them, as
  `initial_grid.m` does.)
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
