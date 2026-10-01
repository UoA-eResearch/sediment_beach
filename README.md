# sediment_beach

MATLAB workflow for hindcasting shoreline change on a New Zealand embayment
("site 2", Bay of Plenty) with the [ShorelineS](https://github.com/Deltares/ShorelineS)
one-line shoreline model, and comparing the model against CoastSat satellite
shorelines.

The workflow:

1. interpolates an Oceanum wave hindcast onto the 25 m depth contour
   (`create_nc.m`);
2. rebuilds the observed CoastSat shoreline of the selected beaches from a PCA
   reconstruction and uses it as the model's initial coastline
   (`initial_grid.m`);
3. runs ShorelineS for 2000-01-01 to 2024-12-30 at a 3-hour time step
   (`hindcast_run.m`);
4. projects modelled and observed shorelines onto the initial coastline and
   plots both as cross-shore change (`project_output.m`).

---

## Requirements

| Requirement | Used for |
|---|---|
| MATLAB R2016b or newer (implicit expansion) | everything |
| Mapping Toolbox (`projcrs`, `projfwd`) | `initial_grid.m`, `create_nc.m` |
| `readtable` (base MATLAB) | `initial_grid.m` |
| `VideoWriter` (base MATLAB) | animation when `S.video=1` |
| [Git LFS](https://git-lfs.com) | the `.mat` and `.nc` data files |

The `.mat` and `.nc` files are stored with Git LFS. After cloning, run
`git lfs pull`. Otherwise these files are small text pointer files and
loading them fails.

The ShorelineS model core (`ShorelineS_functions/ShorelineS.m` and the time
loop) also runs in GNU Octave 8 with the `netcdf` package. The benchmarks
below were made that way. The project scripts need MATLAB, because Octave has
no `readtable`, `projcrs`/`projfwd` or `scatteredInterpolant`.

---

## Quick start

(For the Python version, see [python/README.md](python/README.md).)

```matlab
% from the repository root, in MATLAB
git lfs pull                  % (in a shell) fetch the data files once
create_nc                     % optional: rebuilds input_25_contour.nc + wave_data.mat
hindcast_run                  % runs the model and calls project_output at the end
plot_waves                    % optional: plots the wave time series along the contour
```

`hindcast_run.m` adds `ShorelineS_functions/` to the path relative to its own
location, so it can be run from any working directory.

---

## Repository contents

### Scripts and functions (project code)

| File | Type | What it does |
|---|---|---|
| `hindcast_run.m` | script | **Main entry point.** Sets all ShorelineS parameters in a struct `S`, builds the initial coastline with `initial_grid`, runs `ShorelineS(S)` and then `project_output`. |
| `initial_grid.m` | function | `[x, y, x_initial, y_initial, g] = initial_grid(beach, dy, n_d)`. For each selected beach it reads the transects from the CoastSat workbook (`nzdXXXX.xlsx`, sheet *Transects*), converts transect origins to NZTM2000 (EPSG:2193), and rebuilds shoreline positions as `origin + (A + mint) * [sin(phi), cos(phi)]` from the PCA reconstruction (`A_site2.mat`, `mint_site2.mat`). Transects left out of the PCA (`remove{b}`) are held fixed at their first valid CoastSat observation (sheet *Intersect points*). Beaches are joined with one NaN column between them. Returns the full daily time series (`x`, `y`: time x transects) and the first row as the initial coastline. `dy` and `n_d` are not used here; trimming happens in `project_output`. |
| `project_output.m` | function | `project_output(x, y, file_out, dy, interv, x_initial, y_initial, n_d)`. Loads `<file_out>/output.mat`, selects observations `dy : dy+n_d`, and uses `get_polydistance` (normals of ±500 m around the initial coastline) to compute cross-shore change for the model, for the observations every `interv` days, and for the daily observations. Plots model, CoastSat, and CoastSat with a 30-day moving mean. |
| `create_nc.m` | script | Preprocesses the Oceanum hindcast `oceanum_paper2.nc` (48 nodes, hourly from 1979). Selects 1999-08-01 to 2024-12-30, finds the 25 m depth contour, interpolates Hs, Tp and mean direction (direction through its sin/cos components) linearly onto the contour points, with nearest-neighbour extrapolation, and writes `input_25_contour.nc` in the ShorelineS wave-station format (`time` in seconds since 1999-08-01, `station_x/y` in NZTM, `point_hm0`, `point_tp`, `point_wavdir`). Also saves `wave_data.mat`. Note that `nccreate` fails if `input_25_contour.nc` already exists, so delete it first. |
| `plot_waves.m` | function | Plots direction, Hs and Tp along the contour against time, from `wave_data.mat`. |
| `load_data2.m` | function | Earlier, **unfinished** loader for the CoastSat workbooks. It is not called anywhere and does not run as is (it uses an undefined `input_file` and the loop is not closed). |
| `cmocean.m` | function | Third-party perceptually uniform colormaps (Thyng et al., 2016). |

### Data

| File | Content |
|---|---|
| `nzd0204.xlsx` … `nzd0332.xlsx` | CoastSat exports, one per beach. Sheets: *Transects* (one row per transect: `orientation` in column 3, `land_x`/`land_y` = lon/lat of the landward end in columns 18-19), *Intersect points* (one row per image date, one `"lat,lon"` text cell per transect), *Intersects* (cross-shore distance per transect and date), *Tides*. Site 2 uses 16 workbooks, listed west to east in `initial_grid.m`. |
| `A_site2.mat` | `A{1..16}`: PCA-reconstructed daily cross-shore shoreline anomaly per beach (9284 days x retained transects, from 1999-08-01). |
| `mint_site2.mat` | `mint{1..16}`: mean (offset) per retained transect, added to `A`. |
| `Astd_site2.mat` | `Astd.beach1..16`: standard deviation of the reconstruction (not used by the scripts). |
| `unique_days.mat` | Helper date lists (not used by the scripts). |
| `oceanum_paper2.nc` | Raw Oceanum wave hindcast at 48 nodes: `hs`, `tps`, `dpm`, `botl` (depth), `latitude`, `longitude`, `time` (hours since 1979-02-01). |
| `input_25_contour.nc` | Model wave input made by `create_nc.m`: 94 stations along the 25 m contour, 74,265 3-hourly records. |
| `wave_data.mat` | The same interpolated waves as MATLAB arrays (v7.3). |
| `30_sept_int_w/output.mat` | Stored result of an earlier `hindcast_run` (`O`, `P`, `S`). The benchmark harness takes its initial coastline from here. |

### `ShorelineS_functions/`

A vendored copy of the ShorelineS model (IHE Delft & Deltares, LGPL-2.1).
`ShorelineS(S)` takes the parameter struct `S`, fills in the defaults
(`initialize_defaultvalues.m` documents every parameter), and runs the time
loop. Each step goes through four phases:

1. **Grid**: regrid the coastline (`make_sgrid_mc`), with groynes and structures.
2. **Transport**, per coastline section: interpolate waves to the coast
   (`introduce_wave` → `get_interpolation_on_grid`), refract and shoal them
   (`wave_refraction`, `wave_breakingheight`), find the critical high-angle
   transport angle (`wave_angles` → `get_Sphimax`), compute bulk longshore
   transport (`transport`, KAMP formula here), and handle shadowing
   (`transport_shadow_treat` → `find_shadows_mc`) and boundary conditions.
3. **Coastline change** (`coastline_change`), plus overwash
   (`find_overwash_mc`), merging of sections, and NaN cleanup.
4. **Output**: store to `output.mat` every `S.storageinterval` days
   (`save_shorelines`), and plot every `S.plotinterval` steps (`plot_coast`).

The run settings used in this project (see `hindcast_run.m`):

| Parameter | Value | Meaning |
|---|---|---|
| `reftime` / `endofsimulation` | 2000-01-01 / 2024-12-30 | simulation period |
| `dt`, `tc` | 1/(365·8) yr, 0 | fixed 3-hour time step (≈ 73,000 steps) |
| `wvcfile` | `input_25_contour.nc` | wave time series at 94 contour stations |
| `interpolationmethod` | `alongshoremapping` | how stations map onto the coastline |
| `ddeep`, `dnearshore` | 25 m, 25 m | depth of the wave input and of the nearshore point |
| `ds0` | 75 m | alongshore grid size |
| `trform` | `KAMP` | Kamphuis transport formula |
| `d` | 10 m | active profile height |
| boundary conditions | `periodic` | at both ends |
| `storageinterval` | 30 days | output interval |
| `plotinterval` | 240 steps (30 days) | plot/video interval |

---

## Performance

### How it was measured

There is no MATLAB in the environment used for this work, so all timings come
from **GNU Octave 8.4** (headless, under Xvfb) on a 4-core cloud container.
The harness is `benchmark/bench_run.m`. It runs the exact `hindcast_run.m`
model configuration (same coastline, wave file, grid and physics) for a short
period from 2000-01-01, using the initial coastline stored in
`30_sept_int_w/output.mat`. Absolute times in MATLAB will be lower, because
its JIT runs scalar loops much faster than Octave does, but the hotspots are
the same.

Per-step cost is the difference between a 6-day and a 2-day run (32 steps),
so start-up work (reading the 168 MB NetCDF, first plot) cancels out.

| Run (Octave 8.4, model only, no per-step plotting) | Original | Optimised | Speed-up |
|---|---:|---:|---:|
| 2 simulated days (17 steps), total | 171.0 s | 52.4 s | 3.3× |
| 6 simulated days (49 steps), total | 455.3 s | 104.3 s | 4.4× |
| 30 simulated days (241 steps), total | – | 369.8 s | |
| **Marginal cost per 3-hour time step** | **8.9 s** | **1.4–1.6 s** | **≈ 6×** |
| Extrapolated full 2000–2024 hindcast (73,040 steps) | ≈ 180 h (7.5 days) | ≈ 30 h | ≈ 6× |

Each run was timed alone on the machine. A run uses one core: runtime does
not change with the number of CPUs, because every time step depends on the
previous one and the arrays per step (~600 coastline points) are too small
for multithreading. More cores only help to run several simulations at once.
For the [Python port](python/README.md), see the end of this section.

Plotting adds to that. With the original `hindcast_run.m` setting (plot and
grab a video frame every step), each plotted step cost a further **≈ 8.9 s**
in Octave (2-day run: 306.6 s with plotting vs 164.7 s without). Over a full
run that is about 180 more hours, plus a 73,000-frame video kept in memory.
With `S.plotinterval = 240`, plotting happens about 300 times per run.

| Component (Octave) | Original | Optimised | Speed-up |
|---|---:|---:|---:|
| `get_polydistance`, 609 normals vs. 650-point coastline | 0.78 s/call | 0.006 s/call | 131× |
| ↳ `project_output.m`, full 9,131-day CoastSat series | ≈ 2.1 h (estimated) | 26.5 s (measured) | ≈ 280× |
| `wave_breakingheight`, 5,000 points | 2.93 s | 0.021 s | ≈ 140× |
| `get_intersections`, coastline vs. short line (300 calls) | 0.37 s | 0.11 s | 3.3× |

The 6-day runs with the original and optimised code give **identical
output**: all 57 fields of `O` in `output.mat` are bit-for-bit equal.

### Where the time went (original code)

Profile of the original code over 9 time steps (≈ 112 s in Octave):

| Function | Inclusive time | Why |
|---|---|---|
| `wave_breakingheight` (direct + via `get_Sphimax`) | 35 s | per-point scalar loop with an iterative solve: 240,000 scalar calls to `wave_shoalref`/`get_disper` in 9 steps |
| `plot_coast` | 31 s | full figure redraw (here only one frame; the original setting redraws every step) |
| `introduce_wave` | 19 s | for each of 94 stations, every step: `unique()` over the 74,265-entry time series, `sind`/`cosd` of the whole series, 4 × `interp1` over the whole series |
| `get_intersections` (via `find_shadows_mc`, `find_overwash_mc`) | 19 s | about 26 `repmat` calls per call building every segment-pair matrix, then a scalar loop |

### What was changed

All changes keep the model's behaviour. The equivalence tests in `benchmark/`
compare each function against the original. A 2-day model run gives
**identical coastline positions** (`O.x`, `O.y` and all geometry fields are
bit-for-bit equal). Transport differs only at floating-point round-off level
(`O.QS` max |Δ| = 3·10⁻¹⁰ m³/yr on values of order 10⁶).

| File | Change | Result |
|---|---|---|
| `ShorelineS_functions/wave_breakingheight.m` | The per-point secant iteration is vectorised over all points: points that have converged are masked out, and a vectorised `wave_shoalref_vec` does exactly the same branching as `wave_shoalref`. | ~140× faster for 5,000 points. Results equal to ≤ 4·10⁻¹⁵ (round-off). |
| `ShorelineS_functions/get_intersections.m` | Implicit expansion replaces the ~26 full-size `repmat` copies. Only segments whose bounding box can reach the other polyline are evaluated (NaN-adjacent segments are always kept, so even the original's edge cases stay the same). The `ui`/`uj` loop is vectorised. | Same crossings, indices and coordinates on 3,000 random cases. The fractions `ui`/`uj` occasionally differ by 1 ulp (≤ 9·10⁻¹⁶). ~3× faster for the model's typical coastline-vs-short-line calls. |
| `ShorelineS_functions/introduce_wave.m` | `length(unique(t))>1` becomes `any(t~=t(1))`. Each station is interpolated on the 2 time points around `tnow` (the result is the same for a strictly increasing time axis; otherwise it falls back to the full series). All stations are then interpolated in one `interp1` call per parameter. | Bit-identical output. |
| `ShorelineS_functions/get_polydistance.m` | All normals are intersected with the coastline at once. Only (normal, segment) pairs whose bounding boxes overlap are evaluated, with the same crossing arithmetic. | Bit-identical on the stored model output. ≤ 10⁻¹³ m on random polylines. **131× faster.** This speeds up `project_output.m` (≈ 9,700 calls over the daily CoastSat series) and `save_shorelines`. |
| `ShorelineS_functions/get_Sphimax.m` | The warning state is saved and restored around the 3×3 solves, instead of `warning off`/`warning on` inside the per-point loop. The blanket `warning on` switched on every warning for the rest of the run (in Octave that printed ~500,000 broadcasting notices per 6 simulated days). | Same results; no warning flood. |
| `ShorelineS_functions/make_video.m` | Drops empty frames before `writeVideo`. With `plotinterval>1`, frames are stored at `V(it+1)`, which leaves gaps, and the `try` around `writeVideo` used to swallow the error, so no video was written. | Video works with `plotinterval>1`. |
| `hindcast_run.m` | `S.plotinterval = 240` (one frame per 30 days instead of per 3-hour step). The hard-coded `C:\Users\...` `addpath` becomes a path relative to the script. | Avoids ~73,000 figure redraws, and a video frame array that would need hundreds of GB of memory for the full 25-year run. **This changes the animation**: one frame per 30 days. Set `S.plotinterval=1` for the original behaviour. |
| `create_nc.m` | Interpolation is linear in the node values, so a 94 × 48 weight matrix is built once and each field becomes one matrix product, instead of a `scatteredInterpolant` evaluation for every hour (3 × ~74,000 evaluations) into arrays that grew inside the loop. Hours with missing node values fall back to the original per-hour interpolation. | Same values (≤ 10⁻¹⁵, checked with an equivalent SciPy implementation). Removes the ~220,000 interpolant evaluations and the quadratic array growth. |
| `initial_grid.m` | Only the *Intersect points* columns that are actually used (the transects removed from the PCA) are parsed cell by cell and projected, instead of every cell in the sheet. | For beaches `[2 3]`, parses 159 of 609 columns (26 %). Same output. |
| `project_output.m` | One vectorised `movmean(..., 2, 'omitnan')` replaces the per-transect loop. | Same output. |

### Python port vs. Octave

`python/` contains a Python port of the workflow (see
[python/README.md](python/README.md)). Same machine, each run alone:

| | Octave, original | Octave, optimised | Python |
|---|---:|---:|---:|
| Model, 2 days (17 steps), total | 171.0 s | 52.4 s | 0.9 s |
| Model, 6 days (49 steps), total | 455.3 s | 104.3 s | 2.1 s |
| Model, 30 days (241 steps), total | – | 369.8 s | 9.0 s |
| **Model, cost per time step** | **8.9 s** | **1.4–1.6 s** | **0.035 s** |
| `project_output`, full 9,131-day series | ≈ 2.1 h (est.) | 26.5 s | 25.8 s |
| `initial_grid` (beaches 2, 3) | – ¹ | – ¹ | 4.8 s |
| `create_nc` | – ¹ | – ¹ | 6 s |

¹ Needs `readtable`, `projcrs`/`projfwd` or `scatteredInterpolant`, which
Octave does not have.

FULLRUN_PLACEHOLDER

**Why Python is faster for the model but not for `project_output`.** The
MATLAB model code makes hundreds of thousands of small function calls per
simulated day (`sind`, `interp1`, `unique`, and per-coastline-point calls of
`get_intersections` in `find_shadows_mc` and `find_overwash_mc`). Octave has
no JIT compiler, so each call costs 10–100 µs of interpreter overhead, which
dominates the run time. The Python port does the same arithmetic in batched
NumPy operations: for example, all shadow rays of a coastline section are
tested in one call (276 ms per call in Octave, 3 ms in Python). In
`project_output` both versions are vectorised and spend their time in
compiled array code, so they run at the same speed. MATLAB's JIT makes
function calls and loops much cheaper than Octave's interpreter, so the gap to
MATLAB will be smaller than the gap to Octave (MATLAB was not available to
measure it).

### Reproducing

```matlab
% baseline = the commit before the optimisation
% (shell)  git worktree add ../baseline aa45896
cd benchmark
bench_run('../baseline/ShorelineS_functions', 2, 0, 'out_orig');   % original
bench_run('../ShorelineS_functions',          2, 0, 'out_new');    % optimised
compare_outputs('out_orig/output.mat', 'out_new/output.mat')
test_wave_breakingheight('../baseline/ShorelineS_functions')
test_get_intersections('../baseline/ShorelineS_functions')
test_get_polydistance('../baseline/ShorelineS_functions')
```

In Octave: `pkg load netcdf`, and run headless with
`xvfb-run octave-cli --eval "..."` (the model always draws the first frame).

---

## Notes and possible issues (not changed)

- **Observation window off by one?** Row 1 of `x`/`y` from `initial_grid` is
  1999-08-01. `hindcast_run` sets `dy = days(2000-01-01 − 1999-08-01) = 153`,
  and `project_output` uses rows `dy : dy+n_d`. Row 153 is 1999-12-31, so the
  observations may start one day before the model. If the intent is to start
  on 2000-01-01, the rows would be `dy+1 : dy+1+n_d`.
- **Initial coastline date.** `initial_grid` returns the first row
  (1999-08-01) as the initial coastline, while the model starts on 2000-01-01.
  The `dy` argument is passed in but not used.
- **Stored run.** `30_sept_int_w/output.mat` only has output up to 2004-03
  (52 storage times), so that run probably did not finish the full period.
- `make_video.m` writes `[S.outputdir,'\animation']` with a Windows path
  separator.
- **First CoastSat image skipped.** `initial_grid.m` reads `Range B2:KH517`
  with `readtable`, which uses the first row of the range (row 2, the image
  of 1999-09-29) as variable names. The transects that are not in the PCA are
  therefore held at the position of the *second* image. (The Python port
  confirmed this: it reproduces the stored initial coastline only with this
  behaviour.)
- `create_nc.m` fails if `input_25_contour.nc` already exists (`nccreate`
  does not overwrite).
- The ShorelineS default `S.plotinterval = 1` combined with `S.video = 1`
  keeps a full-resolution frame in memory for every time step. Keep the
  interval coarse for long runs.
