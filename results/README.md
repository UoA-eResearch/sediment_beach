# Results

| File | Content |
|---|---|
| `python_hindcast_2000_2024/output.mat` | Output structure `O` of the full 2000-01-01 → 2024-12-30 hindcast made with the Python port (`python/hindcast_run.py` settings; storage every 30 days, 306 stored times). Same fields as ShorelineS' `output.mat`. Git LFS. |
| `python_hindcast_2000_2024/projection.mat` | `project_output` results for that run: `zg_model` (model cross-shore change, transects × 30-day times), `zg_obs` (CoastSat every 30 days), `zg_obs_d` (CoastSat daily, 9,131 days from 1999-12-31), `sm` (30-day moving mean), `x_initial`/`y_initial`. Positive = seaward. Git LFS. |
| `octave_reference_30days/output.mat` | Output of the (optimised) MATLAB code run in GNU Octave 8.4 for 2000-01-01 → 2000-01-31 (241 steps), the reference the Python port was checked against (coastline bit-identical). Git LFS. |
| `transect_trends.csv` | Linear 2000–2024 trend per transect (model and CoastSat) and the modelled change by 2024. `in_pca = 0` marks transects whose observations are held constant (not in the PCA). |
| `benchmarks.csv` | All runtimes behind the README tables (measured or estimated). |

Figures made from these files: `docs/figures/` (`python python/make_figures.py`).
