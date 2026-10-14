"""Run the hindcast for several parameter sets in parallel (one process per run).

Each parameter set overrides entries of the hindcast settings (hindcast_run.py,
MATLAB field names, case-insensitive) and runs in its own process, writing to
its own folder. The time loop of one run is sequential, so running parameter
sets side by side is how this workload uses several CPUs.

    # 3 x 2 grid (Cartesian product), 4 processes
    python run_ensemble.py --set qscal=0.5,1,2 --set gamma=0.5,0.72 --end 2001-08-01 -j 4

    # parameter sets listed in a file (JSON list of objects, or CSV with a header)
    python run_ensemble.py --runs runs.json --out-dir output/ensemble1

    # both: every row of the file combined with every --set value
    python run_ensemble.py --runs runs.csv --set wave_file=../sinthe_waves.nc,../input_25_contour.nc

Besides the ShorelineS settings (qscal, d, d50, gamma, ds0, dnearshore,
boundaryconditionstart, ...), a parameter set may contain
start, end (dates), wave_file (wave input .nc) and interval (storage interval, days).

Output (in --out-dir):
    run_000_<params>/output.mat, params.json   one folder per run
    run_000_<params>/project_output.png        with --plot
    summary.csv                                one row per run: parameters,
                                               status, run time and the
                                               model-vs-CoastSat statistics
"""
import os

# One BLAS/OpenMP thread per process: the runs are the parallel unit, and the
# per-step arrays are too small for multithreaded BLAS to help.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse                                     # noqa: E402
import csv                                          # noqa: E402
import itertools                                    # noqa: E402
import json                                         # noqa: E402
import re                                           # noqa: E402
import shutil                                       # noqa: E402
import time                                         # noqa: E402
import traceback                                    # noqa: E402
from concurrent.futures import ProcessPoolExecutor, as_completed   # noqa: E402

import numpy as np                                  # noqa: E402

from hindcast_run import HERE, REPO, days_since_coastsat, hindcast_settings   # noqa: E402

# Parameters handled by the runner itself (not ShorelineS settings)
RUN_KEYS = {"start", "end", "wave_file", "interval"}
BASE = dict(start="1999-08-01", end="2024-12-30",
            wave_file=os.path.join(REPO, "input_25_contour.nc"), interval=30)


# ----------------------------------------------------------------------------
# Parameter sets
def _value(text):
    """'0.5' -> 0.5, 'closed' -> 'closed', 'true' -> True, '[1,2]' -> [1, 2]."""
    try:
        return json.loads(text)
    except ValueError:
        return text


def parse_set(items):
    """['qscal=0.5,1', 'gamma=0.72'] -> {'qscal': [0.5, 1], 'gamma': [0.72]}."""
    grid = {}
    for item in items or []:
        if "=" not in item:
            raise ValueError(f"--set expects key=value[,value...], got {item!r}")
        key, vals = item.split("=", 1)
        grid[key.strip()] = [_value(v.strip()) for v in vals.split(",")]
    return grid


def read_runs(path):
    """Parameter sets from a JSON list of objects or a CSV file with a header row."""
    if path.lower().endswith(".json"):
        with open(path) as f:
            runs = json.load(f)
        if isinstance(runs, dict):
            runs = [runs]
    else:
        with open(path, newline="") as f:
            runs = [{k.strip(): _value(v.strip()) for k, v in row.items() if v is not None and v.strip() != ""}
                    for row in csv.DictReader(f)]
    if not all(isinstance(r, dict) for r in runs):
        raise ValueError(f"{path}: expected a list of parameter sets (objects)")
    return runs


def expand(runs=None, grid=None):
    """Every run in `runs` (default: one empty run) combined with every point of `grid`."""
    runs = runs or [{}]
    grid = grid or {}
    keys = list(grid)
    out = []
    for r in runs:
        for combo in itertools.product(*(grid[k] for k in keys)):
            p = dict(r)
            p.update(zip(keys, combo))
            out.append(p)
    return out


def _resolve(path):
    if os.path.isabs(path) or os.path.exists(path):
        return os.path.abspath(path)
    return os.path.join(REPO, path) if os.path.exists(os.path.join(REPO, path)) else os.path.abspath(path)


def split_params(params, base=None):
    """(run options, ShorelineS overrides) of one parameter set."""
    run = dict(BASE if base is None else base)
    s_over = {}
    for k, v in params.items():
        key = k.replace("-", "_").lower()
        if key in RUN_KEYS:
            run[key] = v
        else:
            s_over[key.replace("_", "")] = v
    run["wave_file"] = _resolve(str(run["wave_file"]))
    run["interval"] = float(run["interval"])
    return run, s_over


def label(params):
    """Short folder-name label for a parameter set."""
    parts = []
    for k, v in params.items():
        if k in ("wave_file",):
            v = os.path.splitext(os.path.basename(str(v)))[0]
        parts.append(f"{k}={v}")
    return re.sub(r"[^A-Za-z0-9=.+-]+", "_", "_".join(parts))[:80].strip("_")


def validate(param_sets, base=None):
    """Check every parameter set before any run starts (unknown keys, unported options, files)."""
    from shorelines.model import DEFAULTS, _check_supported, _normalise_keys
    for i, p in enumerate(param_sets):
        run, s_over = split_params(p, base)
        unknown = sorted(k for k in s_over if k not in DEFAULTS)
        if unknown:
            raise ValueError(f"parameter set {i}: unknown setting(s) {unknown}")
        days_since_coastsat(run["start"], run["end"])
        if not os.path.exists(run["wave_file"]):
            raise FileNotFoundError(f"parameter set {i}: wave file {run['wave_file']} not found")
        S = dict(DEFAULTS)
        S.update(_normalise_keys(hindcast_settings(None, None, run["start"], run["end"], run["wave_file"],
                                                   run["interval"])))
        S.update(s_over)
        _check_supported(S)


# ----------------------------------------------------------------------------
# Worker
_GRID = {}


def _init_worker(grid_dir):
    """Load the observed coastlines (from initial_grid, computed once) in each worker."""
    _GRID["x"] = np.load(os.path.join(grid_dir, "x.npy"), mmap_mode="r")
    _GRID["y"] = np.load(os.path.join(grid_dir, "y.npy"), mmap_mode="r")


def _stats(P):
    """Model-vs-CoastSat statistics of project_output's cross-shore distances."""
    zm, zo = P["zg_model"], P["zg_obs"]
    ok = ~np.isnan(zm) & ~np.isnan(zo)
    out = dict(n_compared=int(ok.sum()))
    if ok.sum() > 1:
        d = zm[ok] - zo[ok]
        out.update(rmse_m=float(np.sqrt(np.mean(d ** 2))), bias_m=float(np.mean(d)),
                   corr=float(np.corrcoef(zm[ok], zo[ok])[0, 1]))
    end = zm[:, -1] if zm.shape[1] else np.array([np.nan])
    out["model_mean_abs_change_end_m"] = float(np.nanmean(np.abs(end))) if np.any(~np.isnan(end)) else np.nan
    return out


def run_one(i, params, out_root, base=None, project=True, plot=False):
    """Run parameter set `params` (the i-th) in folder out_root/run_<i>_<label>. Returns a summary row."""
    from project_output import project_output
    from shorelines import ShorelineS

    run, s_over = split_params(params, base)
    out_dir = os.path.join(out_root, f"run_{i:03d}" + (f"_{label(params)}" if params else ""))
    os.makedirs(out_dir, exist_ok=True)
    row = dict(run=i, out_dir=os.path.relpath(out_dir, out_root), status="ok", error="")
    t0 = time.perf_counter()
    try:
        with open(os.path.join(out_dir, "params.json"), "w") as f:
            json.dump(dict(params=params, run=run, shorelines=s_over), f, indent=2, default=str)
        dy, n_d = days_since_coastsat(run["start"], run["end"])
        x, y = _GRID["x"], _GRID["y"]
        x_initial, y_initial = np.array(x[dy]), np.array(y[dy])    # initial_grid: coastline on the start date
        S = hindcast_settings(x_initial, y_initial, run["start"], run["end"], run["wave_file"],
                              run["interval"], out_dir)
        S.update(s_over)
        model = ShorelineS(S, verbose=False)
        O = model.run()
        row["steps"] = model.TIME["it"] + 1
        row["model_s"] = round(time.perf_counter() - t0, 2)
        if project:
            P = project_output(np.asarray(x), np.asarray(y), out_dir, dy, int(run["interval"]),
                               x_initial, y_initial, n_d, O=O, plot=plot)
            row.update(_stats(P))
    except Exception as e:                           # one failed run must not stop the others
        row.update(status="failed", error=f"{type(e).__name__}: {e}")
        with open(os.path.join(out_dir, "error.txt"), "w") as f:
            f.write(traceback.format_exc())
    row["wall_s"] = round(time.perf_counter() - t0, 2)
    return row


# ----------------------------------------------------------------------------
def run_ensemble(param_sets, out_dir, jobs=None, base=None, project=True, plot=False, verbose=True):
    """Run every parameter set (list of dicts) in parallel. Returns the summary rows (also summary.csv)."""
    from initial_grid import initial_grid

    param_sets = [dict(p) for p in param_sets]
    validate(param_sets, base)
    os.makedirs(out_dir, exist_ok=True)
    jobs = max(1, min(jobs or os.cpu_count() or 1, len(param_sets)))

    # Observed coastlines: computed once (dy = 0 keeps the full series), shared via .npy files
    grid_dir = os.path.join(out_dir, "_grid")
    os.makedirs(grid_dir, exist_ok=True)
    x, y, _, _, _ = initial_grid([2, 3], 0, None, verbose=False)
    np.save(os.path.join(grid_dir, "x.npy"), x)
    np.save(os.path.join(grid_dir, "y.npy"), y)
    with open(os.path.join(out_dir, "ensemble.json"), "w") as f:
        json.dump(dict(base=base or BASE, param_sets=param_sets), f, indent=2, default=str)

    if verbose:
        print(f"{len(param_sets)} runs on {jobs} processes -> {out_dir}")
    t0 = time.perf_counter()
    rows = []
    try:
        with ProcessPoolExecutor(max_workers=jobs, initializer=_init_worker, initargs=(grid_dir,)) as ex:
            futs = {ex.submit(run_one, i, p, out_dir, base, project, plot): i for i, p in enumerate(param_sets)}
            for fut in as_completed(futs):
                row = fut.result()
                rows.append(row)
                if verbose:
                    extra = f"rmse {row['rmse_m']:.2f} m" if "rmse_m" in row else row["error"]
                    print(f"  [{len(rows)}/{len(param_sets)}] run {row['run']:3d} {row['status']:6s} "
                          f"{row['wall_s']:8.1f} s  {extra}")
    finally:
        shutil.rmtree(grid_dir, ignore_errors=True)
    rows.sort(key=lambda r: r["run"])
    for r in rows:
        r.update({f"param_{k}": v for k, v in param_sets[r["run"]].items()})
    write_summary(rows, os.path.join(out_dir, "summary.csv"))
    if verbose:
        n_fail = sum(r["status"] != "ok" for r in rows)
        print(f"done in {time.perf_counter() - t0:.1f} s ({n_fail} failed); summary: "
              f"{os.path.join(out_dir, 'summary.csv')}")
    return rows


def write_summary(rows, path):
    first = ["run", "status"]
    params = sorted({k for r in rows for k in r if k.startswith("param_")})
    rest = sorted({k for r in rows for k in r} - set(first) - set(params) - {"error", "out_dir"})
    fields = first + params + rest + ["out_dir", "error"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set", action="append", metavar="KEY=V1,V2,...",
                    help="parameter values; several --set options form a grid (Cartesian product)")
    ap.add_argument("--runs", help="JSON (list of objects) or CSV file with one parameter set per entry/row")
    ap.add_argument("--start", default=BASE["start"], help="default start date of every run")
    ap.add_argument("--end", default=BASE["end"], help="default end date of every run")
    ap.add_argument("--wave-file", default=BASE["wave_file"], help="default wave input")
    ap.add_argument("--interval", type=float, default=BASE["interval"], help="storage interval [days]")
    ap.add_argument("--out-dir", default=os.path.join(HERE, "output", "ensemble"))
    ap.add_argument("-j", "--jobs", type=int, default=None, help="number of processes (default: number of CPUs)")
    ap.add_argument("--no-project", action="store_true", help="skip project_output and the statistics")
    ap.add_argument("--plot", action="store_true", help="save project_output.png for every run")
    a = ap.parse_args(argv)

    base = dict(start=a.start, end=a.end, wave_file=a.wave_file, interval=a.interval)
    try:
        param_sets = expand(read_runs(a.runs) if a.runs else None, parse_set(a.set))
        validate(param_sets, base)
    except (ValueError, FileNotFoundError, NotImplementedError) as e:
        ap.error(str(e))
    rows = run_ensemble(param_sets, a.out_dir, a.jobs, base, project=not a.no_project, plot=a.plot)
    return 0 if all(r["status"] == "ok" for r in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
