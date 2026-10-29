"""Python port of hindcast_run.m: ShorelineS hindcast for site 2, beaches 2 and 3.

    python hindcast_run.py                       # full run
    python hindcast_run.py --end 2000-03-01      # shorter run
    python hindcast_run.py --out-dir my_run --no-project
    python hindcast_run.py --set qscal=0.5 --set gamma=0.72   # change model parameters

The output folder gets output.mat (with the settings used as S, like
ShorelineS) and settings.json (every parameter value used, readable).
To run several parameter sets in parallel, see run_ensemble.py.
"""
import argparse
import datetime as dt
import json
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
COASTSAT_START = dt.date(1999, 8, 1)              # first CoastSat date (row 0 of initial_grid's x, y)

# Defaults of the run options (not ShorelineS settings)
# (wave_file: relative to the working directory, else to the repository)
RUN_DEFAULTS = dict(start="1999-08-01", end="2024-12-30", wave_file="input_25_contour_av.nc", interval=30)
DEFAULT_OUT_DIR = os.path.join(HERE, "output", "27_oct_closed_trselect")

# Parameters that can be changed (--set here, --set/--runs in run_ensemble.py):
# (name, unit, description). The defaults are taken from the code, see default_parameters().
PARAMETERS = [
    # run options
    ("start", "date", "start of the simulation; the initial coastline is the CoastSat/PCA coastline of this date"),
    ("end", "date", "end of the simulation"),
    ("wave_file", "", "wave input (NetCDF stations), e.g. ../sinthe_waves.nc for the synthetic waves"),
    ("interval", "days", "storage interval of the output"),
    # transport
    ("trform", "", "transport formula: KAMP (Kamphuis), CERC (S.b), CERC2, CERC3, MILH (Mil-Homens) or "
     "VR14 (Van Rijn 2014, S.pswell)"),
    ("qscal", "-", "calibration factor of the transport"),
    ("b", "m3/yr", "CERC coefficient (trform CERC)"),
    ("pswell", "%", "percentage swell (trform VR14)"),
    ("d50", "m", "median grain diameter"),
    ("gamma", "-", "breaking coefficient Hs/h"),
    ("porosity", "-", "porosity of the sediment"),
    ("tanbeta", "-", "mean bed slope"),
    ("rhos", "kg/m3", "density of sand"),
    ("rhow", "kg/m3", "density of water"),
    ("g", "m/s2", "gravitational acceleration"),
    ("suppresshighangle", "0/1", "limit the transport angle to the critical high angle"),
    ("twopoints", "0/1", "high-angle correction spread over two downdrift cells"),
    ("maxangle", "deg", "maximum coastline re-orientation between cells (smoothing of sharp angles)"),
    ("relaxationlength", "m", "length over which the transport adapts to a decrease downdrift "
     "(none: no relaxation)"),
    ("sphimax", "", "critical (maximum-transport) wave angle: none = computed every time step, "
     "estimate = computed once"),
    # coastline and grid
    ("ds0", "m", "alongshore grid size"),
    ("d", "m", "active profile height: one value, or [[x, y, d], ...] interpolated along the coast "
     "(give that in a --runs JSON file)"),
    ("phif", "deg", "foreshore orientation: none = smoothed initial coastline, or one angle, or "
     "[[x, y, angle], ...]"),
    ("preserveorientation", "0/1", "move the coastline normal to the initial grid (fails in ShorelineS "
     "for coasts of several sections)"),
    ("boundaryconditionstart", "", "start of open sections: closed, periodic, neumann, fixed (= neumann), "
     "angleconstant or gradient (= angleconstant); closed:Q fixes the transport at Q m3/yr, "
     "angleconstant:A the coastline angle at A degrees"),
    ("boundaryconditionend", "", "end of open sections: as boundaryconditionstart"),
    # waves
    ("ddeep", "m", "water depth of the wave input"),
    ("dnearshore", "m", "depth of the nearshore point (refraction/shoaling down to here)"),
    ("interpolationmethod", "", "wave stations to coast: alongshoremapping or weighted_distance"),
    ("smoothrefrac", "0-1", "smoothing of the coastline angle used only for the refraction from the "
     "nearshore point to breaking (not for the transport)"),
    # time
    ("dt", "yr", "fixed time step (1/(365*8) = 3 hours)"),
    # overwash of narrow spits
    ("spitwidth", "m", "critical spit width below which overwash occurs"),
    ("owscale", "-", "fraction of the overwash deficit moved per time step"),
    ("owtimescale", "yr", "overwash time scale (used instead of owscale when > 0)"),
    ("bheight", "m", "berm height used for overwash"),
    ("spitdsf", "m", "shoreface depth used for the spit width"),
    ("spitdbb", "m", "back-barrier depth used for the spit width"),
]
RUN_KEYS = set(RUN_DEFAULTS)


def hindcast_settings(x_initial, y_initial, start=RUN_DEFAULTS["start"], end=RUN_DEFAULTS["end"],
                      wave_file=RUN_DEFAULTS["wave_file"], interval=RUN_DEFAULTS["interval"], out_dir=DEFAULT_OUT_DIR):
    """The ShorelineS settings (S) of the hindcast. Override entries to change the run."""
    return dict(
        reftime=start, endofsimulation=end,
        xmc=x_initial, ymc=y_initial,                 # initial coastline (NZTM2000)
        wvcfile=resolve_path(str(wave_file)), interpolationmethod="alongshoremapping",
        ddeep=25, dnearshore=25,
        ds0=75,                                       # initial space step [m]
        trform="KAMP",
        boundaryconditionstart="closed", boundaryconditionend="closed",
        d=10,                                         # active profile height [m]
        d50=1.0e-4,                                   # median grain size [m]
        dt=1 / (365 * 8), tc=0,                       # fixed 3-hourly time step [yr]
        storageinterval=interval, outputdir=out_dir,
        suppresshighangle=0,
    )


def full_settings(run, overrides=None, x_initial=None, y_initial=None, out_dir="output"):
    """Every ShorelineS setting of a run: the ShorelineS defaults, then the hindcast
    settings, then `overrides` (setting names in lower case, without underscores)."""
    from shorelines.model import DEFAULTS, _normalise_keys
    S = dict(DEFAULTS)
    S.update(_normalise_keys(hindcast_settings(x_initial, y_initial, run["start"], run["end"],
                                               run["wave_file"], run["interval"], out_dir)))
    S.update(overrides or {})
    return S


def parameter_values(run, S):
    """{name: value} of PARAMETERS (run options and settings) for one run."""
    vals = {}
    for name, _, _ in PARAMETERS:
        vals[name] = run[name] if name in RUN_KEYS else S[name]
    vals["wave_file"] = display_path(resolve_path(str(vals["wave_file"])))
    return vals


def default_parameters():
    """Default value of every entry of PARAMETERS."""
    run = dict(RUN_DEFAULTS)
    return parameter_values(run, full_settings(run))


def display_path(path):
    """Path relative to the repository when inside it (as shown in help texts and summaries)."""
    path = os.path.abspath(str(path))
    rel = os.path.relpath(path, REPO)
    return path if rel.startswith("..") else rel


def _fmt(v):
    if v is None:
        return "none"
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def parameter_help():
    """The parameter table for --help texts: name, default value, unit and meaning."""
    defaults = default_parameters()
    lines = ["parameters (name, default, unit, meaning):"]
    for name, unit, text in PARAMETERS:
        lines.append(f"  {name:23s} {_fmt(defaults[name]):19s} {unit:6s} {text}")
    lines.append("  fixed in this configuration: tc = 0, griddingmethod = 2; any other\n"
                 "  ShorelineS setting the Python port supports can also be given (shorelines/model.py, DEFAULTS).")
    return "\n".join(lines)


def to_json(v):
    """JSON-compatible value (numpy scalars and arrays included)."""
    if isinstance(v, np.ndarray):
        return v.tolist() if v.size <= 20 else f"<array of {v.size} values, see output.mat S>"
    if isinstance(v, (np.floating, np.integer, np.bool_)):
        return v.item()
    if isinstance(v, float) and not np.isfinite(v):
        return str(v)
    return v


def write_settings(path, run, S, parameter_set=None):
    """settings.json/params.json: the parameter values used (PARAMETERS) and every ShorelineS setting."""
    rec = {}
    if parameter_set is not None:
        rec["parameter_set"] = {k: to_json(v) for k, v in parameter_set.items()}
    rec["parameters"] = {k: to_json(v) for k, v in parameter_values(run, S).items()}
    rec["shorelines_settings"] = {k: to_json(v) for k, v in sorted(S.items())}
    with open(path, "w") as f:
        json.dump(rec, f, indent=2)


def parse_assignments(items):
    """['qscal=0.5', 'boundaryconditionstart=periodic'] -> {'qscal': 0.5, ...}."""
    out = {}
    for item in items or []:
        if "=" not in item:
            raise ValueError(f"--set expects key=value, got {item!r}")
        k, v = item.split("=", 1)
        try:
            v = json.loads(v.strip())
        except ValueError:
            v = v.strip()
        out[k.strip()] = v
    return out


def split_params(params, base=None):
    """(run options, ShorelineS overrides) of one parameter set."""
    run = dict(RUN_DEFAULTS if base is None else base)
    s_over = {}
    for k, v in params.items():
        key = k.replace("-", "_").lower()
        if key in RUN_KEYS:
            run[key] = v
        else:
            s_over[key.replace("_", "")] = v
    run["wave_file"] = resolve_path(str(run["wave_file"]))
    run["interval"] = float(run["interval"])
    return run, s_over


def resolve_path(path):
    """A path as given, relative to the working directory, or relative to the repository."""
    if os.path.isabs(path) or os.path.exists(path):
        return os.path.abspath(path)
    in_repo = os.path.join(REPO, path)
    return in_repo if os.path.exists(in_repo) else os.path.abspath(path)


def check_parameters(run, s_over, what="parameters"):
    """Raise ValueError/FileNotFoundError/NotImplementedError for an unusable parameter set."""
    from shorelines.model import DEFAULTS, _check_supported
    unknown = sorted(k for k in s_over if k not in DEFAULTS)
    if unknown:
        raise ValueError(f"{what}: unknown setting(s) {unknown}")
    days_since_coastsat(run["start"], run["end"])
    if not os.path.exists(run["wave_file"]):
        raise FileNotFoundError(f"{what}: wave file {run['wave_file']} not found")
    _check_supported(full_settings(run, s_over))


def days_since_coastsat(start, end):
    """(dy, n_d): days from 1 August 1999 to `start`, and days from `start` to `end`."""
    ref = dt.date.fromisoformat(start)
    if ref < COASTSAT_START:
        raise ValueError(f"start {start} is before the first CoastSat date {COASTSAT_START}")
    return (ref - COASTSAT_START).days, (dt.date.fromisoformat(end) - ref).days


class HelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Keeps the layout of the description/epilog and adds '(default: ...)' to the options that take a value."""

    def _get_help_string(self, action):
        h = action.help or ""
        if "default" in h or action.default in (None, False, [], argparse.SUPPRESS):
            return h
        return h + " (default: %(default)s)"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, epilog=parameter_help(), formatter_class=HelpFormatter)
    ap.add_argument("--start", default=RUN_DEFAULTS["start"], help="start date")
    ap.add_argument("--end", default=RUN_DEFAULTS["end"], help="end date")
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR, help="output folder (default: python/output/27_oct_closed_trselect)")
    ap.add_argument("--wave-file", default=RUN_DEFAULTS["wave_file"],
                    help="wave input (e.g. ../sinthe_waves.nc for the synthetic waves of sinth_waves.m)")
    ap.add_argument("--interval", type=float, default=RUN_DEFAULTS["interval"],
                    help="storage interval of the output [days]")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help="change a parameter (see the list below); repeat for several")
    ap.add_argument("--no-project", action="store_true", help="skip project_output")
    ap.add_argument("--quiet", action="store_true", help="no progress output")
    a = ap.parse_args(argv)

    base = dict(start=a.start, end=a.end, wave_file=a.wave_file, interval=a.interval)
    try:
        run, s_over = split_params(parse_assignments(a.set), base)
        check_parameters(run, s_over, "parameters")
    except (ValueError, FileNotFoundError, NotImplementedError) as e:
        ap.error(str(e))

    from initial_grid import initial_grid
    from project_output import project_output
    from shorelines import ShorelineS

    dy, n_d = days_since_coastsat(run["start"], run["end"])
    t0 = time.perf_counter()
    x, y, x_initial, y_initial, g = initial_grid([2, 3], dy, n_d, verbose=not a.quiet)
    t_grid = time.perf_counter() - t0

    S = hindcast_settings(x_initial, y_initial, run["start"], run["end"], run["wave_file"], run["interval"],
                          a.out_dir)
    S.update(s_over)
    t1 = time.perf_counter()
    model = ShorelineS(S, verbose=not a.quiet)
    os.makedirs(a.out_dir, exist_ok=True)
    write_settings(os.path.join(a.out_dir, "settings.json"), run, model.S)
    O = model.run()
    t_model = time.perf_counter() - t1

    t2 = time.perf_counter()
    if not a.no_project:
        project_output(x, y, a.out_dir, dy, int(run["interval"]), x_initial, y_initial, n_d, O=O)
    t_proj = time.perf_counter() - t2
    print(f"initial_grid {t_grid:.1f}s | ShorelineS {t_model:.1f}s ({model.TIME['it'] + 1} steps) | "
          f"project_output {t_proj:.1f}s")


if __name__ == "__main__":
    main()
