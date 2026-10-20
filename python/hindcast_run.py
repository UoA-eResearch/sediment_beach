"""Python port of hindcast_run.m: ShorelineS hindcast for site 2, beaches 2 and 3.

    python hindcast_run.py                       # full run
    python hindcast_run.py --end 2000-03-01      # shorter run
    python hindcast_run.py --out-dir my_run --no-project

To run several parameter sets in parallel, see run_ensemble.py.
"""
import argparse
import datetime as dt
import os
import time

from initial_grid import initial_grid
from project_output import project_output
from shorelines import ShorelineS

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
COASTSAT_START = dt.date(1999, 8, 1)              # first CoastSat date (row 0 of initial_grid's x, y)


def hindcast_settings(x_initial, y_initial, start="1999-08-01", end="2024-12-30",
                      wave_file=os.path.join(REPO, "input_25_contour_av.nc"), interval=30,
                      out_dir="output"):
    """The ShorelineS settings (S) of the hindcast. Override entries to change the run."""
    return dict(
        reftime=start, endofsimulation=end,
        xmc=x_initial, ymc=y_initial,                 # initial coastline (NZTM2000)
        wvcfile=wave_file, interpolationmethod="alongshoremapping",
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


def days_since_coastsat(start, end):
    """(dy, n_d): days from 1 August 1999 to `start`, and days from `start` to `end`."""
    ref = dt.date.fromisoformat(start)
    if ref < COASTSAT_START:
        raise ValueError(f"start {start} is before the first CoastSat date {COASTSAT_START}")
    return (ref - COASTSAT_START).days, (dt.date.fromisoformat(end) - ref).days


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--start", default="1999-08-01")
    ap.add_argument("--end", default="2024-12-30")
    ap.add_argument("--out-dir", default=os.path.join(HERE, "output", "19_oct_closed_av"))
    ap.add_argument("--wave-file", default=os.path.join(REPO, "input_25_contour.nc"),
                    help="wave input (e.g. ../sinthe_waves.nc for the synthetic waves of sinth_waves.m)")
    ap.add_argument("--interval", type=float, default=30, help="storage interval of the output [days]")
    ap.add_argument("--no-project", action="store_true", help="skip project_output")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    dy, n_d = days_since_coastsat(a.start, a.end)

    t0 = time.perf_counter()
    x, y, x_initial, y_initial, g = initial_grid([2, 3], dy, n_d, verbose=not a.quiet)
    t_grid = time.perf_counter() - t0

    S = hindcast_settings(x_initial, y_initial, a.start, a.end, a.wave_file, a.interval, a.out_dir)
    t1 = time.perf_counter()
    model = ShorelineS(S, verbose=not a.quiet)
    O = model.run()
    t_model = time.perf_counter() - t1

    t2 = time.perf_counter()
    if not a.no_project:
        project_output(x, y, a.out_dir, dy, int(a.interval), x_initial, y_initial, n_d, O=O)
    t_proj = time.perf_counter() - t2
    print(f"initial_grid {t_grid:.1f}s | ShorelineS {t_model:.1f}s ({model.TIME['it'] + 1} steps) | "
          f"project_output {t_proj:.1f}s")


if __name__ == "__main__":
    main()
