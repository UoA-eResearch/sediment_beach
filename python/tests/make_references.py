"""Regenerate the GNU Octave reference results in tests/data from the ShorelineS
submodule (../../ShorelineS, or --shorelines DIR):

    python tests/make_references.py              # all model cases + function cases
    python tests/make_references.py F01 bc_A     # only these model cases
    python tests/make_references.py --functions  # only the single-function cases

Needs GNU Octave with the netcdf package (octave-cli, pkg load netcdf) and the Git
LFS input files. The model cases are defined in reference_cases.py."""
import argparse
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import scipy.io

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import reference_cases as rc  # noqa: E402

OCTAVE_DIR = os.path.join(HERE, "octave")


def _octave(expr, cwd):
    cmd = ["octave-cli", "--no-gui", "-q", "--eval", f"addpath('{OCTAVE_DIR}'); {expr}"]
    if not os.environ.get("DISPLAY"):
        cmd = ["xvfb-run", "-a"] + cmd        # ShorelineS opens figures
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)


def _for_octave(v):
    """Numbers as doubles (an int64 would make Octave compute in integers), lists as cells."""
    if isinstance(v, (bool, int, np.integer)):
        return float(v)
    if isinstance(v, list):
        return np.array([_for_octave(e) for e in v], dtype=object)
    return v


def run_model_case(name, fundir):
    with tempfile.TemporaryDirectory() as tmp:
        S = rc.model_settings(name, "out")        # relative: ShorelineS saves to pwd/outputdir
        S = {k: _for_octave(v) for k, v in S.items()}
        infile = os.path.join(tmp, "S.mat")
        scipy.io.savemat(infile, {"S": S})
        outfile = os.path.join(rc.DATA, f"octave_{name}_output.mat")
        r = _octave(f"run_reference_case('{infile}', '{outfile}', '{fundir}')", tmp)
        ok = r.returncode == 0 and os.path.exists(outfile)
        return name, ok, (r.stdout + r.stderr)[-2000:]


def overwash_case_with_phicxy():
    """Add PHIcxy_mc (orientation of the transport grid, as get_coastline_orientation)
    to overwash_cyclic_case.mat."""
    from shorelines.coastline import _transport_points
    from shorelines.mfuncs import atan2d, get_cyclic, mmod
    p = os.path.join(rc.DATA, "overwash_cyclic_case.mat")
    d = {k: v for k, v in scipy.io.loadmat(p).items() if not k.startswith("__")}
    x, y = d["x_mc"].ravel(), d["y_mc"].ravel()
    xq, yq = _transport_points(x, y, get_cyclic(x, y, 75.0))
    d["PHIcxy_mc"] = mmod(360.0 - atan2d(np.diff(yq), np.diff(xq)), 360)[None, :]
    scipy.io.savemat(p, d)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cases", nargs="*", help="model cases (default: all)")
    ap.add_argument("--functions", action="store_true", help="only the single-function cases")
    ap.add_argument("--shorelines", default=os.path.join(rc.REPO, "ShorelineS"),
                    help="ShorelineS folder (default: the submodule)")
    ap.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    a = ap.parse_args()
    fundir = os.path.abspath(os.path.join(a.shorelines, "functions"))
    if not a.cases:
        overwash_case_with_phicxy()
        r = _octave(f"make_function_references('{rc.DATA}', '{fundir}')", HERE)
        print("functions:", "ok" if r.returncode == 0 else "FAILED\n" + r.stdout + r.stderr)
        if a.functions:
            return
    names = a.cases or list(rc.MODEL_CASES)
    with ThreadPoolExecutor(a.jobs) as ex:
        for name, ok, log in ex.map(lambda n: run_model_case(n, fundir), names):
            print(f"{name}: {'ok' if ok else 'FAILED'}")
            if not ok:
                print(log)


if __name__ == "__main__":
    main()
