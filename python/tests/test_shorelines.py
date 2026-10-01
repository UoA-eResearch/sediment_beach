"""Tests of the Python port against GNU Octave results.

Run from the python/ folder:   python -m pytest tests
The model regression test needs the Git LFS data files (git lfs pull)."""
import os
import sys

import numpy as np
import pytest
import scipy.io

HERE = os.path.dirname(os.path.abspath(__file__))
PYDIR = os.path.dirname(HERE)
REPO = os.path.dirname(PYDIR)
sys.path.insert(0, PYDIR)

from shorelines.coastline import colon  # noqa: E402
from shorelines.geometry import get_intersections, get_polydistance  # noqa: E402
from shorelines.mfuncs import cosd, interp1, mmod, sind  # noqa: E402

DATA = os.path.join(HERE, "data")


def _lfs_available(path):
    return os.path.exists(path) and os.path.getsize(path) > 1000


# ---------------------------------------------------------------------------
def test_degree_trig_exact_multiples():
    x = np.array([-360.0, -180.0, 0.0, 180.0, 360.0, 540.0])
    assert np.all(sind(x) == 0)
    assert np.all(cosd(x + 90) == 0)
    assert np.array_equal(cosd(np.array([0.0, 360.0])), [1.0, 1.0])


def test_mod_matches_matlab_definition():
    x = np.array([-725.5, -360.0, -0.0, 1e-300, 359.9999, 720.0, 1e6 + 0.25])
    r = mmod(x, 360)
    assert np.all((r >= 0) & (r < 360))
    assert r[1] == 0 and r[5] == 0


def test_interp1_nan_outside_and_sorted():
    x = np.array([3.0, 1.0, 2.0])
    y = np.array([30.0, 10.0, 20.0])
    out = interp1(x, y, np.array([0.5, 1.5, 3.0, 3.5]))
    assert np.isnan(out[0]) and np.isnan(out[3])
    assert out[1] == 15.0 and out[2] == 30.0


def test_colon_length():
    for n in (1, 7, 100, 613):
        r = colon(0.0, 1234.5678 / n, 1234.5678)
        assert r.size == n + 1
        assert r[-1] <= 1234.5678


def test_get_intersections_matches_octave():
    """400 random cases (incl. NaN separators, vertical segments, self-crossings)
    computed with get_intersections.m in GNU Octave."""
    C = scipy.io.loadmat(os.path.join(DATA, "get_intersections_cases.mat"))["C"]
    R = scipy.io.loadmat(os.path.join(DATA, "get_intersections_octave.mat"))["R"]
    for k in range(C.shape[0]):
        xi, yi, xj, yj = [C[k, q].ravel() for q in range(4)]
        out = get_intersections(xi, yi) if C[k, 4].ravel()[0] else get_intersections(xi, yi, xj, yj)
        for q in range(8):
            ref = R[k, q].ravel()
            val = np.asarray(out[q]).ravel()
            assert ref.size == val.size, (k, q)
            np.testing.assert_allclose(val, ref, rtol=0, atol=1e-12, err_msg=f"case {k}, output {q}")


def test_get_polydistance_straight_line():
    xr = np.linspace(0, 1000, 11)
    yr = np.zeros_like(xr)
    xc = np.linspace(-100, 1100, 50)
    yc = np.full_like(xc, -40.0)        # coastline 40 m to the right of the reference line
    d, xcr, ycr = get_polydistance(xr, yr, xc, yc, 500)
    np.testing.assert_allclose(d, -40.0, atol=1e-9)
    np.testing.assert_allclose(ycr, -40.0, atol=1e-9)


# ---------------------------------------------------------------------------
@pytest.mark.skipif(not _lfs_available(os.path.join(REPO, "input_25_contour.nc")),
                    reason="Git LFS data not available")
def test_model_2days_matches_octave(tmp_path):
    """2 simulated days (17 time steps) of the hindcast configuration, compared
    with the output of ShorelineS in GNU Octave 8.4."""
    from shorelines import ShorelineS
    r = scipy.io.loadmat(os.path.join(REPO, "30_sept_int_w", "output.mat"), squeeze_me=True,
                         struct_as_record=False)["S"]
    S = dict(reftime="2000-01-01", endofsimulation="2000-01-03", xmc=r.xmc, ymc=r.ymc,
             wvcfile=os.path.join(REPO, "input_25_contour.nc"), interpolationmethod="alongshoremapping",
             ddeep=25, dnearshore=25, ds0=75, trform="KAMP", boundaryconditionstart="periodic",
             boundaryconditionend="periodic", d=10, dt=1 / (365 * 8), tc=0, storageinterval=30,
             outputdir=str(tmp_path), suppresshighangle=0)
    O = ShorelineS(S, verbose=False, save_output=False).run()
    ref = scipy.io.loadmat(os.path.join(DATA, "octave_2day_output.mat"), squeeze_me=False,
                           struct_as_record=False)["O"][0, 0]
    # coastline geometry: bit-identical
    for f in ("x", "y", "x1", "y1", "distx", "n"):
        np.testing.assert_array_equal(O[f], getattr(ref, f), err_msg=f)
    # waves and transport: equal to round-off
    for f in ("QS", "dSds", "HSo", "PHIo", "TP", "HS", "PHI", "HSbr", "PHIbr", "dPHIbr", "hbr", "PHIc", "PHIf"):
        a = getattr(ref, f)
        np.testing.assert_array_equal(np.isnan(O[f]), np.isnan(a), err_msg=f)
        scale = np.nanmax(np.abs(a))
        np.testing.assert_allclose(O[f], a, rtol=0, atol=1e-13 * scale, err_msg=f)
