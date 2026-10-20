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


# ---------------------------------------------------------------------------
# Cyclic (closed) coastline sections
def test_get_clockpoly_matches_octave():
    """300 random polygons, both orientations, open and closed, computed with
    get_clockpoly.m in GNU Octave."""
    from shorelines.mfuncs import get_clockpoly
    C = scipy.io.loadmat(os.path.join(DATA, "clockpoly_cases.mat"))["C"]
    R = scipy.io.loadmat(os.path.join(DATA, "clockpoly_octave.mat"))["R"].ravel()
    got = [get_clockpoly(C[k, 0].ravel(), C[k, 1].ravel()) for k in range(C.shape[0])]
    np.testing.assert_array_equal(got, R)


def test_overwash_closed_section_matches_octave():
    """find_overwash_mc on a thin closed island whose first point is shadowed,
    so that the first/last-point merge of closed sections is exercised."""
    import shorelines.coastline as cl
    from shorelines.mfuncs import get_cyclic, get_one_polygon, get_one_section
    d = {k: v.ravel() for k, v in scipy.io.loadmat(os.path.join(DATA, "overwash_cyclic_case.mat")).items()
         if not k.startswith("__")}
    C = dict(x_mc=d["x_mc"].copy(), y_mc=d["y_mc"].copy(), n_mc=1, ds0=75.0)

    def seg(ii, with_i2=False):
        x, y, _, i1, i2 = get_one_polygon(C["x_mc"], C["y_mc"], ii)
        _, PHI, _, _, _ = get_one_polygon(d["QS_mc"], d["PHIbr_mc"], ii)
        out = (x, y, get_one_section(d["shadow_mc"], ii), PHI, i1, get_cyclic(x, y, 75.0))
        return out + (i2,) if with_i2 else out

    spit = dict(method="default", spitwidth=50.0, owscale=0.1, owtimescale=0.0, Dsf=8.0, Dbb=4.0, bheight=2.0)
    cl.find_overwash_mc(C, spit, dict(dt=1 / (365 * 8)), seg)
    ref = scipy.io.loadmat(os.path.join(DATA, "overwash_cyclic_octave.mat"))
    np.testing.assert_array_equal(C["x_mc"], ref["x_mc"].ravel())
    np.testing.assert_array_equal(C["y_mc"], ref["y_mc"].ravel())
    assert C["y_mc"][0] == C["y_mc"][-1] != d["y_mc"][0]      # the closed ends moved together


@pytest.mark.skipif(not _lfs_available(os.path.join(REPO, "input_25_contour.nc")),
                    reason="Git LFS data not available")
def test_model_cyclic_island_matches_octave(tmp_path):
    """6 days (49 steps) with the two open beaches plus a thin closed island
    (overwash active), compared with ShorelineS in GNU Octave 8.4, daily output."""
    from shorelines import ShorelineS
    c = scipy.io.loadmat(os.path.join(DATA, "cyclic_coastline_input.mat"))
    S = dict(reftime="2000-01-01", endofsimulation="2000-01-07", xmc=c["xmc"].ravel(), ymc=c["ymc"].ravel(),
             wvcfile=os.path.join(REPO, "input_25_contour.nc"), interpolationmethod="alongshoremapping",
             ddeep=25, dnearshore=25, ds0=75, trform="KAMP", boundaryconditionstart="periodic",
             boundaryconditionend="periodic", d=10, dt=1 / (365 * 8), tc=0, storageinterval=1,
             outputdir=str(tmp_path), suppresshighangle=0)
    model = ShorelineS(S, verbose=False, save_output=False)
    O = model.run()
    assert np.nansum(model.COAST["cyclic_mc"]) == 1                 # one closed section
    ref = scipy.io.loadmat(os.path.join(DATA, "octave_cyclic_6day_output.mat"), squeeze_me=False,
                           struct_as_record=False)["O"][0, 0]
    for f in ("x", "y", "x1", "y1", "distx", "n"):
        np.testing.assert_array_equal(O[f], getattr(ref, f), err_msg=f)
    for f in ("QS", "HSo", "PHIo", "TP", "PHI", "PHIbr", "dPHIbr", "PHIc", "PHIf"):
        a = getattr(ref, f)
        np.testing.assert_array_equal(np.isnan(O[f]), np.isnan(a), err_msg=f)
        np.testing.assert_allclose(O[f], a, rtol=0, atol=1e-10 * np.nanmax(np.abs(a)), err_msg=f)


# ---------------------------------------------------------------------------
def test_model_synthetic_waves_closed_matches_octave(tmp_path):
    """The current hindcast_run settings (start 1999-08-01, closed boundaries,
    synthetic waves of sinth_waves.m, d50 = 1e-4, gamma = 0.5): 2 days compared
    with ShorelineS in GNU Octave 8.4."""
    from shorelines import ShorelineS
    wave_file = os.path.join(REPO, "sinthe_waves.nc")
    if not _lfs_available(wave_file):
        pytest.skip("sinthe_waves.nc not available (git lfs pull)")
    c = scipy.io.loadmat(os.path.join(DATA, "synthetic_closed_input.mat"))
    S = dict(reftime="1999-08-01", endofsimulation="1999-08-03", xmc=c["xmc"].ravel(), ymc=c["ymc"].ravel(),
             wvcfile=wave_file, interpolationmethod="alongshoremapping", ddeep=25, dnearshore=25, ds0=75,
             trform="KAMP", boundaryconditionstart="closed", boundaryconditionend="closed", d=10, d50=1e-4,
             gamma=0.5, dt=1 / (365 * 8), tc=0, storageinterval=1, outputdir=str(tmp_path),
             suppresshighangle=0)
    O = ShorelineS(S, verbose=False, save_output=False).run()
    ref = scipy.io.loadmat(os.path.join(DATA, "octave_synthetic_closed_2day_output.mat"), squeeze_me=False,
                           struct_as_record=False)["O"][0, 0]
    for f in ("x", "y", "x1", "y1", "distx", "n"):
        np.testing.assert_array_equal(O[f], getattr(ref, f), err_msg=f)
    # The synthetic directions include waves from the land side, where the
    # angle is clipped near 90 degrees and dividing by its cosine magnifies
    # last-bit differences of libm (a few points at ~3e-11 relative).
    for f in ("QS", "dSds", "HSo", "PHIo", "TP", "HS", "PHI", "HSbr", "PHIbr", "dPHIbr", "hbr"):
        a = getattr(ref, f)
        np.testing.assert_array_equal(np.isnan(O[f]), np.isnan(a), err_msg=f)
        np.testing.assert_allclose(O[f], a, rtol=1e-9, atol=1e-13 * np.nanmax(np.abs(a)), err_msg=f)


# ---------------------------------------------------------------------------
# Parallel parameter sets (run_ensemble.py)
def test_ensemble_parameter_sets(tmp_path):
    import run_ensemble as re_
    grid = re_.parse_set(["qscal=0.5,1", "boundaryconditionstart=closed,periodic"])
    assert grid == {"qscal": [0.5, 1], "boundaryconditionstart": ["closed", "periodic"]}
    runs = re_.expand([{"d": 8}, {"d": 12}], grid)
    assert len(runs) == 8 and runs[0] == {"d": 8, "qscal": 0.5, "boundaryconditionstart": "closed"}
    assert re_.expand(None, None) == [{}]
    p = tmp_path / "runs.csv"
    p.write_text("qscal,wave_file,end\n0.5,sinthe_waves.nc,2000-01-01\n2,,\n")
    rows = re_.read_runs(str(p))
    assert rows == [{"qscal": 0.5, "wave_file": "sinthe_waves.nc", "end": "2000-01-01"}, {"qscal": 2}]
    run, s_over = re_.split_params(rows[0])
    assert run["end"] == "2000-01-01" and run["wave_file"] == os.path.join(REPO, "sinthe_waves.nc")
    assert s_over == {"qscal": 0.5}
    with pytest.raises(ValueError, match="unknown setting"):
        re_.validate([{"qscall": 1}])
    with pytest.raises(NotImplementedError):
        re_.validate([{"trform": "CERC"}])


def test_ensemble_runs_in_parallel(tmp_path):
    """Two 1-day runs in two processes; each writes its own output and a summary row."""
    import csv
    import run_ensemble as re_
    if not _lfs_available(os.path.join(REPO, "PC_lst_site2.mat")):
        pytest.skip("LFS data not available (git lfs pull)")
    base = dict(re_.BASE, end="1999-08-02", interval=1)
    rows = re_.run_ensemble([{"qscal": 0.5}, {"qscal": 2}], str(tmp_path), jobs=2, base=base, verbose=False)
    assert [r["status"] for r in rows] == ["ok", "ok"], rows
    for r in rows:
        assert (tmp_path / r["out_dir"] / "output.mat").exists()
        assert r["steps"] == 9 and np.isfinite(r["rmse_m"])
    import json
    from hindcast_run import PARAMETERS
    with open(tmp_path / "summary.csv") as f:
        reader = csv.DictReader(f)
        summary = list(reader)
    # every parameter is a column with the value used; the varied one comes first
    assert reader.fieldnames[:3] == ["run", "status", "qscal"]
    assert {name for name, _, _ in PARAMETERS} <= set(reader.fieldnames)
    assert [s["qscal"] for s in summary] == ["0.5", "2"]
    assert [s["gamma"] for s in summary] == ["0.5", "0.5"] and summary[0]["end"] == "1999-08-02"
    for r, q in zip(rows, (0.5, 2)):
        rec = json.loads((tmp_path / r["out_dir"] / "params.json").read_text())
        assert rec["parameter_set"] == {"qscal": q}
        assert rec["parameters"]["qscal"] == q and rec["parameters"]["d50"] == 1e-4
        assert rec["shorelines_settings"]["trform"] == "KAMP"
        S = scipy.io.loadmat(tmp_path / r["out_dir"] / "output.mat", squeeze_me=True, struct_as_record=False)["S"]
        assert S.qscal == q and S.gamma == 0.5 and S.boundaryconditionstart == "closed"
    # larger qscal, more coastline change
    assert rows[1]["model_mean_abs_change_end_m"] > rows[0]["model_mean_abs_change_end_m"]


def test_merge_crossing_sections_matches_octave():
    """merge_coastlines_mc on 400 random pairs of crossing sections (two open
    beaches, a beach and an island in either order, two islands, wavy lines
    crossing many times), compared with ShorelineS in GNU Octave 8.4."""
    import contextlib
    import io
    from shorelines.coastline import merge_coastlines_mc
    cases = scipy.io.loadmat(os.path.join(DATA, "merge_crossing_input.mat"))["cases"]
    ref = scipy.io.loadmat(os.path.join(DATA, "merge_crossing_octave.mat"))
    n_merged = 0
    for k in range(cases.shape[0]):
        x, y = cases[k, 0].ravel(), cases[k, 1].ravel()
        with contextlib.redirect_stdout(io.StringIO()):
            C = merge_coastlines_mc(dict(x_mc=x.copy(), y_mc=y.copy(), ds0=75.0))
        np.testing.assert_array_equal(C["x_mc"], ref["out"][k, 0].ravel(), err_msg=f"case {k}")
        np.testing.assert_array_equal(C["y_mc"], ref["out"][k, 1].ravel(), err_msg=f"case {k}")
        n_merged += C["x_mc"].size != x.size
    assert n_merged > 250

    # the merge of beaches 2 and 3 on 2002-06-10 in the hindcast with the
    # synthetic waves (input recorded in the Python run, output from Octave)
    ev = scipy.io.loadmat(os.path.join(DATA, "merge_event_synthetic_waves.mat"))["ev"]
    ref = scipy.io.loadmat(os.path.join(DATA, "merge_event_synthetic_waves_octave.mat"))["out"]
    for k in range(ev.shape[0]):
        C = merge_coastlines_mc(dict(x_mc=ev[k, 0].ravel().copy(), y_mc=ev[k, 1].ravel().copy(), ds0=75.0))
        np.testing.assert_array_equal(C["x_mc"], ref[k, 0].ravel())
        np.testing.assert_array_equal(C["y_mc"], ref[k, 1].ravel())
        assert np.isnan(C["x_mc"]).sum() < np.isnan(ev[k, 0]).sum()      # one section fewer


def test_help_lists_every_parameter_with_its_default():
    from hindcast_run import PARAMETERS, default_parameters, parameter_help
    text = parameter_help()
    defaults = default_parameters()
    assert defaults["qscal"] == 1.0 and defaults["gamma"] == 0.5 and defaults["d50"] == 1e-4
    assert defaults["boundaryconditionstart"] == "closed" and defaults["wave_file"] == "input_25_contour.nc"
    lines = text.splitlines()
    for name, _, _ in PARAMETERS:
        line = next(ln for ln in lines if ln.split()[:1] == [name])
        assert line.split()[1] == f"{defaults[name]:.6g}" if isinstance(defaults[name], float) \
            else line.split()[1] == str(defaults[name])
