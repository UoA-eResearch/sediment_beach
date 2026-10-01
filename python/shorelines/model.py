"""ShorelineS time loop (port of ShorelineS.m) for the configuration used in
this project: one or more open coastline sections, wave time series at
stations along a depth contour, Kamphuis transport, fixed time step, no
structures, dunes, mud, nourishments or diffraction.

Options outside that scope raise ``NotImplementedError`` instead of being
silently ignored."""
import os
import time

import numpy as np
import scipy.io

from . import coastline as cl
from .geometry import find_shadows_mc
from .mfuncs import (atan2d, cosd, datenum, datestr, get_cyclic, get_one_polygon, m_pow,
                     get_one_section, mmod, sind)
from .transport import (get_smoothangles, get_upwindcorrection, transport,
                        transport_boundary_condition, transport_shadow_treat)
from .waves import (introduce_wave, read_wave_stations, wave_angles,
                    wave_breakingheight, wave_refraction)

# Defaults of initialize_defaultvalues.m that matter for this configuration
DEFAULTS = dict(
    ddeep=25.0, dnearshore=8.0, interpolationmethod="weighted_distance",
    xmc=None, ymc=None, ds0=100.0, griddingmethod=2, d=10.0, phif=None,
    maxangle=60.0, preserveorientation=0, trform="CERC", qscal=1.0,
    d50=2.0e-4, porosity=0.4, tanbeta=0.03, rhos=2650.0, rhow=1025.0, g=9.81,
    gamma=0.72, relaxationlength=None, suppresshighangle=0, tc=1, dt=0.0,
    reftime="2020-01-01", endofsimulation="2040-01-01", twopoints=1,
    smoothfac=0.0, smoothrefrac=0.0, boundaryconditionstart="Fixed",
    boundaryconditionend="Fixed", wvcfile="", sphimax=None,
    spitmethod="default", spitwidth=50.0, spitheadwidth=200.0, owscale=0.1,
    owtimescale=0.0, bheight=2.0, storageinterval=50.0, outputdir="Output",
    # features that are not ported (must stay off)
    xhard=None, ldbstructures="", xrevet=None, ldbrevetments="", xperm=None,
    perm=0, dune=0, mud=0, nourish=0, fnourish=0, diffraction=0, channel=0,
    tideinteraction=False, waveinteraction=False, xsedlim=None, ldbsedlim="",
    ccslr=None, cchs=None, ccdir=None, ldbcoastline="", transmission=0,
    xyout=None, xyprofiles=None,
)
# spitdsf/spitdbb are computed from the *default* S.d in ShorelineS
DEFAULTS["spitdsf"] = DEFAULTS["d"] * 0.8
DEFAULTS["spitdbb"] = 0.5 * DEFAULTS["spitdsf"]


def _normalise_keys(S0):
    return {k.replace("_", "").lower(): v for k, v in S0.items()}


def _check_supported(S):
    def on(v):
        return v not in (None, "", 0, False) and not (isinstance(v, (list, np.ndarray)) and np.size(v) == 0)
    for key in ("xhard", "ldbstructures", "xrevet", "ldbrevetments", "xperm", "perm", "dune", "mud",
                "nourish", "fnourish", "diffraction", "channel", "tideinteraction", "waveinteraction",
                "xsedlim", "ldbsedlim", "ccslr", "cchs", "ccdir", "ldbcoastline", "transmission",
                "xyout", "xyprofiles", "phif", "relaxationlength", "sphimax"):
        if on(S.get(key)):
            raise NotImplementedError(f"S.{key} is not supported by the Python port")
    if S["tc"] != 0 or S["dt"] <= 0:
        raise NotImplementedError("only a fixed time step (S.tc = 0, S.dt > 0) is ported")
    if not S["wvcfile"] or not str(S["wvcfile"]).endswith(".nc"):
        raise NotImplementedError("wave input must be a NetCDF station file (S.wvcfile)")
    if S["trform"].upper() != "KAMP":
        raise NotImplementedError("only the KAMP transport formula is ported")
    if S["griddingmethod"] != 2:
        raise NotImplementedError("only griddingmethod 2 is ported")


def _add_vector(a, b):
    """addVECTOR of save_shorelines.m: append column b, padding with NaN."""
    b = np.asarray(b, dtype=float).ravel()[:, None]
    if a is None or a.size == 0:
        return b.copy()
    if b.shape[0] > a.shape[0]:
        a = np.vstack([a, np.full((b.shape[0] - a.shape[0], a.shape[1]), np.nan)])
    elif a.shape[0] > b.shape[0]:
        b = np.vstack([b, np.full((a.shape[0] - b.shape[0], 1), np.nan)])
    return np.hstack([a, b])


_COAST_COLLECT = ["s", "ds", "x1", "y1", "xq", "yq", "xq1", "yq1", "PHIc", "dPHIc", "PHIcxy",
                  "PHIcxy0", "PHIf", "h0", "cyclic", "clockwise"]
_WAVE_COLLECT = ["PHIo", "PHItdp", "PHIbr", "HSo", "HStdp", "HSbr", "dPHIo", "dPHItdp", "dPHIbr",
                 "TP", "hbr", "dPHIcrit", "diff"]
_TRANSP_COLLECT = ["QS", "QSmax", "shadowS", "shadowS_h", "shadowS_hD", "shadow", "shadow_h",
                   "ivals", "idrev"]
_O_FIELDS = ["it", "dt", "tc", "nt", "timenum", "n", "x", "y", "distx", "distQS", "dSds",
             "wberm", "xdune", "ydune", "qs", "ql", "qw", "R", "SWL", "PHIcxy", "Bf", "Bm", "Bfm",
             "n1", "x1", "y1", "distx1", "distQS1", "PHIc", "QS", "HSo", "PHIo", "TP", "PHIf",
             "HS", "PHI", "dPHI", "HSbr", "PHIbr", "dPHIbr", "hbr", "xhard", "yhard", "nhard",
             "xnour", "ynour", "nnour", "x_fnour", "y_fnour", "n_fnour", "V_fnour_t", "q_fnour_t",
             "x_groyne", "y_groyne"]


def _collect(store, src, fields, i_mc):
    for f in fields:
        v = np.atleast_1d(np.asarray(src[f], dtype=float))
        if i_mc == 1:
            store[f + "_mc"] = v.copy()
        else:
            store[f + "_mc"] = np.concatenate([store[f + "_mc"], [np.nan], v])


class ShorelineS:
    """One model run. ``run()`` executes the time loop and returns the output
    structure ``O`` (dict of arrays, same fields as ShorelineS' output.mat)."""

    def __init__(self, S0, verbose=True, save_output=True, record_steps=False):
        S = dict(DEFAULTS)
        S.update(_normalise_keys(S0))
        _check_supported(S)
        self.S = S
        self.verbose = verbose
        self.save_output = save_output
        self.record_steps = record_steps
        self.steps = []
        self.timings = {}

        # ---- time
        t0 = datenum(S["reftime"])
        self.TIME = dict(dt=S["dt"], tc=S["tc"], timenum0=t0, tnow=t0, tend=datenum(S["endofsimulation"]),
                         it=-1, itout=0, tprev=t0, nt=0)

        # ---- coastline (prepare_coastline.m)
        x_mc = np.asarray(S["xmc"], dtype=float).ravel()
        y_mc = np.asarray(S["ymc"], dtype=float).ravel()
        ok = np.flatnonzero(~np.isnan(x_mc))
        x_mc = x_mc[ok[0]:ok[-1] + 1]
        y_mc = y_mc[ok[0]:ok[-1] + 1]
        idnan = np.flatnonzero(np.isnan(x_mc))
        drop = idnan[:-1][np.diff(idnan) == 1] if idnan.size > 1 else np.zeros(0, int)
        keep = np.setdiff1d(np.arange(x_mc.size), drop)
        x_mc, y_mc = x_mc[keep], y_mc[keep]
        self.COAST = dict(
            x_mc=x_mc, y_mc=y_mc, n_mc=int(np.sum(np.isnan(x_mc))) + 1,
            h0input=S["d"], ds0=float(S["ds0"]), twopoints=S["twopoints"], smoothfac=S["smoothfac"],
            smoothrefrac=S["smoothrefrac"], PHIf0=None, tanbeta=S["tanbeta"],
            griddingmethod=S["griddingmethod"], maxangle=S["maxangle"],
            preserveorientation=S["preserveorientation"], PHIc0bnd=np.array([np.nan, np.nan]),
            boundaryconditionstart=S["boundaryconditionstart"],
            boundaryconditionend=S["boundaryconditionend"], xq_mc=np.zeros(0), yq_mc=np.zeros(0),
        )

        # ---- waves (prepare_waveconditions.m)
        t_read = time.perf_counter()
        WVC = read_wave_stations(S["wvcfile"])
        DT0 = np.median(np.diff(WVC["timenum"]))
        factor = max(int(np.ceil(S["dt"] * 365 / DT0)), 1)
        if factor > 2:
            raise NotImplementedError("aggregation of wave time series (time step > 2x wave interval)")
        self.timings["read_waves"] = time.perf_counter() - t_read
        self.WAVE = dict(WVC=WVC, ddeep=S["ddeep"], dnearshore=S["dnearshore"], gamma=S["gamma"],
                         interpolationmethod=S["interpolationmethod"], sphimax=S["sphimax"])

        # ---- transport (prepare_transport.m)
        self.TRANSP = dict(trform=S["trform"], qscal=S["qscal"], d50=S["d50"], porosity=S["porosity"],
                           tanbeta=S["tanbeta"], rhos=S["rhos"], rhow=S["rhow"], g=S["g"],
                           gamma=S["gamma"], twopoints=S["twopoints"],
                           suppresshighangle=S["suppresshighangle"],
                           relaxationlength=S["relaxationlength"],
                           boundaryconditionstart=S["boundaryconditionstart"],
                           boundaryconditionend=S["boundaryconditionend"])

        # ---- spit / overwash (prepare_spit.m)
        owscale = S["owscale"] if S["owtimescale"] <= 0 else 0.0
        self.SPIT = dict(method=S["spitmethod"], spitwidth=S["spitwidth"], owscale=owscale,
                         owtimescale=S["owtimescale"], Dsf=S["spitdsf"], Dbb=S["spitdbb"],
                         bheight=S["bheight"])

        self.O = {f: np.zeros((0, 0)) for f in _O_FIELDS}

    # ------------------------------------------------------------------
    def _segment(self, i_mc):
        """get_segmentdata.m: data of section i_mc after the transport phase."""
        C, W, T = self.COAST, self.WAVE, self.TRANSP
        x, y, _, i1, i2 = get_one_polygon(C["x_mc"], C["y_mc"], i_mc)
        seg = dict(x=x, y=y, i1=i1, i2=i2, i_mc=i_mc)
        seg["s"] = get_one_section(C["s_mc"], i_mc)
        seg["ds"] = get_one_section(C["ds_mc"], i_mc)
        seg["h0"] = get_one_section(C["h0_mc"], i_mc)
        seg["PHIcxy0"] = get_one_section(C["PHIcxy0_mc"], i_mc)
        QS, PHI, _, _, _ = get_one_polygon(T["QS_mc"], W["PHIbr_mc"], i_mc)
        seg["QS"] = QS
        seg["PHI"] = PHI
        seg["shadow"] = get_one_section(T["shadow_mc"], i_mc)
        seg["cyclic"] = get_cyclic(x, y, C["ds0"])
        seg["n"] = x.size
        seg["nq"] = QS.size
        return seg

    def _store(self):
        """save_shorelines.m (without output grids / profiles)."""
        S, TIME, C, W, T, O = self.S, self.TIME, self.COAST, self.WAVE, self.TRANSP, self.O
        storedata = ((TIME["tnow"] - TIME["timenum0"]) >= TIME["itout"] * S["storageinterval"]
                     or (TIME["tnow"] + TIME["dt"]) >= TIME["tend"])
        if not storedata:
            return
        x_mc, y_mc = C["x_mc"], C["y_mc"]
        x_mc1, y_mc1 = C["x1_mc"], C["y1_mc"]

        def dist(x, y):
            dx = x[1:] - x[:-1]
            dy = y[1:] - y[:-1]
            dx[np.isnan(dx)] = 0
            dy[np.isnan(dy)] = 0
            d = np.concatenate([[0.0], np.cumsum(m_pow(dx * dx + dy * dy, 0.5))])
            return d, (d[:-1] + d[1:]) / 2

        distx, distQS = dist(x_mc, y_mc)
        distx1, distQS1 = dist(x_mc1, y_mc1)
        O["it"] = _add_vector(O["it"], TIME["it"])
        O["dt"] = _add_vector(O["dt"], TIME["dt"])
        O["tc"] = _add_vector(O["dt"], TIME["tc"])        # (as in ShorelineS)
        O["nt"] = _add_vector(O["nt"], TIME["nt"])
        O["timenum"] = _add_vector(O["timenum"], TIME["tnow"])
        O["n"] = _add_vector(O["n"], np.sum(np.isnan(x_mc)) + 1)
        O["x"] = _add_vector(O["x"], x_mc)
        O["y"] = _add_vector(O["y"], y_mc)
        O["distx"] = _add_vector(O["distx"], distx)
        O["distQS"] = _add_vector(O["distQS"], distQS)
        O["dSds"] = _add_vector(O["dSds"], C["dSds_mc"])
        O["n1"] = _add_vector(O["n1"], np.sum(np.isnan(x_mc1)) + 1)
        O["x1"] = _add_vector(O["x1"], x_mc1)
        O["y1"] = _add_vector(O["y1"], y_mc1)
        O["distx1"] = _add_vector(O["distx1"], distx1)
        O["distQS1"] = _add_vector(O["distQS1"], distQS1)
        O["PHIc"] = _add_vector(O["PHIc"], C["PHIc_mc"])
        O["QS"] = _add_vector(O["QS"], T["QS_mc"])
        O["HSo"] = _add_vector(O["HSo"], W["HSo_mc"])
        O["PHIo"] = _add_vector(O["PHIo"], W["PHIo_mc"])
        O["TP"] = _add_vector(O["TP"], W["TP_mc"])
        O["PHIf"] = _add_vector(O["PHIf"], C["PHIf_mc"])
        O["HS"] = _add_vector(O["HS"], W["HStdp_mc"])
        O["PHI"] = _add_vector(O["PHI"], W["PHItdp_mc"])
        O["dPHI"] = _add_vector(O["dPHI"], W["dPHItdp_mc"])
        O["HSbr"] = _add_vector(O["HSbr"], W["HSbr_mc"])
        O["PHIbr"] = _add_vector(O["PHIbr"], W["PHIbr_mc"])
        O["dPHIbr"] = _add_vector(O["dPHIbr"], W["dPHIbr_mc"])
        O["hbr"] = _add_vector(O["hbr"], W["hbr_mc"])
        for f in ("xhard", "yhard", "xnour", "ynour", "x_groyne", "y_groyne"):
            O[f] = _add_vector(O[f], np.zeros(0))
        O["nhard"] = _add_vector(O["nhard"], 0)
        O["nnour"] = _add_vector(O["nnour"], 0)
        TIME["itout"] += 1
        if self.save_output:
            self.save()

    def save(self, fname="output.mat"):
        outdir = self.S["outputdir"]
        os.makedirs(outdir, exist_ok=True)
        O = dict(self.O)
        O["outputdir"] = outdir
        O["storageinterval"] = self.S["storageinterval"]
        scipy.io.savemat(os.path.join(outdir, fname), {"O": O}, do_compression=False)

    # ------------------------------------------------------------------
    def step(self):
        C, W, T, TIME, S = self.COAST, self.WAVE, self.TRANSP, self.TIME, self.S
        TIME["it"] += 1
        TIME["nt"] = TIME["it"]

        # PHASE 1: transport, per coastline section
        for i_mc in range(1, C["n_mc"] + 1):
            cl.make_sgrid_mc(C, TIME, i_mc)
            introduce_wave(W, TIME, C)
            cl.get_coastline_orientation(C)
            cl.get_foreshore_orientation(C)
            cl.get_active_profile(C)
            wave_refraction(W, C)
            W["diff"] = np.zeros(C["nq"])
            W["dPHIo"] = atan2d(sind(C["PHIc"] - W["PHIo"]), cosd(C["PHIc"] - W["PHIo"]))
            W["dPHItdp"] = atan2d(sind(C["PHIcs"] - W["PHItdp"]), cosd(C["PHIcs"] - W["PHItdp"]))
            wave_breakingheight(W, T)
            W["PHIbr"] = mmod(C["PHIc"] - W["dPHIbr"], 360)
            wave_angles(C, W, T)
            T["QS"] = transport(T, W)
            transport_shadow_treat(C, W, T, find_shadows_mc)
            T["idrev"] = np.zeros(T["QS"].size)
            get_smoothangles(C, T)
            get_upwindcorrection(C, W, T)
            transport_boundary_condition(T, C)
            # collect_variables
            C["x1"], C["y1"] = C["x"], C["y"]
            C["xq1"], C["yq1"] = C["xq"], C["yq"]
            _collect(C, C, _COAST_COLLECT, i_mc)
            _collect(W, W, _WAVE_COLLECT, i_mc)
            _collect(T, T, _TRANSP_COLLECT, i_mc)
            W["diff_mc"] = np.nan_to_num(W["diff_mc"], nan=0.0).astype(bool)
            if i_mc == C["n_mc"]:
                C["x_mc"], C["y_mc"] = C["x1_mc"].copy(), C["y1_mc"].copy()
                C["xq_mc"], C["yq_mc"] = C["xq1_mc"].copy(), C["yq1_mc"].copy()

        # PHASE 2: coastline change
        for i_mc in range(1, C["n_mc"] + 1):
            seg = self._segment(i_mc)
            C.update(i_mc=i_mc, x=seg["x"], y=seg["y"], n=seg["n"], nq=seg["nq"], h0=seg["h0"],
                     cyclic=seg["cyclic"], s=seg["s"], ds=seg["ds"], PHIcxy0=seg["PHIcxy0"])
            T["QS"] = seg["QS"]
            cl.coastline_change(C, T, TIME)

        # PHASE 3: overwash, merging, cleanup
        def seg_overwash(ii):
            s = self._segment(ii)
            return s["x"], s["y"], s["shadow"], s["PHI"], s["i1"], s["cyclic"]

        cl.find_overwash_mc(C, self.SPIT, TIME, seg_overwash)
        cl.get_transportpoints(C)                    # get_reconnectedgroynes (no groynes)
        for i_mc in range(1, C["n_mc"] + 1):
            cl.merge_coastlines(C, i_mc)
        cl.merge_coastlines_mc(C)
        cl.cleanup_nans(C)

        # PHASE 4: storage
        self._store()
        if self.record_steps:
            self.steps.append(dict(x_mc=C["x_mc"].copy(), y_mc=C["y_mc"].copy(), QS=T["QS_mc"].copy(),
                                   tnow=TIME["tnow"]))
        TIME["tprev"] = TIME["tnow"]
        TIME["tnow"] = TIME["tnow"] + TIME["dt"] * 365

    def run(self):
        TIME = self.TIME
        if self.verbose:
            print("  Loop over time")
        t0 = time.perf_counter()
        while TIME["tnow"] < TIME["tend"] or (TIME["tnow"] == TIME["tend"] and TIME["tc"] == 0):
            self.step()
            if self.verbose:
                W = self.WAVE
                h = W["HStdp_mc"][~np.isnan(W["HStdp_mc"])]
                p = W["PHItdp_mc"][~np.isnan(W["PHItdp_mc"])]
                phiavg = mmod(atan2d(np.median(sind(p)), np.median(cosd(p))), 360)
                print(f"   {datestr(TIME['tnow'])}  Hs_tdp={np.median(h):2.1f}, Dir_tdp={float(phiavg):1.0f}")
        self.timings["time_loop"] = time.perf_counter() - t0
        if self.save_output:
            self.save()
        return self.O


def run(S, **kwargs):
    """Run ShorelineS with settings S (dict, MATLAB field names). Returns (O, model)."""
    model = ShorelineS(S, **kwargs)
    O = model.run()
    return O, model
