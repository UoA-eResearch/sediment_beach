"""Wave input, interpolation along the coast, refraction, breaking and the
critical (maximum-transport) wave angle. Ports of introduce_wave.m,
get_interpolation_on_grid.m, wave_refraction.m, wave_breakingheight.m,
wave_angles.m and get_Sphimax.m for the configuration used in this project."""
import numpy as np
import scipy.linalg

from .mfuncs import (acosd, asind, atan2d, cosd, datenum, get_disper, interp1,
                     interpNANs, m_hypot, m_pow, mmod, sind)
from .transport import transport

_EPS = np.finfo(float).eps


# ---------------------------------------------------------------------------
# Wave input
# ---------------------------------------------------------------------------
def read_wave_stations(ncfile):
    """Read a ShorelineS wave-station NetCDF file (station_x/y, point_hm0,
    point_tp, point_wavdir, time in '<unit> since <date>')."""
    import netCDF4

    with netCDF4.Dataset(ncfile) as d:
        x = np.asarray(d["station_x"][:].filled(np.nan), dtype=float)
        y = np.asarray(d["station_y"][:].filled(np.nan), dtype=float)

        def var(name):
            v = d[name]
            a = np.asarray(v[:].filled(np.nan) if np.ma.isMaskedArray(v[:]) else v[:], dtype=float)
            # Python reads (station, time) when the MATLAB dimension order is {time, station}
            if a.shape[0] != x.size:
                a = a.T
            return np.ascontiguousarray(a)

        hm0 = var("point_hm0")
        tp = var("point_tp")
        wd = var("point_wavdir")
        units = d["time"].units
        tt = np.asarray(d["time"][:], dtype=float)
    reftime = datenum(units[14:])
    timenum = reftime + tt / 24 / 60 / 60
    return {"x": x, "y": y, "timenum": timenum, "Hs": hm0, "Tp": tp, "Dir": wd}


# ---------------------------------------------------------------------------
# get_interpolation_on_grid
# ---------------------------------------------------------------------------
def _dist(xq, yq, xw, yw):
    dxq = xq - xw
    dyq = yq - yw
    return m_pow(dxq * dxq + dyq * dyq, 0.5)


def interp_weighted_distance(xq, yq, xw, yw, var1, var2):
    """'weighted_distance' method: inverse-distance weights of the 2 nearest points."""
    xq = np.asarray(xq, dtype=float)
    yq = np.asarray(yq, dtype=float)
    xw = np.asarray(xw, dtype=float).ravel()
    yw = np.asarray(yw, dtype=float).ravel()
    nq = xq.size
    var1i = {k: np.zeros((v.shape[0], nq)) for k, v in var1.items()}
    var2i = {k: np.zeros((v.shape[0], nq)) for k, v in var2.items()}
    distw = np.zeros(nq)
    for ii in range(nq):
        dist = _dist(xq[ii], yq[ii], xw, yw)
        ids = np.argsort(dist, kind="stable")[:2]
        dists = dist[ids]
        distw[ii] = np.min(dists)
        if ids.size >= 2:
            wght0 = 1 - (dists / max(dists[0] + dists[1], _EPS))
            wght = wght0 / max(wght0[0] + wght0[1], _EPS)
            for k, v in var1.items():
                p = v[:, ids] * wght
                var1i[k][:, ii] = p[:, 0] + p[:, 1]
            for k, v in var2.items():
                c = cosd(v[:, ids]) * wght
                s = sind(v[:, ids]) * wght
                var2i[k][:, ii] = mmod(atan2d(s[:, 0] + s[:, 1], c[:, 0] + c[:, 1]), 360)
        elif ids.size == 1:
            for k, v in var1.items():
                var1i[k][:, ii] = v[:, ids[0]]
            for k, v in var2.items():
                var2i[k][:, ii] = mmod(v[:, ids[0]], 360)
    return var1i, var2i, distw


def interp_alongshore_mapping(xq, yq, xw, yw, var1, var2):
    """'alongshore_mapping' method: map each data point to its nearest grid
    point, drop points that are far offshore of their neighbours, and
    interpolate linearly along the coast."""
    xq = np.asarray(xq, dtype=float).ravel()
    yq = np.asarray(yq, dtype=float).ravel()
    xw = np.asarray(xw, dtype=float).ravel()
    yw = np.asarray(yw, dtype=float).ravel()
    nq = xq.size
    var1 = {k: np.asarray(v, dtype=float) for k, v in var1.items()}
    var2 = {k: np.asarray(v, dtype=float) for k, v in var2.items()}
    dxq = np.diff(xq)
    dyq = np.diff(yq)
    distq = np.concatenate([[0.0], np.cumsum(m_pow(dxq * dxq + dyq * dyq, 0.5))])
    nw = xw.size
    distw = np.full(nq, np.nan)
    # distance of every data point to every grid point (nw x nq)
    D = _dist(xq[None, :], yq[None, :], xw[:, None], yw[:, None])
    dmin = np.nanmin(D, axis=1)
    idGRID = np.argmax(D == dmin[:, None], axis=1)        # first grid point at minimum distance
    dist0 = distq[idGRID]
    # distw(id) = distance of the (last) data point mapped to grid point id
    u, last = np.unique(idGRID[::-1], return_index=True)
    distw[u] = dmin[nw - 1 - last]
    distw = interpNANs(distw)
    _, idu = np.unique(dist0, return_index=True)
    IDS = np.sort(idu)                   # intersect(sort-permutation, idu) = sorted idu
    dist0 = dist0[IDS]
    idGRID = idGRID[IDS]
    for k in var1:
        var1[k] = var1[k][:, IDS]
    for k in var2:
        var2[k] = var2[k][:, IDS]
    xw = xw[IDS]
    yw = yw[IDS]
    Dw = D[IDS]
    # data points mapped to the same grid point: keep only the closest one(s)
    dmin_w = np.nanmin(Dw, axis=1)
    dist3 = dmin_w.copy()
    rem = []
    vals, counts = np.unique(idGRID, return_counts=True)
    for v in vals[counts > 1]:
        ID = np.flatnonzero(idGRID == v)
        dist2 = dmin_w[ID]
        rem.extend(ID[dist2 != np.min(dist2)])
        dist3[ID] = np.min(dist2)
    idDATA = np.setdiff1d(np.arange(idGRID.size), np.asarray(rem, dtype=int))
    dist0 = dist0[idDATA]
    for k in var1:
        var1[k] = var1[k][:, idDATA]
    for k in var2:
        var2[k] = var2[k][:, idDATA]
    idGRID = idGRID[idDATA]

    # drop data points that lie much further offshore than their neighbours
    dx0 = m_pow((xq[-1] - xq[0]) ** 2 + (yq[-1] - yq[0]) ** 2, 0.5)
    for _ in range(10):
        dxend = distq[-1] - dist0[-1]
        dx1 = dist0[0] + dxend
        dxv = np.concatenate([[dx0 + dx1], np.diff(dist0), [dx0 + dx1]])
        dcross = np.concatenate([[dist3[0] - dist3[-1]], np.diff(dist3), [dist3[0] - dist3[-1]]])
        dist3p = np.concatenate([[dist3[-1]], dist3, [dist3[0]]])
        m = idGRID.size
        xx = np.arange(m)
        drop = (((dcross[xx] > dxv[xx]) & (dist3p[xx + 1] > 2 * dist3p[xx]))
                | ((-dcross[xx + 1] > dxv[xx + 1]) & (dist3p[xx + 1] > 2 * dist3p[xx + 2])))
        keep = np.flatnonzero(~drop)
        dist0 = dist0[keep]
        dist3 = dist3p[keep + 1]
        idGRID = idGRID[keep]
        for k in var1:
            var1[k] = var1[k][:, keep]
        for k in var2:
            var2[k] = var2[k][:, keep]

    var1i, var2i = {}, {}
    if dist0.size >= 2:
        for k, v in var1.items():
            var1i[k] = np.vstack([interpNANs(interp1(dist0, v[g], distq)) for g in range(v.shape[0])])
        for k, v in var2.items():
            rows = []
            for g in range(v.shape[0]):
                cphi = interpNANs(interp1(dist0, cosd(v[g]), distq))
                sphi = interpNANs(interp1(dist0, sind(v[g]), distq))
                rows.append(mmod(atan2d(sphi, cphi), 360))
            var2i[k] = np.vstack(rows)
    elif dist0.size == 1:
        for k, v in var1.items():
            var1i[k] = np.repeat(v, nq, axis=1)
        for k, v in var2.items():
            var2i[k] = np.repeat(mmod(v, 360), nq, axis=1)
    else:
        for k, v in var1.items():
            var1i[k] = np.zeros((v.shape[0], nq))
        for k, v in var2.items():
            var2i[k] = np.zeros((v.shape[0], nq))
    return var1i, var2i, distw


def interpolate_on_grid(method, xq, yq, xw, yw, var1, var2):
    if ("weight" in method) or ("distance" in method):
        return interp_weighted_distance(xq, yq, xw, yw, var1, var2)
    if ("alongshore" in method) or ("mapping" in method):
        return interp_alongshore_mapping(xq, yq, xw, yw, var1, var2)
    raise NotImplementedError(f"interpolation method '{method}'")


# ---------------------------------------------------------------------------
# introduce_wave
# ---------------------------------------------------------------------------
def introduce_wave(WAVE, TIME, COAST):
    """Offshore wave conditions (HSo, PHIo, TP) on the transport points of the
    current coastline section, interpolated in time and along the coast."""
    WVC = WAVE["WVC"]
    t = WVC["timenum"]
    tnow = TIME["tnow"]
    if np.any(t != t[0]):
        if tnow >= np.max(t) or tnow <= np.min(t):
            # recycle the wave time series (as introduce_wave.m)
            dt = t - tnow
            mdt = mmod(dt, 365.25)
            idt = int(np.flatnonzero(mdt < np.min(mdt) + 1 / 24)[0])
            dtrewind = dt[idt]
            if np.min(mdt) > 30 or (np.max(t) - np.min(t)) < 120:
                dtrewind = dt[0]
            WVC["timenum"] = t = t - dtrewind
            WAVE["dtrewind"] = WAVE.get("dtrewind", 0.0) + dtrewind
        it = np.flatnonzero(t <= tnow)
        if it.size and it[-1] < t.size - 1 and np.all(np.diff(t) > 0):
            i = it[-1]
            iw = [i, i + 1]
        else:
            iw = slice(None)
        tw = t[iw]
        Hs = interp1(tw, WVC["Hs"][:, iw].T, tnow)
        Tp = interp1(tw, WVC["Tp"][:, iw].T, tnow)
        D = WVC["Dir"][:, iw].T
        Dir = mmod(atan2d(interp1(tw, sind(D), tnow), interp1(tw, cosd(D), tnow)), 360)
    else:
        raise NotImplementedError("wave climates (time-independent wave input)")
    Hs = np.atleast_1d(Hs)[None, :]
    Tp = np.atleast_1d(Tp)[None, :]
    Dir = np.atleast_1d(Dir)[None, :]
    var1N, var2N, distw = interpolate_on_grid(WAVE["interpolationmethod"], COAST["xq"], COAST["yq"],
                                              WVC["x"], WVC["y"], {"Hs": Hs, "Tp": Tp}, {"Dir": Dir})
    eps = 0.01
    WAVE["HSo"] = np.fmax(var1N["Hs"][0], eps) * 1.0          # CC.HScor = 1
    WAVE["PHIo"] = mmod(var2N["Dir"][0] + 0.0, 360.0)          # CC.PHIcor = 0
    WAVE["TP"] = np.fmax(var1N["Tp"][0], eps)
    WAVE["dist"] = distw
    WAVE["i_mc"] = COAST["i_mc"]
    return WAVE


# ---------------------------------------------------------------------------
# Refraction, breaking
# ---------------------------------------------------------------------------
def _real_asind(z):
    """real(asind(z)) of MATLAB/Octave, also for |z| > 1."""
    z = np.asarray(z, dtype=float)
    with np.errstate(invalid="ignore"):
        out = asind(np.clip(z, -1.0, 1.0))
    big = (np.pi / 2) * 180.0 / np.pi
    out = np.where(z > 1, big, out)
    out = np.where(z < -1, -big, out)
    return out


def wave_refraction(WAVE, COAST):
    """Refract the offshore waves to the nearshore depth (Snell + shoaling)."""
    PHIf = COAST["PHIf"]
    relANGLE = mmod(WAVE["PHIo"] - PHIf + 180, 360) - 180
    off = np.abs(relANGLE) >= 89.999
    relANGLE = np.where(off, np.sign(relANGLE) * 89.999, relANGLE)
    _, c_deep, _, n_deep = get_disper(WAVE["ddeep"], WAVE["TP"])
    _, c_tdp, _, n_tdp = get_disper(WAVE["dnearshore"], WAVE["TP"])
    WAVE["PHItdp"] = mmod(PHIf + _real_asind(c_tdp / c_deep * sind(relANGLE)), 360)
    relANGLEtdp = mmod(WAVE["PHItdp"] - PHIf + 180, 360) - 180
    with np.errstate(invalid="ignore", divide="ignore"):
        arg = n_deep * c_deep * cosd(relANGLE) / (n_tdp * c_tdp * cosd(relANGLEtdp))
    WAVE["HStdp"] = WAVE["HSo"] * np.fmax(np.sqrt(np.fmax(arg, 0.0)), 0.0)
    return WAVE


def _shoalref_vec(hbr, tper, gamma, hstdp, ctdp, ntdp, sinPHIw, cosPHIw):
    _, cbr, _, nbr = get_disper(hbr, tper)
    dPHIbr = acosd(cosPHIw)
    hstdpbr = hstdp.copy()
    with np.errstate(invalid="ignore", divide="ignore"):
        in1 = np.abs(cbr / ctdp * sinPHIw) < 1
        dPHIbr[in1] = asind(cbr[in1] / ctdp[in1] * sinPHIw[in1])
        cosbr = cosd(dPHIbr)
        in2 = in1 & ((nbr * cbr * cosbr) > 0) & (np.abs(cosbr) > 1e-3)
        hstdpbr[in2] = hstdp[in2] * np.sqrt(ntdp[in2] * ctdp[in2] * cosPHIw[in2] / (nbr[in2] * cbr[in2] * cosbr[in2]))
    sel = in1 & ~in2
    dPHIbr[sel] = acosd(cosPHIw[sel])
    hbrnew = hstdpbr / gamma
    return hbrnew, dPHIbr, hbrnew - hbr, cbr, nbr


def wave_breakingheight(WAVE, TRANSP):
    """Breaking wave height and angle by secant iteration of shoaling and
    refraction from the nearshore point to the breaker line (vectorised)."""
    eps = 1e-5
    if TRANSP["suppresshighangle"] == 1:
        raise NotImplementedError("suppresshighangle")
    if TRANSP["trform"].upper() in ("RAY", "CERC", "CERC2"):
        raise NotImplementedError(f"transport formula {TRANSP['trform']}")
    TP = np.asarray(WAVE["TP"], dtype=float)
    HStdp = np.asarray(WAVE["HStdp"], dtype=float)
    gamma = WAVE["gamma"]
    dPHItdp = np.asarray(WAVE["dPHItdp"], dtype=float)
    _, ctdp, _, ntdp = get_disper(WAVE["dnearshore"], TP)
    cosPHI = np.fmax(cosd(dPHItdp), eps)
    sinPHI = sind(dPHItdp)
    hbr1 = np.fmax(HStdp / gamma, eps)
    hbr2, _, err1, _, _ = _shoalref_vec(hbr1, TP, gamma, HStdp, ctdp, ntdp, sinPHI, cosPHI)
    hbr2b = np.fmax(hbr2, eps)
    hbr3, dPHIbr0, err2, cbr0, nbr0 = _shoalref_vec(hbr2b, TP, gamma, HStdp, ctdp, ntdp, sinPHI, cosPHI)
    hbrnew = hbr3.copy()
    act = np.flatnonzero(~(np.abs(err2) < eps))
    for _ in range(3, 11):
        if act.size == 0:
            break
        with np.errstate(invalid="ignore", divide="ignore"):
            hn = np.fmax(hbr1[act] - err1[act] * (hbr2[act] - hbr1[act]) / (err2[act] - err1[act]), eps)
        _, dPHIa, erra, cbra, nbra = _shoalref_vec(hn, TP[act], gamma, HStdp[act], ctdp[act], ntdp[act],
                                                   sinPHI[act], cosPHI[act])
        hbrnew[act] = hn
        dPHIbr0[act] = dPHIa
        cbr0[act] = cbra
        nbr0[act] = nbra
        with np.errstate(invalid="ignore"):
            cont = np.abs(erra) > eps
        ac = act[cont]
        hbr1[ac] = hbr2[ac]
        err1[ac] = err2[ac]
        hbr2[ac] = hn[cont]
        err2[ac] = erra[cont]
        act = ac
    hbrnew[np.isnan(hbrnew)] = 0
    dPHIbr0[np.isnan(dPHIbr0)] = 0
    WAVE["HSbr"] = hbrnew * gamma
    WAVE["dPHIbr"] = dPHIbr0
    WAVE["hbr"] = hbrnew
    WAVE["cbr"] = cbr0
    WAVE["nbr"] = nbr0
    return WAVE


# ---------------------------------------------------------------------------
# Critical wave angle (maximum transport)
# ---------------------------------------------------------------------------
def _solve3(A, b):
    """A[k] \\ b for a batch of 3x3 systems (LAPACK LU, as MATLAB's mldivide)."""
    try:
        return np.linalg.solve(A, b[..., None])[..., 0]
    except np.linalg.LinAlgError:
        out = np.empty_like(b)
        for k in range(A.shape[0]):
            with np.errstate(all="ignore"):
                lu, piv = scipy.linalg.lu_factor(A[k], check_finite=False)
                out[k] = scipy.linalg.lu_solve((lu, piv), b[k], check_finite=False)
        return out


def _mmin(a):
    """MATLAB min of a vector (NaNs ignored unless all are NaN)."""
    return np.nan if np.all(np.isnan(a)) else np.nanmin(a)


def get_Sphimax(WAVE, TRANSP):
    eps = 1e-4
    nq = WAVE["HStdp"].size
    W1 = dict(WAVE)
    W1["dPHItdp"] = np.full(nq, 35.0)
    W1["HStdp"] = np.fmax(WAVE["HStdp"], eps)
    W1 = wave_breakingheight(W1, TRANSP)
    QS1 = transport(TRANSP, W1)
    dPHI1 = W1["dPHItdp"].copy()
    W2 = dict(WAVE)
    W2["dPHItdp"] = np.full(nq, 45.0)
    W2["HStdp"] = np.fmax(WAVE["HStdp"], eps)
    W2 = wave_breakingheight(W2, TRANSP)
    QS2 = transport(TRANSP, W2)
    dPHI2 = W2["dPHItdp"].copy()
    Wm = dict(W1)
    dPHIm = 0.5 * (dPHI1 + dPHI2)
    err = np.full(nq, 1e3)
    it = 0
    dPHIcrit = np.full(nq, np.nan)
    dPHIcritbr = np.full(nq, np.nan)
    QSmax = np.zeros(nq)
    while (_mmin(err) > eps or it <= 1) and it < 10:
        it += 1
        Wm["dPHItdp"] = dPHIm
        Wm = wave_breakingheight(Wm, TRANSP)
        QSm = transport(TRANSP, Wm)
        dPHImold = dPHIm.copy()
        inn = np.flatnonzero(err > eps)
        if inn.size:
            # A = [d1^2 d1 1; d2^2 d2 1; dm^2 dm 1]; ShorelineS uses a = A\B with
            # B=[QS1;QS2;QSm] (3 x nq) and then a(1), a(2), i.e. the first column
            # of the solution (the transport of the first grid point).
            d1, d2, dm = dPHI1[inn], dPHI2[inn], dPHIm[inn]
            A = np.stack([np.stack([d1 * d1, d1, np.ones_like(d1)], -1),
                          np.stack([d2 * d2, d2, np.ones_like(d2)], -1),
                          np.stack([dm * dm, dm, np.ones_like(dm)], -1)], 1)
            b = np.broadcast_to(np.array([QS1[0], QS2[0], QSm[0]]), (inn.size, 3)).copy()
            a = _solve3(A, b)
            with np.errstate(invalid="ignore", divide="ignore"):
                dPHIm[inn] = -a[:, 1] / (2 * a[:, 0])
        m = QS1 > QS2
        dPHI2[m] = dPHImold[m]
        m = QS1 <= QS2
        dPHI1[m] = dPHImold[m]
        m = QS1 > QS2
        QS2[m] = QSm[m]
        m = QS1 <= QS2
        QS1[m] = QSm[m]
        dPHIm[dPHIm > 179.99] = 179.99
        err = np.abs(dPHIm - dPHImold)
        dPHIm[dPHIm == 999] = np.nan
        dPHIcrit = dPHIm.copy()
        dPHIcritbr = Wm["dPHIbr"].copy()
        QSmax = QSm
    WAVE["dPHIcrit"] = interpNANs(dPHIcrit)
    WAVE["dPHIcritbr"] = interpNANs(dPHIcritbr)
    TRANSP["QSmax"] = interpNANs(QSmax)
    return WAVE, TRANSP


def wave_angles(COAST, WAVE, TRANSP):
    nq = COAST["nq"]
    _, WAVE["ctdp"], _, WAVE["ntdp"] = get_disper(np.full(nq, WAVE["dnearshore"]), WAVE["TP"])
    if WAVE.get("sphimax") not in (None, []) and np.size(WAVE["sphimax"]):
        raise NotImplementedError("S.sphimax")
    WAVE, TRANSP = get_Sphimax(WAVE, TRANSP)
    refraclimit = 180.0
    d = WAVE["dPHItdp"]
    d[(d - WAVE["dPHIo"]) < -refraclimit] = -refraclimit
    d[(d - WAVE["dPHIo"]) > refraclimit] = refraclimit
    d = WAVE["dPHIbr"]
    d[(d - WAVE["dPHIo"]) < -refraclimit] = -refraclimit
    d[(d - WAVE["dPHIo"]) > refraclimit] = refraclimit
    return WAVE, TRANSP


__all__ = ["read_wave_stations", "interpolate_on_grid", "introduce_wave", "wave_refraction",
           "wave_breakingheight", "get_Sphimax", "wave_angles", "acosd", "m_hypot"]
