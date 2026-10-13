"""Coastline grid, orientation and update. Ports of make_sgrid_mc.m
(griddingmethod 2), get_coastline_orientation.m, get_foreshore_orientation.m,
coastline_change.m, get_transportpoints.m, find_overwash_mc.m,
merge_coastlines.m, merge_coastlines_mc.m and cleanup_nans.m."""
import math

import numpy as np

from .geometry import get_intersections, lines_with_crossing
from .mfuncs import (atan2d, cosd, get_clockpoly, get_cyclic, get_mod,
                     get_nansremoved, get_one_polygon, get_one_section,
                     insert_section, interp1, m_hypot, mmod, sind)
from .waves import interpolate_on_grid

_DBL_EPS = np.finfo(float).eps


def msum(a):
    """MATLAB/Octave sum (sequential accumulation)."""
    a = np.asarray(a, dtype=float).ravel()
    return float(np.cumsum(a)[-1]) if a.size else 0.0


def _tfloor(x, ct):
    """Octave's tolerant floor used to count range elements."""
    q = 1.0
    if x < 0.0:
        q = 1.0 - ct
    rmax = q / (2.0 - ct)
    t1 = 1.0 + math.floor(x)
    t1 = (ct / q) * (-t1 if t1 < 0.0 else t1)
    t1 = rmax if rmax < t1 else t1
    t1 = ct if ct > t1 else t1
    t1 = math.floor(x + t1)
    if x <= 0.0 or (t1 - x) < rmax:
        return t1
    return t1 - 1.0


def colon(base, inc, limit):
    """Octave range ``base:inc:limit`` (positive increment)."""
    ct = 3.0 * _DBL_EPS
    n = int(_tfloor((limit - base + inc) / inc, ct))
    n = max(n, 0)
    if n == 0:
        return np.zeros(0)
    r = base + np.arange(n) * inc
    final = base + (n - 1) * inc
    if final > limit:
        final = limit
    r[-1] = final
    r[0] = base
    return r


# ---------------------------------------------------------------------------
def make_sgrid_mc(COAST, TIME, i_mc):
    """Regrid coastline section i_mc to the target grid size ds0."""
    eps = 0.1
    sqrt2 = math.sqrt(2.0)
    x, y, _, _, _ = get_one_polygon(COAST["x_mc"], COAST["y_mc"], i_mc)
    clockwise = get_clockpoly(x, y)
    COAST["gridchange"] = 0
    if COAST["griddingmethod"] != 2:
        raise NotImplementedError("griddingmethod 1")
    ds = m_hypot(np.diff(x), np.diff(y))
    s = np.concatenate([[0.0], np.cumsum(ds)])
    ds0 = COAST["ds0"]
    if np.max(ds) > 2 * ds0 or TIME["it"] == 0:
        COAST["gridchange"] = 1
        IDunique = np.concatenate([[True], ds > eps])
        x0, y0, s0 = x[IDunique], y[IDunique], s[IDunique]
        ns = math.ceil(s0[-1] / ds0)
        ds1 = s0[-1] / ns
        s = colon(0.0, ds1, s0[-1])
        x = interp1(s0, x0, s)
        y = interp1(s0, y0, s)
    snew = list(s)
    i = 2
    while i <= len(snew):
        ds2 = snew[i - 1] - snew[i - 2]
        if ds2 < ds0 / sqrt2 and snew[-1] >= ds0 / sqrt2:
            COAST["gridchange"] = 1
            if i > 2 and i < len(snew):
                del snew[i - 1]
            elif i == 2 and ds2 < ds0 / 4 and len(snew) > 2:
                del snew[1]
            elif ds2 < ds0 / 4:
                snew = snew[:i - 2] + [snew[-1]]
            else:
                i += 1
        elif ds2 > ds0 * sqrt2:
            COAST["gridchange"] = 1
            snew = snew[:i - 1] + [0.5 * (snew[i - 2] + snew[i - 1])] + snew[i - 1:]
            i += 1
        else:
            i += 1
    snew = np.asarray(snew, dtype=float)
    sf = COAST["smoothfac"]
    if snew.size > 2:
        snew[1:-1] = sf * snew[:-2] + (1.0 - 2 * sf) * snew[1:-1] + sf * snew[2:]
    if snew.size != s.size or msum(s) != msum(snew):
        if clockwise == 1 or np.min(s) < ds0 / 5:
            x = interp1(s, x, snew)
            y = interp1(s, y, snew)
            s = snew
    cyclic = get_cyclic(x, y, COAST["ds0"])
    if cyclic:
        x[-1] = x[0]
        y[-1] = y[0]
    d0 = np.diff(s)
    if d0.size:
        if cyclic:
            ds = np.concatenate([[(d0[0] + d0[-1]) / 2], (d0[:-1] + d0[1:]) / 2, [(d0[0] + d0[-1]) / 2]])
        else:
            ds = np.concatenate([d0[:1], (d0[:-1] + d0[1:]) / 2, d0[-1:]])
        clockwise = get_clockpoly(x, y)
        if not cyclic:
            clockwise = 1
        xq, yq = _transport_points(x, y, cyclic)
    else:
        ds = np.zeros(0)
        xq = np.zeros(0)
        yq = np.zeros(0)
        clockwise = 0
        cyclic = False
    COAST["i_mc"] = i_mc
    COAST["x"] = x
    COAST["y"] = y
    COAST["n"] = x.size
    COAST["s"] = s
    COAST["ds"] = ds
    COAST["xq"] = xq
    COAST["yq"] = yq
    COAST["dsq"] = np.concatenate([ds[:1], (ds[:-1] + ds[1:]) / 2, ds[-1:]])
    COAST["nq"] = xq.size
    COAST["cyclic"] = cyclic
    COAST["clockwise"] = clockwise
    COAST["x_mc"], COAST["y_mc"] = insert_section(x, y, COAST["x_mc"], COAST["y_mc"], i_mc)
    COAST["n_mc"] = int(np.sum(np.isnan(COAST["x_mc"]))) + 1
    return COAST


def _transport_points(x, y, cyclic):
    xq = (x[:-1] + x[1:]) / 2
    yq = (y[:-1] + y[1:]) / 2
    if cyclic:
        return np.concatenate([xq[-1:], xq, xq[:1]]), np.concatenate([yq[-1:], yq, yq[:1]])
    # no groynes: BNDgroyne == 0 at both ends
    xq = np.concatenate([[1.5 * x[0] - 0.5 * x[1]], xq, [1.5 * x[-1] - 0.5 * x[-2]]])
    yq = np.concatenate([[1.5 * y[0] - 0.5 * y[1]], yq, [1.5 * y[-1] - 0.5 * y[-2]]])
    return xq, yq


def get_transportpoints(COAST, i_mc_range=None):
    full = i_mc_range is None
    if full:
        i_mc_range = range(1, COAST["n_mc"] + 1)
        COAST["xq_mc"] = np.zeros(0)
        COAST["yq_mc"] = np.zeros(0)
    for i_mc in i_mc_range:
        x, y, _, _, _ = get_one_polygon(COAST["x_mc"], COAST["y_mc"], i_mc)
        if x.size > 1:
            xq, yq = _transport_points(x, y, get_cyclic(x, y, COAST["ds0"]))
        else:
            xq, yq = x, y
        if COAST["xq_mc"].size == 0:
            COAST["xq_mc"], COAST["yq_mc"] = xq, yq
        elif full and len(i_mc_range) > 1:
            COAST["xq_mc"] = np.concatenate([COAST["xq_mc"], [np.nan], xq])
            COAST["yq_mc"] = np.concatenate([COAST["yq_mc"], [np.nan], yq])
        else:
            COAST["xq_mc"], COAST["yq_mc"] = insert_section(xq, yq, COAST["xq_mc"], COAST["yq_mc"], i_mc)
    return COAST


# ---------------------------------------------------------------------------
def get_smoothdata_angle(VAR):
    """get_smoothdata(VAR,'angle',1): one pass of a [1 2 1]/4 vector smoothing."""
    s = sind(VAR)
    c = cosd(VAR)
    savg = np.concatenate([s[:1], (s[:-1] + s[1:]) / 2, s[-1:]])
    cavg = np.concatenate([c[:1], (c[:-1] + c[1:]) / 2, c[-1:]])
    s2 = (savg[:-1] + savg[1:]) / 2
    c2 = (cavg[:-1] + cavg[1:]) / 2
    return mmod(atan2d(s2, c2), 360)


def get_coastline_orientation(COAST):
    PHIc = mmod(360.0 - atan2d(np.diff(COAST["y"]), np.diff(COAST["x"])), 360)
    COAST["PHIcxy"] = mmod(360.0 - atan2d(np.diff(COAST["yq"]), np.diff(COAST["xq"])), 360)
    if np.isnan(COAST["PHIc0bnd"][0]):
        COAST["PHIc0bnd"][0] = PHIc[0]
    if np.isnan(COAST["PHIc0bnd"][1]):
        COAST["PHIc0bnd"][1] = PHIc[-1]
    if COAST["gridchange"] == 1:
        COAST["PHIcxy0"] = COAST["PHIcxy"]
    else:
        COAST["PHIcxy0"] = get_one_section(COAST["PHIcxy0_mc"], COAST["i_mc"])
    if COAST["cyclic"]:
        PHIc = np.concatenate([PHIc[-1:], PHIc, PHIc[:1]])
    else:
        PHIc = np.concatenate([PHIc[:1], PHIc, PHIc[-1:]])
        for key in ("boundaryconditionstart", "boundaryconditionend"):
            if COAST[key].lower() in ("angleconstant", "gradient"):
                raise NotImplementedError(f"{key} {COAST[key]}")
    if COAST["smoothrefrac"] > 0:
        raise NotImplementedError("S.smoothrefrac")
    COAST["PHIc"] = PHIc
    COAST["PHIcs"] = PHIc
    return COAST


def get_foreshore_orientation(COAST):
    """Foreshore orientation PHIf: smoothed coastline angle of the first time
    step, mapped onto the current transport points."""
    PHIf0 = COAST["PHIf0"]
    if PHIf0 is None or (isinstance(PHIf0, np.ndarray) and PHIf0.ndim == 2 and PHIf0.shape[1] == 4):
        nl = COAST["PHIc"].size
        xq = COAST["xq"]
        yq = COAST["yq"]
        if COAST["cyclic"]:
            raise NotImplementedError("cyclic coastline sections")
        PHIfsmooth = mmod(get_smoothdata_angle(COAST["PHIc"]), 360)
        idx = np.arange(nl)
        rows = np.stack([xq[idx], yq[idx], PHIfsmooth, PHIfsmooth], axis=1)
        if COAST["i_mc"] == 1:
            PHIf0 = rows
        else:
            PHIf0 = np.vstack([PHIf0, rows])
        if COAST["i_mc"] == COAST["n_mc"]:
            PHIf0 = PHIf0[:, :3]
        COAST["PHIf0"] = PHIf0
    elif np.isscalar(PHIf0):
        raise NotImplementedError("S.phif scalar")
    PHIf0 = COAST["PHIf0"]
    method = "weighted_distance" if COAST["cyclic"] else "alongshore_mapping"
    _, var2i, _ = interpolate_on_grid(method, COAST["xq"], COAST["yq"], PHIf0[:, 0], PHIf0[:, 1],
                                      {}, {"data": PHIf0[:, 2][None, :]})
    COAST["PHIf"] = var2i["data"][0]
    return COAST


def get_active_profile(COAST):
    COAST["h0"] = np.full(COAST["n"], float(COAST["h0input"]))
    return COAST


# ---------------------------------------------------------------------------
def coastline_change(COAST, TRANSP, TIME):
    """Move the coastline points normal to the coast by -dQS/ds / h0 * dt."""
    eps = 1e-3
    n = COAST["n"]
    nq = COAST["nq"]
    x = COAST["x"]
    y = COAST["y"]
    QS = TRANSP["QS"]
    h0 = COAST["h0"]
    if COAST["cyclic"]:
        i = np.arange(1, n + 1)
        im1 = np.array([get_mod(k - 1, n) for k in i]) - 1
        ip1 = np.array([get_mod(k + 1, n) for k in i]) - 1
        im1q = np.array([get_mod(k, nq) for k in i]) - 1
        ip1q = np.array([get_mod(k + 1, nq) for k in i]) - 1
    else:
        i = np.arange(n)
        im1 = np.maximum(i - 1, 0)
        ip1 = np.minimum(i + 1, n - 1)
        im1q = i
        ip1q = np.minimum(i + 1, nq - 1)
    ds = m_hypot(x[ip1] - x[im1], y[ip1] - y[im1]) / 2
    ds = np.fmax(ds, eps)
    dSds = (QS[ip1q] - QS[im1q]) / ds
    dndt = -dSds / h0
    rate_density = np.zeros(n)   # no nourishments
    q_tot = np.zeros(n)          # no shoreface nourishments
    SLRo = 0.0
    dn = (dndt + (rate_density / h0) - SLRo / COAST["tanbeta"] + (q_tot / h0)) * TIME["dt"]
    if COAST["preserveorientation"] == 1:
        raise NotImplementedError("S.preserveorientation")
    dx = -dn * (y[ip1] - y[im1]) / (2 * ds)
    dy = dn * (x[ip1] - x[im1]) / (2 * ds)
    if not COAST["cyclic"]:
        dx[0] = dn[0] * sind(COAST["PHIc0bnd"][0])
        dy[0] = dn[0] * cosd(COAST["PHIc0bnd"][0])
        dx[-1] = dn[-1] * sind(COAST["PHIc0bnd"][1])
        dy[-1] = dn[-1] * cosd(COAST["PHIc0bnd"][1])
    xn = x + dx
    yn = y + dy
    if COAST["cyclic"]:
        x1 = (xn[-1] + xn[0]) / 2
        y1 = (yn[-1] + yn[0]) / 2
        xn[0] = xn[-1] = x1
        yn[0] = yn[-1] = y1
    if not np.any(~np.isnan(xn)):
        xn = np.array([-1e10])
        yn = np.array([-1e10])
        dSds = np.array([-1e10])
    COAST["ds"] = ds
    COAST["dSds"] = dSds
    COAST["dndt"] = dndt
    COAST["x"] = xn
    COAST["y"] = yn
    if COAST["i_mc"] == 1:
        COAST["dSds_mc"] = dSds
    else:
        COAST["dSds_mc"] = np.concatenate([COAST["dSds_mc"], [np.nan], dSds])
    COAST["x_mc"], COAST["y_mc"] = insert_section(xn, yn, COAST["x_mc"], COAST["y_mc"], COAST["i_mc"])
    COAST = get_transportpoints(COAST, [COAST["i_mc"]])
    return COAST


# ---------------------------------------------------------------------------
def find_overwash_mc(COAST, SPIT, TIME, seg):
    """Overwash of narrow spits (find_overwash_mc.m, 'default' method).

    ``seg(i_mc)`` returns the section data (x, y, shadow, PHI, i1, cyclic)."""
    eps = 5.0
    x_mc = COAST["x_mc"]
    y_mc = COAST["y_mc"]
    dxt = np.zeros(x_mc.size)
    dyt = np.zeros(x_mc.size)
    if SPIT["method"].lower() != "default":
        raise NotImplementedError(f"spit method {SPIT['method']}")
    for ii in range(1, COAST["n_mc"] + 1):
        x, y, shadow, PHI, i1, cyclic = seg(ii)
        n = x.size
        if n < 3:
            return COAST
        phiw_coast = mmod(atan2d(sind(0.5 * (PHI[:-1] + PHI[1:])), cosd(0.5 * (PHI[:-1] + PHI[1:]))), 360.0)
        idx = np.arange(n)
        if cyclic:
            im1 = idx - 1
            im1[im1 < 0] = n - 2
            ip1 = idx + 1
            ip1[ip1 > n - 1] = 1
        else:
            im1 = np.maximum(idx - 1, 0)
            ip1 = np.minimum(idx + 1, n - 1)
        phic = 360 - atan2d(y[ip1] - y[im1], x[ip1] - x[im1])
        sp, cp = sind(phic), cosd(phic)
        x_nl = np.stack([x - eps * sp, x - SPIT["spitwidth"] * sp], axis=1)
        y_nl = np.stack([y - eps * cp, y - SPIT["spitwidth"] * cp], axis=1)
        active = ~np.asarray(shadow, dtype=bool)
        cand = np.flatnonzero(active)
        if cand.size == 0:
            continue
        hit = cand[lines_with_crossing(x_mc, y_mc, x_nl[cand], y_nl[cand])]
        for i in hit:   # few points: use the exact crossing routine, in loop order
            xcr, ycr, _, _, ind_mc, _, ui, _ = get_intersections(x_mc, y_mc, x_nl[i], y_nl[i])
            if xcr.size == 0 or np.any(np.isnan(xcr)) or np.any(np.isnan(ycr)):
                continue
            dist = m_hypot(xcr - x[i], ycr - y[i])
            closest = int(np.nanargmin(dist)) if np.any(~np.isnan(dist)) else 0
            distance = dist[closest]
            j = int(ind_mc[closest])
            weight_j = 1 - ui[closest]
            weight_jp1 = ui[closest]
            if weight_j == 1:
                jp1 = j
                weight_jp1 = 1
            else:
                jp1 = j + 1
            if SPIT["owtimescale"] > 0:
                dn_sea = -TIME["dt"] / SPIT["owtimescale"] * (SPIT["spitwidth"] - distance) * max(float(cosd(phiw_coast[i] - phic[i])), 0.0)
            else:
                dn_sea = -SPIT["owscale"] * (SPIT["spitwidth"] - distance) * max(float(cosd(phiw_coast[i] - phic[i])), 0.0)
            dn_back = -dn_sea * (SPIT["Dsf"] + SPIT["bheight"]) / (SPIT["Dbb"] + SPIT["bheight"])
            ds_sea = m_hypot(x[ip1[i]] - x[im1[i]], y[ip1[i]] - y[im1[i]])
            ddx = -dn_sea * (y[ip1[i]] - y[im1[i]]) / ds_sea
            ddy = dn_sea * (x[ip1[i]] - x[im1[i]]) / ds_sea
            k = i + i1 - 1
            dxt[k] += ddx
            dyt[k] += ddy
            j0, jp10 = j - 1, jp1 - 1
            ds_back = m_hypot(x_mc[j0 + 1] - x_mc[j0], y_mc[j0 + 1] - y_mc[j0])
            dxt[j0] = -dn_back * (y_mc[j0 + 1] - y_mc[j0]) / ds_back * weight_j
            dyt[j0] = dn_back * (x_mc[j0 + 1] - x_mc[j0]) / ds_back * weight_j
            dxt[jp10] = -dn_back * (y_mc[jp10] - y_mc[j0]) / ds_back * weight_jp1
            dyt[jp10] = dn_back * (x_mc[jp10] - x_mc[j0]) / ds_back * weight_jp1
    for ii in range(1, COAST["n_mc"] + 1):
        _, _, shadow, _, i1, cyclic = seg(ii)
        if cyclic:
            raise NotImplementedError("cyclic coastline sections")
    COAST["x_mc"] = x_mc + dxt
    COAST["y_mc"] = y_mc + dyt
    return COAST


# ---------------------------------------------------------------------------
def merge_coastlines(COAST, i_mc):
    """Remove a self-intersecting loop of section i_mc (split into 2 sections)."""
    x, y, _, _, _ = get_one_polygon(COAST["x_mc"], COAST["y_mc"], i_mc)
    x = x.copy()
    y = y.copy()
    if get_cyclic(x, y, COAST["ds0"]):
        x1 = 0.5 * (x[0] + x[-1])
        y1 = 0.5 * (y[0] + y[-1])
        x[0] = x[-1] = x1
        y[0] = y[-1] = y1
        COAST["x_mc"], COAST["y_mc"] = insert_section(x, y, COAST["x_mc"], COAST["y_mc"], i_mc)
    xx, yy, indi, indj, _, _, _, _ = get_intersections(x, y)
    xx0 = np.unique(xx)
    if xx0.size == 2:
        ind = np.sort(np.concatenate([indi[:2], indj[:2]]))
        f = lambda v: int(math.floor(v))
        c = lambda v: int(math.ceil(v))
        xnew = np.concatenate([x[:f(ind[0])], [xx[0]], x[c(ind[3]) - 1:], [np.nan],
                               [xx[1]], x[c(ind[1]) - 1:f(ind[2])], [xx[1]]])
        ynew = np.concatenate([y[:f(ind[0])], [yy[0]], y[c(ind[3]) - 1:], [np.nan],
                               [yy[1]], y[c(ind[1]) - 1:f(ind[2])], [yy[1]]])
        COAST["x_mc"], COAST["y_mc"] = insert_section(xnew, ynew, COAST["x_mc"], COAST["y_mc"], i_mc)
    COAST = get_transportpoints(COAST, [i_mc])
    return COAST


def merge_coastlines_mc(COAST):
    """Merge crossing coastline sections. Only the cases without crossings
    (or with a single crossing, which ShorelineS leaves unchanged) are ported."""
    eps = 0.1
    _, _, n_mc, _, _ = get_one_polygon(COAST["x_mc"], COAST["y_mc"], 1)
    ds0 = COAST["ds0"]
    for i_mc in range(1, n_mc):
        for j_mc in range(i_mc + 1, n_mc + 1):
            xi, yi, _, _, _ = get_one_polygon(COAST["x_mc"], COAST["y_mc"], i_mc)
            xj, yj, _, _, _ = get_one_polygon(COAST["x_mc"], COAST["y_mc"], j_mc)
            if xi.size > 2 and xj.size > 2 and (abs(msum(xi) - msum(xj)) != 0 and abs(msum(yi) - msum(yj)) != 0):
                cyclici = m_hypot(xi[-1] - xi[0], yi[-1] - yi[0]) < eps
                cyclicj = m_hypot(xj[-1] - xj[0], yj[-1] - yj[0]) < eps
                if cyclici and not cyclicj:
                    xi, xj = xj, xi
                    yi, yj = yj, yi
                    cyclici, cyclicj = cyclicj, cyclici
                cwi = get_clockpoly(xi, yi) if cyclici else 1
                cwj = get_clockpoly(xj, yj) if cyclicj else 1
                if cwi < 1 and xi[-1] != xi[0] and yi[-1] != yi[0]:
                    xi = np.append(xi, xi[0])
                    yi = np.append(yi, yi[0])
                if cwj < 1 and xj[-1] != xj[0] and yj[-1] != yj[0]:
                    xj = np.append(xj, xj[0])
                    yj = np.append(yj, yj[0])
                xcr, _, _, _, _, _, _, _ = get_intersections(xi, yi, xj, yj)
                if xcr.size > 1 and xcr.size != 3:
                    raise NotImplementedError("merging of crossing coastline sections")
                xnewi, ynewi, xnewj, ynewj = xi, yi, xj, yj
            elif xi.size <= 1 and xj.size > 2:
                xnewi, ynewi, xnewj, ynewj = [-1e10], [-1e10], xj, yj
            elif xi.size > 2 and xj.size <= 1:
                xnewi, ynewi, xnewj, ynewj = xi, yi, [-1e10], [-1e10]
            elif xi.size <= 2 and xj.size > 2:
                xnewi, ynewi, xnewj, ynewj = [-1e10], [-1e10], xj, yj
            elif xi.size > 1 and xj.size <= 2:
                xnewi, ynewi, xnewj, ynewj = xi, yi, [-1e10], [-1e10]
            elif abs(msum(xi) - msum(xj)) == 0 and abs(msum(yi) - msum(yj)) == 0:
                xnewi, ynewi, xnewj, ynewj = xi, yi, [-1e10], [-1e10]
            elif xi.size > 1 and yi.size > 1 and xj.size > 1 and yj.size > 1:
                xnewi, ynewi, xnewj, ynewj = xi, yi, xj, yj
            else:
                xnewi, ynewi, xnewj, ynewj = [-1e10], [-1e10], [-1e10], [-1e10]
            COAST["x_mc"], COAST["y_mc"] = insert_section(xnewi, ynewi, COAST["x_mc"], COAST["y_mc"], i_mc)
            COAST["x_mc"], COAST["y_mc"] = insert_section(xnewj, ynewj, COAST["x_mc"], COAST["y_mc"], j_mc)
            COAST["x_mc"], COAST["y_mc"] = get_nansremoved(COAST["x_mc"], COAST["y_mc"])
    for i_mc in range(1, n_mc + 1):
        xi, yi, _, _, _ = get_one_polygon(COAST["x_mc"], COAST["y_mc"], i_mc)
        s = 0.0
        for k in range(1, xi.size):
            s = s + m_hypot(xi[k] - xi[k - 1], yi[k] - yi[k - 1])
        if s < ds0:
            COAST["x_mc"], COAST["y_mc"] = insert_section([-1e10], [-1e10], COAST["x_mc"], COAST["y_mc"], i_mc)
    x_mc = COAST["x_mc"]
    y_mc = COAST["y_mc"]
    iddummy = np.flatnonzero(x_mc == -1e10) + 1
    for ii in range(iddummy.size - 1, -1, -1):
        d = iddummy[ii]
        if d >= x_mc.size - 1:
            x_mc, y_mc = x_mc[:d - 2], y_mc[:d - 2]
        elif d == 1:
            x_mc, y_mc = x_mc[d + 1:], y_mc[d + 1:]
        else:
            x_mc = np.concatenate([x_mc[:d - 1], x_mc[d + 1:]])
            y_mc = np.concatenate([y_mc[:d - 1], y_mc[d + 1:]])
    COAST["x_mc"], COAST["y_mc"] = x_mc, y_mc
    return COAST


def cleanup_nans(COAST):
    """Remove redundant NaNs and (nearly) duplicate points."""
    eps = 0.1
    x = list(COAST["x_mc"])
    y = list(COAST["y_mc"])
    isn = math.isnan
    i = 1
    while i <= len(x):
        L = len(x)
        if i == 1 and isn(x[0]):
            x, y = x[1:], y[1:]
        elif i == L and isn(x[i - 1]):
            x, y = x[:-1], y[:-1]
        elif 1 < i < L and isn(x[i - 2]) and isn(x[i]):
            x = x[:i - 1] + x[i + 1:]
            y = y[:i - 1] + y[i + 1:]
        elif 1 < i < L - 1 and isn(x[i - 2]) and isn(x[i + 1]):
            x = x[:i - 1] + x[i + 2:]
            y = y[:i - 1] + y[i + 2:]
        elif i == L - 1 and isn(x[i - 1]):
            x, y = x[:i - 1], y[:i - 1]
        elif isn(x[i - 1]) and isn(x[i]):
            x = x[:i] + x[i + 1:]
            y = y[:i] + y[i + 1:]
        elif i > 1 and float(m_hypot(x[i - 1] - x[i - 2], y[i - 1] - y[i - 2])) < eps:
            x = x[:i - 1] + x[i:]
            y = y[:i - 1] + y[i:]
            i = 1
        elif (i > 2 and float(m_hypot(x[i - 1] - x[i - 3], y[i - 1] - y[i - 3])) < eps
              and not isn(x[i - 3]) and not isn(x[i - 2]) and not isn(x[i - 1])):
            x = x[:i - 1] + x[i:]
            y = y[:i - 1] + y[i:]
            i = 1
        elif i == 2 and isn(x[i - 1]):
            x, y = x[i:], y[i:]
        else:
            i += 1
    COAST["x_mc"] = np.asarray(x, dtype=float)
    COAST["y_mc"] = np.asarray(y, dtype=float)
    COAST["n_mc"] = int(np.sum(np.isnan(COAST["x_mc"]))) + 1
    COAST = get_transportpoints(COAST)
    return COAST
