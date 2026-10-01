"""Polyline geometry: crossings, shadow rays and cross-shore distances.

``get_intersections`` reproduces ShorelineS' ``get_intersections.m`` exactly
(including its handling of NaN-separated polylines). The ``*_batch`` functions
evaluate many short lines against one polyline at once; they give the same
answers as calling ``get_intersections`` in a loop.
"""
import numpy as np

from .mfuncs import cosd, m_hypot, sind

EPS_CROSS = 1e-5  # tolerance used by get_intersections.m


def _pair_crossings(xi, yi, dx1, dy1, xj, yj, dx2, dy2,
                    ximin, ximax, yimin, yimax, xjmin, xjmax, yjmin, yjmax):
    """Crossing points of segments i and j (any broadcastable shapes, e.g.
    column vectors for i and row vectors for j, or flat arrays of pairs).

    Returns xc, yc (NaN where the segments do not cross), identical to the
    arithmetic in get_intersections.m."""
    (xi, yi, dx1, dy1, xj, yj, dx2, dy2, ximin, ximax, yimin, yimax,
     xjmin, xjmax, yjmin, yjmax) = np.broadcast_arrays(xi, yi, dx1, dy1, xj, yj, dx2, dy2, ximin, ximax,
                                                       yimin, yimax, xjmin, xjmax, yjmin, yjmax)
    with np.errstate(divide="ignore", invalid="ignore"):
        rc1 = dy1 / dx1
        rc2 = dy2 / dx2
        y1r = yi - xi * rc1
        y2r = yj - xj * rc2
        xc = (y2r - y1r) / (rc1 - rc2)
        yc = rc1 * xc + y1r
        both = (dx1 != 0) & (dx2 != 0)
        xc = np.where(both, xc, np.nan)
        yc = np.where(both, yc, np.nan)
        id2 = (dx1 == 0) & (dx2 != 0)
        if np.any(id2):
            xc = np.where(id2, xi, xc)
            yc = np.where(id2, rc2 * xi + y2r, yc)
        id3 = (dx1 != 0) & (dx2 == 0)
        if np.any(id3):
            xc = np.where(id3, xj, xc)
            yc = np.where(id3, rc1 * xj + y1r, yc)
        out = (xc < np.fmax(ximin, xjmin) - EPS_CROSS) | (xc > np.fmin(ximax, xjmax) + EPS_CROSS)
        xc = np.where(out, np.nan, xc)
        yc = np.where(out, np.nan, yc)
        out = (yc < np.fmax(yimin, yjmin) - EPS_CROSS) | (yc > np.fmin(yimax, yjmax) + EPS_CROSS)
        xc = np.where(out, np.nan, xc)
        yc = np.where(out, np.nan, yc)
    return xc, yc


def _seg_arrays(x, y):
    """Start point, increments and bounding box of the segments of a polyline."""
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    return dict(x=x[:-1], y=y[:-1], dx=x[1:] - x[:-1], dy=y[1:] - y[:-1],
                xmin=np.fmin(x[:-1], x[1:]), xmax=np.fmax(x[:-1], x[1:]),
                ymin=np.fmin(y[:-1], y[1:]), ymax=np.fmax(y[:-1], y[1:]),
                nan=np.isnan(x[:-1] + x[1:] + y[:-1] + y[1:]))


def _candidate_pairs(A, B):
    """Index pairs (a, b) of segments whose bounding boxes overlap (with a
    tolerance larger than the crossing tolerance), plus all pairs that involve
    a NaN-adjacent segment (which get_intersections.m can report as crossings)."""
    tol = 4 * EPS_CROSS
    with np.errstate(invalid="ignore"):
        ov = ((A["xmin"][:, None] <= B["xmax"][None, :] + tol) & (A["xmax"][:, None] >= B["xmin"][None, :] - tol)
              & (A["ymin"][:, None] <= B["ymax"][None, :] + tol) & (A["ymax"][:, None] >= B["ymin"][None, :] - tol))
    ov |= A["nan"][:, None] | B["nan"][None, :]
    return np.nonzero(ov)


def _crossings_of_pairs(A, B, ia, ib):
    return _pair_crossings(A["x"][ia], A["y"][ia], A["dx"][ia], A["dy"][ia],
                           B["x"][ib], B["y"][ib], B["dx"][ib], B["dy"][ib],
                           A["xmin"][ia], A["xmax"][ia], A["ymin"][ia], A["ymax"][ia],
                           B["xmin"][ib], B["xmax"][ib], B["ymin"][ib], B["ymax"][ib])


def get_intersections(xi0, yi0, xj0=None, yj0=None):
    """All crossings of polyline i with polyline j (or with itself).

    Returns xcr, ycr, indc, inds, indi, indj, ui, uj as in get_intersections.m
    (indices 1-based)."""
    removesameindex = xj0 is None
    xi0 = np.asarray(xi0, dtype=float).ravel()
    yi0 = np.asarray(yi0, dtype=float).ravel()
    if removesameindex:
        xj0, yj0 = xi0, yi0
    else:
        xj0 = np.asarray(xj0, dtype=float).ravel()
        yj0 = np.asarray(yj0, dtype=float).ravel()
    m0 = xi0.size - 1
    n0 = xj0.size - 1
    empty = np.zeros(0)
    if m0 < 1 or n0 < 1:
        return empty, empty, empty, empty, empty, empty, empty, empty

    ximin_a = np.fmin(xi0[:-1], xi0[1:])
    ximax_a = np.fmax(xi0[:-1], xi0[1:])
    yimin_a = np.fmin(yi0[:-1], yi0[1:])
    yimax_a = np.fmax(yi0[:-1], yi0[1:])
    xjmin_a = np.fmin(xj0[:-1], xj0[1:])
    xjmax_a = np.fmax(xj0[:-1], xj0[1:])
    yjmin_a = np.fmin(yj0[:-1], yj0[1:])
    yjmax_a = np.fmax(yj0[:-1], yj0[1:])

    ri = np.arange(m0)
    cj = np.arange(n0)
    if not removesameindex:
        # bounding-box prefilter (same as the optimised get_intersections.m)
        tol = 4 * EPS_CROSS
        nani = np.isnan(xi0[:-1] + xi0[1:] + yi0[:-1] + yi0[1:])
        nanj = np.isnan(xj0[:-1] + xj0[1:] + yj0[:-1] + yj0[1:])
        nani = nani | (np.any(nanj) & ((xi0[1:] - xi0[:-1]) == 0))
        nanj = nanj | (np.any(nani) & ((xj0[1:] - xj0[:-1]) == 0))
        with np.errstate(invalid="ignore"):
            ri = np.flatnonzero(nani | ((ximin_a <= np.nanmax(xjmax_a) + tol) & (ximax_a >= np.nanmin(xjmin_a) - tol)
                                        & (yimin_a <= np.nanmax(yjmax_a) + tol) & (yimax_a >= np.nanmin(yjmin_a) - tol)))
            if ri.size == 0:
                cj = np.zeros(0, dtype=int)
            else:
                cj = np.flatnonzero(nanj | ((xjmin_a <= np.nanmax(ximax_a[ri]) + tol) & (xjmax_a >= np.nanmin(ximin_a[ri]) - tol)
                                            & (yjmin_a <= np.nanmax(yimax_a[ri]) + tol) & (yjmax_a >= np.nanmin(yimin_a[ri]) - tol)))
    if ri.size == 0 or cj.size == 0:
        return empty, empty, empty, empty, empty, empty, empty, empty

    xi = xi0[ri][:, None]
    yi = yi0[ri][:, None]
    xj = xj0[cj][None, :]
    yj = yj0[cj][None, :]
    dx1 = xi0[ri + 1][:, None] - xi
    dy1 = yi0[ri + 1][:, None] - yi
    dx2 = xj0[cj + 1][None, :] - xj
    dy2 = yj0[cj + 1][None, :] - yj
    xc, yc = _pair_crossings(xi, yi, dx1, dy1, xj, yj, dx2, dy2,
                             ximin_a[ri][:, None], ximax_a[ri][:, None], yimin_a[ri][:, None], yimax_a[ri][:, None],
                             xjmin_a[cj][None, :], xjmax_a[cj][None, :], yjmin_a[cj][None, :], yjmax_a[cj][None, :])

    if removesameindex:
        flat_x = xc.ravel(order="F")
        flat_y = yc.ravel(order="F")
        mn = m0 * n0
        for start in (1, 2, m0 + 1):
            idx = np.arange(start, mn + 1, m0 + 1) - 1
            flat_x[idx] = np.nan
            flat_y[idx] = np.nan
        xc = flat_x.reshape((m0, n0), order="F")
        yc = flat_y.reshape((m0, n0), order="F")
        # only crossings of the first line, not the mirrored ones
        I, J = np.meshgrid(np.arange(m0), np.arange(n0), indexing="ij")
        # mat0=meshgrid(1:m0,1:n0) has mat0(r,c)=c, so idmat=mat0>mat0' is c > r
        idmat = J > I
        xc[idmat] = np.nan
        yc[idmat] = np.nan
        if xi0[0] == xi0[-1] and yi0[0] == yi0[-1] and xc.size:
            xc[-1, 0] = np.nan
            yc[-1, 0] = np.nan

    # column-major order of the crossings, as find() in MATLAB
    sel = ~np.isnan(xc)
    jj, ii = np.nonzero(sel.T)
    xcr = xc[ii, jj]
    ycr = yc[ii, jj]
    indi = ri[ii]          # 0-based segment index
    indj = cj[jj]
    with np.errstate(divide="ignore", invalid="ignore"):
        dxi = xi0[indi + 1] - xi0[indi]
        dyi = yi0[indi + 1] - yi0[indi]
        dxj = xj0[indj + 1] - xj0[indj]
        dyj = yj0[indj + 1] - yj0[indj]
        ui = ((xcr - xi0[indi]) * dxi + (ycr - yi0[indi]) * dyi) / (dxi * dxi + dyi * dyi)
        uj = ((xcr - xj0[indj]) * dxj + (ycr - yj0[indj]) * dyj) / (dxj * dxj + dyj * dyj)
    ui = np.fmin(np.fmax(ui, 0.0), 1.0)
    uj = np.fmin(np.fmax(uj, 0.0), 1.0)
    indi1 = indi + 1.0
    indj1 = indj + 1.0
    _, idu = np.unique(indi1 + ui, return_index=True)
    xcr, ycr, indi1, indj1, ui, uj = xcr[idu], ycr[idu], indi1[idu], indj1[idu], ui[idu], uj[idu]
    _, idu = np.unique(indj1 + uj, return_index=True)
    xcr, ycr, indi1, indj1, ui, uj = xcr[idu], ycr[idu], indi1[idu], indj1[idu], ui[idu], uj[idu]
    return xcr, ycr, indi1 + ui, indj1 + uj, indi1, indj1, ui, uj


def any_crossing_batch(x_mc, y_mc, xl, yl):
    """For each polyline row k of (xl, yl) (shape [K, P]), whether it crosses
    polyline (x_mc, y_mc). Equivalent to ``~isempty(get_intersections(x_mc,y_mc,xl(k,:),yl(k,:)))``."""
    A = _seg_arrays(x_mc, y_mc)
    K, P = xl.shape
    hit = np.zeros(K, dtype=bool)
    for p in range(P - 1):
        x1, x2, y1, y2 = xl[:, p], xl[:, p + 1], yl[:, p], yl[:, p + 1]
        B = dict(x=x1, y=y1, dx=x2 - x1, dy=y2 - y1, xmin=np.fmin(x1, x2), xmax=np.fmax(x1, x2),
                 ymin=np.fmin(y1, y2), ymax=np.fmax(y1, y2), nan=np.isnan(x1 + x2 + y1 + y2))
        ia, ib = _candidate_pairs(A, B)
        if ia.size == 0:
            continue
        xc, _ = _crossings_of_pairs(A, B, ia, ib)
        hit[ib[~np.isnan(xc)]] = True
    return hit


def lines_with_crossing(x_mc, y_mc, x_nl, y_nl):
    """Indices k of the 2-point lines (rows of x_nl, y_nl) that cross polyline
    (x_mc, y_mc), i.e. for which get_intersections(x_mc,y_mc,x_nl(k,:),y_nl(k,:))
    is not empty."""
    return np.flatnonzero(any_crossing_batch(x_mc, y_mc, np.asarray(x_nl), np.asarray(y_nl)))


def get_polydistance(Xr, Yr, Xc, Yc, Lcrit=500.0, nargout=3):
    """Cross-shore distance from reference line (Xr, Yr) to coastline (Xc, Yc)
    along normals of +-Lcrit (vectorised port of get_polydistance.m)."""
    Xr = np.asarray(Xr, dtype=float).ravel()
    Yr = np.asarray(Yr, dtype=float).ravel()
    Xc = np.asarray(Xc, dtype=float).ravel()
    Yc = np.asarray(Yc, dtype=float).ravel()
    nr = Xr.size
    dx = Xr[1:] - Xr[:-1]
    dy = Yr[1:] - Yr[:-1]
    dx2 = np.concatenate([dx[:1], (dx[1:] + dx[:-1]) / 2, dx[-1:]])
    dy2 = np.concatenate([dy[:1], (dy[1:] + dy[:-1]) / 2, dy[-1:]])
    xcr = np.full(nr, np.nan)
    ycr = np.full(nr, np.nan)
    dmin = np.full(nr, np.nan)
    xcm = np.full(nr, np.nan)
    ycm = np.full(nr, np.nan)
    dmax = np.full(nr, np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        Lxy = np.sqrt(dx2 * dx2 + dy2 * dy2)
        Xn = np.stack([Xr - Lcrit / Lxy * dy2, Xr + Lcrit / Lxy * dy2], axis=1)
        Yn = np.stack([Yr + Lcrit / Lxy * dx2, Yr - Lcrit / Lxy * dx2], axis=1)
    nseg = Xc.size - 1
    if nseg >= 1:
        xj = Xc[:-1][None, :]
        yj = Yc[:-1][None, :]
        dxc = (Xc[1:] - Xc[:-1])[None, :]
        dyc = (Yc[1:] - Yc[:-1])[None, :]
        jxmin = np.fmin(Xc[:-1], Xc[1:])[None, :]
        jxmax = np.fmax(Xc[:-1], Xc[1:])[None, :]
        jymin = np.fmin(Yc[:-1], Yc[1:])[None, :]
        jymax = np.fmax(Yc[:-1], Yc[1:])[None, :]
        nchunk = max(1, int(2e6 // nseg))
        for i0 in range(0, nr, nchunk):
            ii = np.arange(i0, min(i0 + nchunk, nr))
            xi = Xn[ii, 0][:, None]
            yi = Yn[ii, 0][:, None]
            dx1 = Xn[ii, 1][:, None] - xi
            dy1 = Yn[ii, 1][:, None] - yi
            xc, yc = _pair_crossings(xi, yi, dx1, dy1, xj, yj, dxc, dyc,
                                     np.fmin(xi, Xn[ii, 1][:, None]), np.fmax(xi, Xn[ii, 1][:, None]),
                                     np.fmin(yi, Yn[ii, 1][:, None]), np.fmax(yi, Yn[ii, 1][:, None]),
                                     jxmin, jxmax, jymin, jymax)
            valid = ~np.isnan(xc)
            with np.errstate(invalid="ignore", divide="ignore"):
                ui = ((xc - xi) * dx1 + (yc - yi) * dy1) / (dx1 * dx1 + dy1 * dy1)
            ui = np.fmin(np.fmax(ui, 0.0), 1.0)
            umin_a = np.where(valid, ui, np.inf)
            jmin = np.argmin(umin_a, axis=1)
            umin = umin_a[np.arange(ii.size), jmin]
            has = np.isfinite(umin)
            rows = np.flatnonzero(has)
            xcr[ii[has]] = xc[rows, jmin[has]]
            ycr[ii[has]] = yc[rows, jmin[has]]
            dmin[ii[has]] = (0.5 - umin[has]) * Lcrit * 2
            if nargout > 3:
                umax_a = np.where(valid, ui, -np.inf)
                jmax = np.argmax(umax_a, axis=1)
                umax = umax_a[np.arange(ii.size), jmax]
                xcm[ii[has]] = xc[rows, jmax[has]]
                ycm[ii[has]] = yc[rows, jmax[has]]
                dmax[ii[has]] = (0.5 - umax[has]) * Lcrit * 2
    if nargout > 3:
        return dmin, xcr, ycr, dmax, xcm, ycm
    return dmin, xcr, ycr


def find_shadows_mc(xq, yq, x_mc, y_mc, PHI, PHItdp, PHIbr, distw):
    """Which transport points are in the shadow of the coastline (port of
    find_shadows_mc.m, returning only shadowS). Each point casts a 3-segment
    ray (breaker zone, surf zone, offshore) towards the waves."""
    xq = np.asarray(xq, dtype=float)
    yq = np.asarray(yq, dtype=float)
    nq = xq.size
    if nq == 0:
        return np.zeros(0, dtype=bool)
    PHI = np.broadcast_to(np.asarray(PHI, dtype=float), (nq,))
    lenSURF = np.full(nq, 250.0)
    lenBREAKER = np.full(nq, 50.0)
    length = np.full(nq, 5 * m_hypot(np.max(xq) - np.min(xq), np.max(yq) - np.min(yq)))
    if distw is not None and np.size(distw):
        distw = np.asarray(distw, dtype=float)
        lenBREAKER = np.minimum(lenBREAKER, distw * 0.2)
        lenSURF = np.minimum(lenSURF, distw * 0.8)
        length = np.maximum(distw - lenBREAKER - lenSURF, 1)
    sbr, cbr = sind(PHIbr), cosd(PHIbr)
    stdp, ctdp = sind(PHItdp), cosd(PHItdp)
    so, co = sind(PHI), cosd(PHI)
    xw = np.stack([xq + 1 * sbr,
                   xq + lenBREAKER * sbr,
                   xq + lenBREAKER * sbr + lenSURF * stdp,
                   xq + lenBREAKER * sbr + lenSURF * stdp + length * so], axis=1)
    yw = np.stack([yq + 1 * cbr,
                   yq + lenBREAKER * cbr,
                   yq + lenBREAKER * cbr + lenSURF * ctdp,
                   yq + lenBREAKER * cbr + lenSURF * ctdp + length * co], axis=1)
    return any_crossing_batch(x_mc, y_mc, xw, yw)
