"""Merging of crossing coastline sections: the part of merge_coastlines_mc.m
that handles sections crossing each other (e.g. two beaches growing into each
other, or an island growing into the coast).

The code follows the MATLAB function statement by statement, including its
quirks, so that results can be compared with ShorelineS exactly:

- indices are 1-based (as in MATLAB) and converted only when slicing;
- ``eps`` is doubled for the rest of the call when tracing a section fails
  (``ff > 50``), and ``IDcyclic`` keeps entries from earlier section pairs;
  both are kept in a per-call state dict;
- after the tracing loop, ``xi``/``yi``/``xj``/``yj`` are the sections with the
  crossing points inserted (MATLAB reuses the variable names).
"""
import math

import numpy as np

from .mfuncs import get_clockpoly, m_hypot, m_pow, mmod, msum

_DUMMY = -1e10


def inpolygon(x, y, xv, yv):
    """GNU Octave's inpolygon (points on the boundary count as inside)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xv = np.asarray(xv, dtype=float).ravel()
    yv = np.asarray(yv, dtype=float).ravel()
    inside = np.zeros(x.shape, dtype=bool)
    on = np.zeros(x.shape, dtype=bool)
    npol = xv.size
    j = npol - 1
    with np.errstate(invalid="ignore"):
        for i in range(npol):
            dxv = xv[j] - xv[i]
            dyv = yv[j] - yv[i]
            distance = dxv * (y - yv[i]) - (x - xv[i]) * dyv
            idx1 = ((((yv[i] <= y) & (y < yv[j])) | ((yv[j] <= y) & (y < yv[i])))
                    & (0 < distance * dyv))
            inside[idx1] = ~inside[idx1]
            idx2 = ((((yv[i] <= y) & (y <= yv[j])) | ((yv[j] <= y) & (y <= yv[i])))
                    & (((xv[i] <= x) & (x <= xv[j])) | ((xv[j] <= x) & (x <= xv[i])))
                    & ((distance == 0) | (dxv == 0)))
            on[idx2] = True
            j = i
    return inside | on


def _frac_inside(xp, yp, xv, yv):
    """sum(inpolygon(xp, yp, xv, yv)) / length(xp) (NaN for an empty xp)."""
    n = np.size(xp)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.float64(np.sum(inpolygon(xp, yp, xv, yv))) / np.float64(n)


def _rowmax(a):
    """MATLAB max(a, [], 2): ignores NaN; NaN only if the whole row is NaN."""
    a = np.atleast_2d(np.asarray(a, dtype=float))
    out = np.full(a.shape[0], np.nan)
    for r in range(a.shape[0]):
        v = a[r][~np.isnan(a[r])]
        if v.size:
            out[r] = v.max()
    return out


def _find1(mask):
    """MATLAB find(mask) (1-based)."""
    return [int(k) + 1 for k in np.flatnonzero(mask)]


def _setdiff(a, b):
    """MATLAB setdiff (sorted, unique)."""
    return sorted(set(a) - set(b))


def _dist(xa, ya, xb, yb):
    """((xa-xb).^2+(ya-yb).^2).^0.5"""
    dx = xa - xb
    dy = ya - yb
    return m_pow(dx * dx + dy * dy, 0.5)


def merge_crossing_sections(st, xi, yi, xj, yj, cyclici, cyclicj, cwi, cwj, xcr, ycr, indi, indj):
    """Body of merge_coastlines_mc.m for one pair of sections with more than one
    crossing (and not 3). Returns (xnewi, ynewi, xnewj, ynewj); a section equal
    to [-1e10] is removed by the caller. `st` holds eps and IDcyclic."""
    xcr = np.asarray(xcr, dtype=float).ravel()
    ycr = np.asarray(ycr, dtype=float).ravel()
    indi_f = np.asarray(indi, dtype=float).ravel()
    indj_f = np.asarray(indj, dtype=float).ravel()

    # fraction of the vertex at which the crossings take place
    indi = np.fix(indi_f).astype(np.int64)
    umin = mmod(np.abs(indi_f), 1)
    indj = np.fix(indj_f).astype(np.int64)
    vmin = mmod(np.abs(indj_f), 1)
    P = np.vstack([xcr, ycr])
    ipi = list(range(1, indi.size + 1))
    IDS = []
    id2 = st.get("id2")                              # a function-level variable in MATLAB

    if indi.size:
        # add the crossing points to both sections
        xi2, yi2, xj2, yj2 = xi.copy(), yi.copy(), xj.copy(), yj.copy()
        idi2 = np.argsort(indi + umin, kind="stable")
        idj2 = np.argsort(indj + vmin, kind="stable")
        rr = 0
        for i in idi2:
            k = indi[i] + rr
            xi2 = np.concatenate([xi2[:k], [P[0, i]], xi2[k:]])
            yi2 = np.concatenate([yi2[:k], [P[1, i]], yi2[k:]])
            rr += 1
            indi[i] = indi[i] + rr
        rr = 0
        for j in idj2:
            k = indj[j] + rr
            xj2 = np.concatenate([xj2[:k], [xcr[j]], xj2[k:]])
            yj2 = np.concatenate([yj2[:k], [ycr[j]], yj2[k:]])
            rr += 1
            indj[j] = indj[j] + rr

        if cyclicj or cyclici:
            o = np.argsort(indi, kind="stable")
            indi = indi[o]
            indj = indj[o]
            P = P[:, o]

        # Trace the new sections: start from each crossing, follow one section in
        # positive direction to the next crossing, switch to the other section, ...
        xnew2, ynew2 = [], []
        ids2 = np.zeros(indi.size * indi.size)       # zeros(length(indi)): n x n, linear indexing
        for pp in range(1, indi.size + 1):
            xi, yi, xj, yj = xi2, yi2, xj2, yj2
            id1 = ipi[pp - 1]
            ids = []
            kk = 1
            kkstart = math.nan
            if not cyclicj and cyclici:
                kk = 0
            xnew = np.zeros(0)
            ynew = np.zeros(0)
            circle = False
            ff = 0
            while not circle:
                if mmod(kk, 2) == 1:
                    indk0 = indi[id1 - 1]
                    indk, xk, yk, cyclick = indi, xi, yi, cyclici
                else:
                    indk0 = indj[id1 - 1]
                    indk, xk, yk, cyclick = indj, xj, yj, cyclicj
                if not ids:
                    id1start = id1
                # start and end at xk(1) and xk(end) when both sections are open
                if not cyclici and not cyclicj:
                    if not ids and pp == 1:
                        indk0 = 1
                        kkstart = mmod(kk, 2)

                if indk.size:
                    larger = indk[indk > indk0]
                    smaller = indk[indk < indk0]
                    if larger.size:
                        # next stretch of xk (positive direction) up to the next crossing
                        id2 = int(np.flatnonzero(indk == larger.min())[0]) + 1
                        seg = slice(indk0 - 1, indk[id2 - 1])
                        if id2 <= P.shape[1]:
                            xnew = np.concatenate([xnew, xk[seg], [P[0, id2 - 1]]])
                            ynew = np.concatenate([ynew, yk[seg], [P[1, id2 - 1]]])
                        else:
                            xnew = np.concatenate([xnew, xk[seg]])
                            ynew = np.concatenate([ynew, yk[seg]])
                    elif smaller.size:
                        # last stretch of xk and first stretch up to the first crossing
                        id2 = int(np.flatnonzero(indk == smaller.min())[0]) + 1
                        a_x, a_y = xk[indk0:], yk[indk0:]
                        b_x, b_y = xk[:indk[id2 - 1]], yk[:indk[id2 - 1]]
                        if not cyclick and mmod(kk, 2) == kkstart:
                            xnew = np.concatenate([xnew, a_x])
                            ynew = np.concatenate([ynew, a_y])
                        elif not cyclick:
                            xnew = np.concatenate([xnew, a_x, [np.nan], b_x, [xcr[id2 - 1]]])
                            ynew = np.concatenate([ynew, a_y, [np.nan], b_y, [ycr[id2 - 1]]])
                        else:
                            xnew = np.concatenate([xnew, a_x, b_x, [xcr[id2 - 1]]])
                            ynew = np.concatenate([ynew, a_y, b_y, [ycr[id2 - 1]]])
                    else:
                        print("Special case for merging (e.g. with only 1 intersection point)")
                        indi = np.zeros(0, dtype=np.int64)
                        id2 = 1

                # next crossing, other section
                id1 = int(mmod(id1, indi.size)) + 1 if indi.size else id1 + 1   # mod(x, 0) = x
                kk += 1

                if (ids and np.any(ids2[np.asarray(ids) - 1] == 2)) or xnew.size == 0:
                    xnew = np.array([_DUMMY])
                    ynew = np.array([_DUMMY])
                    circle = True
                # stopping criterion ('isempty(id2) & ...' can never hold: id2 is set above)
                if xnew.size < 2:
                    circle = True
                    ids.append(id2)
                elif (_dist(xnew[0], ynew[0], xnew[-1], ynew[-1]) < st["eps"]
                      or _dist(xnew[-1], ynew[-1], xj[-1], yj[-1]) < st["eps"]
                      or _dist(xnew[-1], ynew[-1], xi[-1], yi[-1]) < st["eps"]
                      or (_dist(xnew[-1], ynew[-1], xcr[id1start - 1], ycr[id1start - 1]) < st["eps"]
                          and (cyclici or cyclicj))):
                    circle = True
                    if ids:
                        ids2[np.asarray(ids) - 1] += 1
                else:
                    ids.append(id2)
                ff += 1

                # tracing failed: redo with a larger eps (which stays for the rest of the call)
                if ff > 50:
                    st["eps"] = 2 * st["eps"]
                    cyclici = bool(m_hypot(xi[-1] - xi[0], yi[-1] - yi[0]) < st["eps"])
                    cyclicj = bool(m_hypot(xj[-1] - xj[0], yj[-1] - yj[0]) < st["eps"])
                    xi, yi, xj, yj = xi2, yi2, xj2, yj2
                    id1 = ipi[pp - 1]
                    ids = []
                    kk = 1
                    kkstart = math.nan
                    if not cyclicj and cyclici:
                        kk = 0
                    xnew = np.zeros(0)
                    ynew = np.zeros(0)
                    circle = False
                    ff = 0

            # join the two halves of an open section into one line
            nanpos = np.flatnonzero(np.isnan(xnew))
            if not cyclici and not cyclicj and nanpos.size:
                pass                                 # no action on NaNs if 2 open sections are used
            else:
                while True:
                    nanpos = np.flatnonzero(np.isnan(xnew))
                    if not nanpos.size:
                        break
                    k = nanpos[-1]
                    xnew = np.concatenate([xnew[k + 1:], xnew[:k]])
                    ynew = np.concatenate([ynew[k + 1:], ynew[:k]])
                if cyclick and not (not cyclici or not cyclicj):
                    xnew = np.append(xnew, xnew[0])
                    ynew = np.append(ynew, ynew[0])

            # remove sequential duplicate points
            dx = np.diff(xnew)
            dy = np.diff(ynew)
            keep = np.concatenate([[True], (dx * dx + dy * dy) != 0])
            xnew = xnew[keep]
            ynew = ynew[keep]

            # split open sections at the remaining NaN
            idnan = _find1(np.isnan(xnew))
            idele = [0] + idnan + [xnew.size + 1]
            for zz in range(len(idnan) + 1):
                if xnew[0] != _DUMMY:
                    xnew2.append(xnew[idele[zz]:idele[zz + 1] - 1])
                    ynew2.append(ynew[idele[zz]:idele[zz + 1] - 1])

        # ---- which of the traced sections to keep
        IDclock = [get_clockpoly(xnew2[g], ynew2[g]) for g in range(len(xnew2))]
        IDcw = [g + 1 for g, c in enumerate(IDclock) if c >= 0]
        IDccw = [g + 1 for g, c in enumerate(IDclock) if c < 0]

        IDcyclic = st["IDcyclic"]                    # persists over section pairs (as in MATLAB)
        IDopen = []
        for pp in range(1, len(xnew2) + 1):
            while len(IDcyclic) < pp:
                IDcyclic.append(0)
            IDcyclic[pp - 1] = 1
            if xnew2[pp - 1].size > 1:
                v = xnew2[pp - 1]
                w = ynew2[pp - 1]
                IDcyclic[pp - 1] = int(m_hypot(v[-1] - v[0], w[-1] - w[0]) < st["eps"])
            IDopen = [k + 1 for k, c in enumerate(IDcyclic) if c == 0]

        if IDopen:
            # two open sections (and rare line + island/lake cases)
            IDS = list(IDopen)
            if (len(IDopen) == 1 and IDccw and cwi != -1 and cwj != -1 and (not cyclici or not cyclicj)
                    and (cyclici or cyclicj)):
                IDS = IDopen + _setdiff(IDccw, IDopen)
            elif len(IDopen) == 1 and IDccw and cwi != -1 and not cyclici and (cwj == -1 and cyclicj):
                IDS = IDopen + _setdiff(IDcw, IDopen)
            elif len(IDopen) == 1 and IDccw and cwj != -1 and not cyclicj and (cwi == -1 and cyclici):
                IDS = IDopen + _setdiff(IDcw, IDopen)
            elif len(IDopen) == 1 and not cyclici and (cwj == -1 and cyclicj):
                IDS = IDopen + _setdiff(IDcw, IDopen)
        else:
            # closed sections: polygons inside other new polygons
            n2 = len(xnew2)
            inpoly1 = np.zeros((n2, n2))
            for pp in range(n2):
                for gg in range(n2):
                    if pp != gg:
                        inpoly1[pp, gg] = _frac_inside(xnew2[pp], ynew2[pp], xnew2[gg], ynew2[gg])
                        if ((msum(xnew2[pp]) - msum(xnew2[gg]) == 0 and msum(ynew2[pp]) - msum(ynew2[gg]) == 0
                             and pp < gg) or xnew2[pp].size <= 1 or np.sum(np.isnan(xnew2[gg])) > 0):
                            inpoly1[pp, gg] = 0          # keep the first of two identical polygons
            IDinsidenew = _find1(_rowmax(inpoly1) > 0.9999999) if n2 else []
            # inside one of the original sections
            inpoly2 = np.zeros((n2, 2))
            for pp in range(n2):
                inpoly2[pp, 0] = _frac_inside(xnew2[pp], ynew2[pp], xi, yi)
                inpoly2[pp, 1] = _frac_inside(xnew2[pp], ynew2[pp], xj, yj)
            IDinsideorignal = _find1(_rowmax(inpoly2) > 0.9999999) if n2 else []

            # outer contour (not inside another contour)
            if cwi >= 0 and cwj >= 0:
                IDouter = _setdiff(IDcw, IDinsideorignal)
                if not IDouter and indi.size:
                    L1L2 = np.array([xnew2[0].size, xnew2[1].size])
                    IDouter = _setdiff(IDcw, _find1(L1L2 == L1L2.min()))
                elif not IDouter:
                    IDouter = list(IDcw)
            else:
                IDouter = list(IDcw)
            # inner contours (counter-clockwise): lakes
            if cwi >= 0 and cwj >= 0:
                IDinner = list(IDccw)
            elif cwi >= 0 and cwj == -1:
                IDinner = []
            elif cwi == -1 and cwj >= 0:
                IDinner = []
            else:
                IDinner = _setdiff(IDccw, IDinsidenew)
            IDS = list(IDouter) + list(IDinner)

            # a small closed section crossing an open one twice: keep the largest pieces
            if len(IDS) > 2:
                if len(IDinner) > 1:
                    len1 = np.array([xnew2[k - 1].size for k in IDinner])
                    IDinner = [IDinner[int(np.flatnonzero(len1 == len1.max())[0])]]
                if len(IDouter) >= 2:
                    len1 = np.array([m_hypot(xnew2[k - 1][-1] - xnew2[k - 1][0],
                                             ynew2[k - 1][-1] - ynew2[k - 1][0]) for k in IDouter])
                    IDouter = [IDouter[int(np.flatnonzero(len1 == len1.max())[0])]]
                IDS = list(IDouter) + list(IDinner)

    st["id2"] = id2

    # ---- output
    if not indi.size:
        return xi, yi, xj, yj                        # nothing happened
    if len(IDS) == 1:
        return xnew2[IDS[0] - 1], ynew2[IDS[0] - 1], np.array([_DUMMY]), np.array([_DUMMY])
    if len(IDS) == 2:
        if IDS[1] == IDS[0]:
            IDS[1] = int(mmod(IDS[1], 2)) + 1
        return xnew2[IDS[0] - 1], ynew2[IDS[0] - 1], xnew2[IDS[1] - 1], ynew2[IDS[1] - 1]
    print("exceptional case")
    return xi, yi, xj, yj
