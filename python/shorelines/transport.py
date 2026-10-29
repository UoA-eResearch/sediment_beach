"""Longshore sediment transport and its corrections. Ports of transport.m
(CERC, CERC2, CERC3, KAMP, MILH, VR14), transport_shadow_treat.m, get_smoothangles.m,
get_upwindcorrection.m and transport_boundary_condition.m."""
import numpy as np

from .mfuncs import cosd, m_pow, mmod, sind


TRANSPORT_FORMULAS = ("CERC", "CERC1", "CERC2", "CERC3", "KAMP", "MILH", "VR14")


def transport(TRANSP, WAVE, field="QS"):
    """Bulk longshore transport [m3/yr] at the transport points (transport.m):
    CERC/CERC1 (simple CERC with S.b), CERC2 (CERC with implicit refraction),
    CERC3 (CERC with breaking waves, Vitousek & Barnard 2015), KAMP (Kamphuis
    1991), MILH (Mil-Homens et al. 2013) or VR14 (Van Rijn 2014).
    `field` is 'QSmax' when called for the critical angle (no suppresshighangle factor)."""
    T = TRANSP
    trform = T["trform"].upper()
    if trform not in TRANSPORT_FORMULAS:
        raise NotImplementedError(f"transport formula {T['trform']} (ported: {', '.join(TRANSPORT_FORMULAS)})")
    HStdp = np.asarray(WAVE["HStdp"], dtype=float)
    dHS = np.zeros(HStdp.size)                       # only non-zero with diffraction (not ported)
    # (constants multiplied in MATLAB's left-to-right order, for identical rounding)
    with np.errstate(invalid="ignore", divide="ignore"):
        if trform in ("CERC", "CERC1"):
            dPHItdp = WAVE["dPHItdp"]
            QS = T["qscal"] * T["b"] * m_pow(HStdp, 2.5) * (sind(2 * dPHItdp) - 2 * cosd(dPHItdp) * dHS)
            QS[np.abs(dPHItdp) > 90] = 0
        elif trform == "CERC2":
            dPHItdp = WAVE["dPHItdp"]
            k = 0.39
            b1 = k * (T["rhow"] * m_pow(T["g"], 0.5) / (16 * np.sqrt(T["gamma"]) * (T["rhos"] - T["rhow"])
                                                        * (1 - T["porosity"])))
            b2 = b1 * m_pow(m_pow(T["gamma"] * T["g"], 0.5) / (2 * np.pi), 0.2)
            QS = (T["qscal"] * 365 * 24 * 60 * 60 * b2 * m_pow(HStdp, 12 / 5) * m_pow(WAVE["TP"], 1 / 5)
                  * (m_pow(cosd(dPHItdp), 6 / 5) * sind(dPHItdp) - 2 * cosd(dPHItdp) * dHS))
            QS[np.abs(dPHItdp) > 90] = 0
        elif trform == "CERC3":
            HSbr, dPHIbr = WAVE["HSbr"], WAVE["dPHIbr"]
            k = 0.35
            b3 = 1 / 16 * k * T["rhow"] / ((T["rhos"] - T["rhow"]) * (1 - T["porosity"]))
            QS = (T["qscal"] * 365 * 24 * 60 * 60 * b3 * m_pow(T["g"], 0.5) * m_pow(T["gamma"], -0.52) * m_pow(HSbr, 2.5)
                  * (sind(2 * dPHIbr) - 2 * cosd(dPHIbr) * dHS))
            QS[np.abs(dPHIbr) > 90] = 0
        elif trform == "KAMP":
            HSbr, dPHIbr, TP = WAVE["HSbr"], WAVE["dPHIbr"], WAVE["TP"]
            QSkampmass = (T["qscal"] * 2.33 * T["rhos"] / (T["rhos"] - T["rhow"]) * m_pow(TP, 1.5)
                          * m_pow(T["tanbeta"], 0.75) * m_pow(T["d50"], -0.25) * m_pow(HSbr, 2)
                          * (m_pow(np.abs(sind(2 * dPHIbr)), 0.6) * np.sign(dPHIbr)
                             - (2 / T["tanbeta"]) * cosd(dPHIbr) * dHS))
            QS = T["qscal"] * 365 * 24 * 60 * 60 * (QSkampmass / T["rhos"]) / (1.0 - T["porosity"])
            QS[np.abs(dPHIbr) > 90] = 0
        elif trform == "MILH":
            HSbr, dPHIbr, TP = WAVE["HSbr"], WAVE["dPHIbr"], WAVE["TP"]
            QSmilhmass = (T["qscal"] * 0.15 * T["rhos"] / (T["rhos"] - T["rhow"]) * m_pow(TP, 0.89)
                          * m_pow(T["tanbeta"], 0.86) * m_pow(T["d50"], -0.69) * m_pow(HSbr, 2.75)
                          * (m_pow(np.abs(sind(2 * dPHIbr)), 0.5) * np.sign(dPHIbr)
                             - (2 / T["tanbeta"]) * cosd(dPHIbr) * dHS))
            QS = T["qscal"] * 365 * 24 * 60 * 60 * (QSmilhmass / T["rhos"]) / (1.0 - T["porosity"])
            QS[np.abs(dPHIbr) > 90] = 0
        else:                                        # VR14
            HSbr, dPHIbr = WAVE["HSbr"], WAVE["dPHIbr"]
            vtide = 0
            kswell = 0.015 * T["pswell"] + (1 - 0.01 * T["pswell"])
            vwave = 0.3 * sind(2 * dPHIbr) * m_pow(T["g"] * HSbr, 0.5)
            vtotal = vwave + vtide
            QSvr14mass = (T["qscal"] * 0.0006 * kswell * T["rhos"] * m_pow(T["tanbeta"], 0.4)
                          * m_pow(T["d50"], -0.6) * m_pow(HSbr, 2.6) * np.abs(vtotal) * np.sign(dPHIbr))
            QS = 365 * 24 * 60 * 60 * (QSvr14mass / T["rhos"]) / (1.0 - T["porosity"])
            QS[np.abs(dPHIbr) > 90] = 0
    QS = np.asarray(QS, dtype=float)
    if field == "QS" and T["suppresshighangle"] == 1:
        idx = WAVE["dPHItdp"] == WAVE["dPHIcrit"]
        QS[idx] = 1.01 * QS[idx]
    return QS


def transport_shadow_treat(COAST, WAVE, TRANSP, find_shadows_mc):
    """Set the transport to zero in the shadow of other parts of the coast."""
    QS = TRANSP["QS"].copy()
    with np.errstate(invalid="ignore"):
        shadowS_90deg = np.abs(WAVE["dPHItdp"]) > 2 * WAVE["dPHIcrit"]
    if TRANSP["trform"].upper() == "CERC":
        WAVE["PHIbr"] = WAVE["PHItdp"]                  # CERC uses the nearshore angle
    shadowS = find_shadows_mc(COAST["xq"], COAST["yq"], COAST["x_mc"], COAST["y_mc"],
                              WAVE["PHIo"], WAVE["PHItdp"], WAVE["PHIbr"], WAVE["dist"])
    shadowS = shadowS | shadowS_90deg
    shadowS_h = np.zeros(QS.size, dtype=bool)
    shadow = np.maximum(shadowS[:-1], shadowS[1:])
    shadowS_hD = np.zeros(COAST["x"].size, dtype=bool)
    dsh = np.diff(shadowS.astype(float))
    idsh = np.flatnonzero((dsh[:-1] - dsh[1:]) == 2)
    shadowS = shadowS.copy()
    shadowS[idsh + 1] = False
    diff = np.asarray(WAVE["diff"], dtype=bool)
    shadowS = shadowS & ~diff
    shadowS_h = shadowS_h & ~diff
    QS[shadowS] = 0
    QS[shadowS_h] = 0
    with np.errstate(invalid="ignore"):
        idpos = np.flatnonzero((QS > 0) & np.concatenate([[False], shadow]))
        QS[idpos] = np.fmax(QS[idpos - 1], 0)
        idneg = np.flatnonzero((QS < 0) & np.concatenate([shadow, [False]]))
        QS[idneg] = np.fmin(QS[idneg + 1], 0)
    TRANSP["QS"] = QS
    TRANSP["shadowS"] = shadowS
    TRANSP["shadowS_h"] = shadowS_h
    TRANSP["shadowS_hD"] = shadowS_hD
    TRANSP["shadow"] = shadowS[:-1] & shadowS[1:]
    TRANSP["shadow_h"] = shadowS_h[:-1] & shadowS_h[1:]
    return TRANSP


def get_smoothangles(COAST, TRANSP):
    """Smooth the transport at sharp coastline angles."""
    COAST["dPHIc"] = mmod(np.diff(COAST["PHIc"]) + 180, 360) - 180
    if COAST["clockwise"]:
        dPHIc = COAST["dPHIc"]
        maxangle = COAST["maxangle"]
        idrev = TRANSP["idrev"]
        ok = (idrev[:-1] == 0) & (idrev[1:] == 0)
        with np.errstate(invalid="ignore"):
            idc1 = np.flatnonzero((dPHIc > maxangle * 0.5) & ok)
            idc2 = np.flatnonzero((dPHIc < -maxangle * 0.5) & ok)
        idc = np.union1d(idc1, idc2)
        idc = np.setdiff1d(idc, np.flatnonzero(TRANSP["shadow"]))
        if not COAST["cyclic"]:
            idc = np.setdiff1d(idc, [0, dPHIc.size - 1])
        QS0 = TRANSP["QS"]
        QS1 = QS0.copy()
        QS2 = QS0.copy()
        if idc.size:
            factor = np.fmin(np.fmax(np.abs(2 * np.abs(dPHIc[idc]) - maxangle) / maxangle, 0), 1)
            ids = dPHIc[idc] * (QS0[idc] - QS0[idc + 1]) < 0
            factor[ids] = 0
            QSfac1 = ((QS0[idc] + QS0[idc + 1]) / 2 <= 0).astype(float)
            QS1[idc] = (1 - factor * QSfac1) * QS0[idc] + factor * QSfac1 * QS0[idc + 1]
            QSfac2 = ((QS0[idc] + QS0[idc + 1]) / 2 >= 0).astype(float)
            QS2[idc + 1] = (1 - factor * QSfac2) * QS0[idc + 1] + factor * QSfac2 * QS0[idc]
            TRANSP["QS"] = (QS1 + QS2) / 2
    return COAST, TRANSP


def _get_mod(x, y):
    z = x % y
    return y if z == 0 else z


def get_upwindcorrection(COAST, WAVE, TRANSP):
    """Upwind correction for high-angle wave instability: where the wave angle
    exceeds the critical angle, the transport is set to the maximum transport
    of the upwind point (sequential sweep in both directions). With
    S.relaxationlength [m], the transport also adapts gradually (over that
    length) to a decrease downdrift (relaxationmethod 2 of ShorelineS)."""
    from .coastline import colon                     # (coastline imports this module's users)
    relax = TRANSP.get("relaxationlength")
    relax = None if relax is None or np.size(relax) == 0 else float(np.asarray(relax).ravel()[0])
    maxangle = max(COAST["maxangle"], 60)
    nq = COAST["nq"]
    dsq = COAST["dsq"]
    QS = TRANSP["QS"].copy()
    QS0 = TRANSP["QS"].copy()
    QS1 = QS0.copy()
    QS2 = QS0.copy()
    QSmax = TRANSP["QSmax"]
    dPHItdp = WAVE["dPHItdp"]
    dPHIcrit = WAVE["dPHIcrit"]
    ivals = np.zeros(nq)
    shadowS = TRANSP["shadowS"]
    shadowS_h = TRANSP["shadowS_h"]
    idrev = TRANSP["idrev"]
    twopoints = TRANSP["twopoints"]
    if COAST["clockwise"] == 1:
        dPHIcr1 = mmod(dPHItdp - dPHIcrit + 180, 360) - 180
        dPHIcr2 = mmod(dPHItdp + dPHIcrit + 180, 360) - 180
        for cw in (1, -1):
            # 1-based indices as in MATLAB
            iirange = list(range(1, nq + 1)) if COAST["cyclic"] else list(range(2, nq))
            dPHIcr = dPHIcr1
            if cw == -1:
                dPHIcr = dPHIcr2
                iirange = iirange[::-1]
            cwd = cw * dPHIcr
            for i in iirange:
                if relax is None and not cwd[i - 1] > 0:
                    continue                         # (nothing to do at this point)
                if COAST["cyclic"]:
                    im1 = _get_mod(i - 1 * cw, nq)
                    ip1 = _get_mod(i + 1 * cw, nq)
                    ip2 = _get_mod(i + 2 * cw, nq)
                elif cw == 1:
                    im1 = max(i - 1, 1)
                    ip1 = min(i + 1, nq)
                    ip2 = min(i + 2, nq)
                else:
                    ip1 = max(i - 1, 1)
                    ip2 = max(i - 2, 1)
                    im1 = min(i + 1, nq)
                if relax is not None:
                    fac0 = dsq[min(i, dsq.size) - 1] / relax
                    fac1 = colon(1.0, -fac0, 0.0)
                    ind = np.minimum(np.maximum(i - cw * np.round(colon(0.0, 1.0, 1 / fac0)), 1), nq).astype(int)
                    if QS0[im1 - 1] > 0 and QS0[i - 1] >= 0 and cw == 1:
                        QS1[i - 1] = _mmax(QS0[i - 1] + np.fmax(QS1[ind - 1] - QS0[i - 1], 0) * fac1)
                    elif QS0[im1 - 1] < 0 and QS0[i - 1] <= 0 and cw == -1:
                        QS2[i - 1] = _mmin(QS0[i - 1] + np.fmin(QS2[ind - 1] - QS0[i - 1], 0) * fac1)
                    if cw == -1:
                        QS[i - 1] = QS1[i - 1] + QS2[i - 1] - QS0[i - 1]
                cond = (cwd[i - 1] > 0
                        and (cwd[im1 - 1] <= 0 or (cw * QS[im1 - 1] > 0 and (cwd[im1 - 1] - cwd[i - 1]) < -maxangle))
                        and not shadowS[im1 - 1]
                        and (shadowS_h.size == 0 or (not shadowS_h[ip1 - 1] and not shadowS_h[ip2 - 1]))
                        and idrev[i - 1] == 0 and idrev[min(max(i - cw, 1), nq) - 1] == 0)
                if cond:
                    vals = 1.05 * QSmax[[im1 - 1, i - 1, ip1 - 1]]
                    QS[i - 1] = cw * (np.nan if np.all(np.isnan(vals)) else np.nanmax(vals))
                    cf = 1.1
                    if twopoints == 1:
                        QS[ip1 - 1] = QS[ip1 - 1] + cf * cw * max(0.0, (cw * QS[i - 1] - cw * QS[ip1 - 1]) / 2)
                    elif twopoints == 2:
                        QS[ip1 - 1] = QS[ip1 - 1] + cf * cw * max(0.0, (cw * QS[i - 1] - cw * QS[ip1 - 1]) / 2)
                        QS[ip2 - 1] = QS[ip2 - 1] + cf * cw * max(0.0, (cw * QS[ip1 - 1] - cw * QS[ip2 - 1]) / 2)
                    ivals[i - 1] = cw
                    QSk = QS1 if cw == 1 else QS2
                    QSk[i - 1] = QS[i - 1]
                    QSk[ip1 - 1] = QS[ip1 - 1]
                    QSk[ip2 - 1] = QS[ip2 - 1]
    TRANSP["QS"] = QS
    TRANSP["ivals"] = ivals
    return TRANSP


def _mmax(a):
    """MATLAB max of a vector (NaNs ignored unless all are NaN)."""
    a = np.asarray(a, dtype=float)
    return np.nan if a.size == 0 or np.all(np.isnan(a)) else float(np.nanmax(a))


def _mmin(a):
    a = np.asarray(a, dtype=float)
    return np.nan if a.size == 0 or np.all(np.isnan(a)) else float(np.nanmin(a))


def transport_boundary_condition(TRANSP, COAST):
    QS0 = TRANSP["QS"].copy()
    QS = TRANSP["QS"]
    if not COAST["cyclic"]:
        # 'Angleconstant'/'Gradient' act on the coastline angle (get_coastline_orientation)
        bs, vs = TRANSP["boundaryconditionstart"]
        be, ve = TRANSP["boundaryconditionend"]
        bs, be = bs.lower(), be.lower()
        if bs == "periodic":
            QS[0] = (QS0[-2] + QS0[1]) / 2
        elif bs == "closed":
            QS[0] = 0 if np.isnan(vs) else vs      # {'Closed', Q}: fixed transport Q [m3/yr]
        elif bs in ("neumann", "fixed"):
            QS[0] = QS[1]
        if be == "periodic":
            QS[-1] = (QS0[-2] + QS0[1]) / 2
        elif be == "closed":
            QS[-1] = 0 if np.isnan(ve) else ve
        elif be in ("neumann", "fixed"):
            QS[-1] = QS[-2]
    TRANSP["QS"] = QS
    return TRANSP
