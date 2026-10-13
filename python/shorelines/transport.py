"""Longshore sediment transport and its corrections. Ports of transport.m
(Kamphuis formula), transport_shadow_treat.m, get_smoothangles.m,
get_upwindcorrection.m and transport_boundary_condition.m."""
import numpy as np

from .mfuncs import cosd, m_pow, mmod, sind


def transport(TRANSP, WAVE):
    """Bulk longshore transport [m3/yr] at the transport points (Kamphuis 1991)."""
    if TRANSP["trform"].upper() != "KAMP":
        raise NotImplementedError(f"transport formula {TRANSP['trform']} (only KAMP is ported)")
    HSbr = WAVE["HSbr"]
    dPHIbr = WAVE["dPHIbr"]
    TP = WAVE["TP"]
    dHS = np.zeros(HSbr.size)
    T = TRANSP
    with np.errstate(invalid="ignore"):
        QSkampmass = (T["qscal"] * 2.33 * T["rhos"] / (T["rhos"] - T["rhow"]) * m_pow(TP, 1.5)
                      * m_pow(T["tanbeta"], 0.75) * m_pow(T["d50"], -0.25) * m_pow(HSbr, 2)
                      * (m_pow(np.abs(sind(2 * dPHIbr)), 0.6) * np.sign(dPHIbr)
                         - (2 / T["tanbeta"]) * cosd(dPHIbr) * dHS))
        QS = T["qscal"] * 365 * 24 * 60 * 60 * (QSkampmass / T["rhos"]) / (1.0 - T["porosity"])
        QS[np.abs(dPHIbr) > 90] = 0
    return QS


def transport_shadow_treat(COAST, WAVE, TRANSP, find_shadows_mc):
    """Set the transport to zero in the shadow of other parts of the coast."""
    QS = TRANSP["QS"].copy()
    with np.errstate(invalid="ignore"):
        shadowS_90deg = np.abs(WAVE["dPHItdp"]) > 2 * WAVE["dPHIcrit"]
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
    of the upwind point (sequential sweep in both directions)."""
    if TRANSP.get("relaxationlength") not in (None, []) and np.size(TRANSP["relaxationlength"]):
        raise NotImplementedError("S.relaxationlength")
    maxangle = max(COAST["maxangle"], 60)
    nq = COAST["nq"]
    QS = TRANSP["QS"].copy()
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
                if not cwd[i - 1] > 0:
                    continue
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
                cond = ((cwd[im1 - 1] <= 0 or (cw * QS[im1 - 1] > 0 and (cwd[im1 - 1] - cwd[i - 1]) < -maxangle))
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
    TRANSP["QS"] = QS
    TRANSP["ivals"] = ivals
    return TRANSP


def transport_boundary_condition(TRANSP, COAST):
    QS0 = TRANSP["QS"].copy()
    QS = TRANSP["QS"]
    if not COAST["cyclic"]:
        bs = TRANSP["boundaryconditionstart"].lower()
        be = TRANSP["boundaryconditionend"].lower()
        if bs == "periodic":
            QS[0] = (QS0[-2] + QS0[1]) / 2
        elif bs == "closed":
            QS[0] = 0
        elif bs in ("neumann", "fixed"):
            QS[0] = QS[1]
        if be == "periodic":
            QS[-1] = (QS0[-2] + QS0[1]) / 2
        elif be == "closed":
            QS[-1] = 0
        elif be in ("neumann", "fixed"):
            QS[-1] = QS[-2]
    TRANSP["QS"] = QS
    return TRANSP
