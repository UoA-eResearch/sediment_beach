"""Why the model and the CoastSat/PCA shorelines change differently over time.

    python compare_model_coastsat.py

Prints the statistics quoted in the README section "Why ShorelineS does not
reproduce the CoastSat pattern". Uses the committed results
(results/python_hindcast_2000_2024/), the wave input (input_25_contour.nc)
and the CoastSat workbooks (Intersects and Tides sheets).
"""
import os
import warnings

import netCDF4
import numpy as np
import pandas as pd
import scipy.io

warnings.filterwarnings("ignore", category=RuntimeWarning)
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
RES = os.path.join(REPO, "results", "python_hindcast_2000_2024")

BEACHES = {"beach 2 (nzd0207)": (np.arange(0, 293), "nzd0207.xlsx"),
           "beach 3 (nzd0217)": (np.arange(294, 609), "nzd0217.xlsx")}
NOT_PCA = [(0, 60), (232, 293), (294, 332)]


def pca_mask(n):
    m = np.ones(n, bool)
    m[293] = False
    for a, b in NOT_PCA:
        m[a:b] = False
    return m


def detrend_rows(z):
    t = np.arange(z.shape[1])
    out = np.full_like(z, np.nan)
    for i in range(z.shape[0]):
        ok = ~np.isnan(z[i])
        if ok.sum() > 10:
            out[i, ok] = z[i, ok] - np.polyval(np.polyfit(t[ok], z[i, ok], 1), t[ok])
    return out


def highpass(a, k=24):
    """Remove a k-sample running mean (k=24 x 30 days = 2 years)."""
    out = np.full_like(a, np.nan)
    ker = np.ones(k) / k
    for i in range(a.shape[0]):
        if np.isfinite(a[i]).sum() > k:
            out[i] = a[i] - np.convolve(np.nan_to_num(a[i]), ker, "same")
    return out[:, k // 2:-k // 2]


def main():
    P = scipy.io.loadmat(os.path.join(RES, "projection.mat"))
    zm, zo = P["zg_model"], P["zg_obs"]          # every 30 days, transects x time
    pca = pca_mask(zm.shape[0])

    print("1. Share of the variance that is the whole beach moving together (PCA transects)")
    for name, (ix, _) in BEACHES.items():
        ix = ix[pca[ix]]
        z = zo[ix] - np.nanmean(zo[ix], axis=1, keepdims=True)
        frac = np.nanvar(np.broadcast_to(np.nanmean(z, axis=0), z.shape)) / np.nanvar(z)
        print(f"   {name}: CoastSat {frac * 100:.0f}%")

    print("2. Beach-average cross-shore position (sediment conservation)")
    for name, (ix, _) in BEACHES.items():
        m_model = np.nanmean(zm[ix], axis=0)
        m_obs = np.nanmean(zo[ix[pca[ix]]], axis=0)
        print(f"   {name}: model (all transects) std {np.nanstd(m_model):.2f} m | "
              f"CoastSat (PCA transects) std {np.nanstd(m_obs):.1f} m, range {np.nanmin(m_obs):.0f}..{np.nanmax(m_obs):.0f} m")

    print("3. Variability on time scales < ~2 years (detrended, 2-year running mean removed)")
    for name, (ix, _) in BEACHES.items():
        ix = ix[pca[ix]]
        hm, ho = highpass(detrend_rows(zm[ix])), highpass(detrend_rows(zo[ix]))
        r = [np.corrcoef(hm[i], ho[i])[0, 1] for i in range(hm.shape[0])
             if np.all(np.isfinite(hm[i])) and np.all(np.isfinite(ho[i]))]
        rm = np.corrcoef(np.nanmean(hm, axis=0), np.nanmean(ho, axis=0))[0, 1]
        print(f"   {name}: std model {np.nanmedian(np.nanstd(hm, axis=1)):.2f} m, CoastSat "
              f"{np.nanmedian(np.nanstd(ho, axis=1)):.2f} m; correlation per transect {np.median(r):.2f}, "
              f"beach-average {rm:.2f}")

    print("4. Seasonal cycle of the raw CoastSat positions (Intersects sheet), with and without the")
    print("   part that correlates with the tide at image time")
    for name, (_, f) in BEACHES.items():
        X = pd.read_excel(os.path.join(REPO, f), sheet_name="Intersects")
        T = pd.read_excel(os.path.join(REPO, f), sheet_name="Tides")
        X["dates"] = pd.to_datetime(X["dates"], utc=True)
        T["dates"] = pd.to_datetime(T["dates"], utc=True)
        D = pd.merge_asof(X.sort_values("dates"), T.sort_values("dates"), on="dates",
                          direction="nearest", tolerance=pd.Timedelta("1h"))
        V = D.drop(columns=["dates", "tide"]).apply(pd.to_numeric, errors="coerce")
        V = V.loc[:, V.notna().mean() > 0.5]
        bm = (V - V.mean()).mean(axis=1)
        t = (D["dates"] - D["dates"].iloc[0]).dt.days.values.astype(float)
        ok = (bm.notna() & D["tide"].notna()).values
        y = bm[ok].values - np.polyval(np.polyfit(t[ok], bm[ok].values, 1), t[ok])
        tide = D["tide"][ok].values
        slope = np.polyfit(tide, y, 1)[0]
        mon = D["dates"][ok].dt.month.values
        clim = np.array([y[mon == m].mean() for m in range(1, 13)])
        yres = y - slope * (tide - tide.mean())
        clim2 = np.array([yres[mon == m].mean() for m in range(1, 13)])
        tcl = np.array([tide[mon == m].mean() for m in range(1, 13)])
        print(f"   {name}: corr(position, tide) {np.corrcoef(tide, y)[0, 1]:+.2f}")
        print("      position anomaly J..D (m):", " ".join(f"{v:+.1f}" for v in clim),
              f"| peak-to-peak {np.ptp(clim):.1f} m")
        print("      tide effect removed  (m):", " ".join(f"{v:+.1f}" for v in clim2),
              f"| peak-to-peak {np.ptp(clim2):.1f} m")
        print("      tide at image time   (m):", " ".join(f"{v:+.2f}" for v in tcl))

    print("5. Wave climate at the 25 m contour by month")
    with netCDF4.Dataset(os.path.join(REPO, "input_25_contour.nc")) as d:
        t = d["time"][:] / 86400.0
        hs = d["point_hm0"][:].filled(np.nan)
        tp = d["point_tp"][:].filled(np.nan)
        wd = d["point_wavdir"][:].filled(np.nan)
    mon = (np.datetime64("1999-08-01") + t.astype("timedelta64[D]")).astype("datetime64[M]").astype(int) % 12 + 1
    hmean = np.nanmean(hs, axis=0)
    hmax = np.nanmax(hs, axis=0)
    s = np.nanmean(np.sin(np.deg2rad(wd)), axis=0)
    c = np.nanmean(np.cos(np.deg2rad(wd)), axis=0)
    print("   month  mean Hs  Tp    dir   % time Hs>2 m")
    for m in range(1, 13):
        k = mon == m
        direc = np.rad2deg(np.arctan2(s[k].mean(), c[k].mean())) % 360
        print(f"   {m:5d}  {hmean[k].mean():6.2f} {np.nanmean(tp[:, k]):5.1f} {direc:5.0f}   {np.mean(hmax[k] > 2) * 100:5.1f}")


if __name__ == "__main__":
    main()
