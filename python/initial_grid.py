"""Python port of initial_grid.m.

Reconstructs CoastSat shoreline positions from the PCA model and restores the
transects that were removed from the PCA (held at their first valid CoastSat
observation). Beaches are joined with one NaN between them.

    x, y, x_initial, y_initial, g = initial_grid([2, 3], dy, n_d)

x, y            time x transects (full daily time series, NaN separators)
x_initial/y_... first row = initial coastline (NZTM2000 / EPSG:2193)
g               [x_initial, y_initial] as an (n x 2) array
"""
import os

import numpy as np
import pandas as pd
import scipy.io
from pyproj import Transformer

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

FILE_NAMES = ["nzd0204.xlsx", "nzd0207.xlsx", "nzd0217.xlsx", "nzd0222.xlsx", "nzd0220.xlsx",
              "nzd0226.xlsx", "nzd0229.xlsx", "nzd0231.xlsx", "nzd0234.xlsx", "nzd0239.xlsx",
              "nzd0240.xlsx", "nzd0238.xlsx", "nzd0236.xlsx", "nzd0233.xlsx", "nzd0230.xlsx",
              "nzd0227.xlsx"]


def _r(a, b):
    """MATLAB a:b (1-based, inclusive)."""
    return list(range(a, b + 1))


# Transects removed from the PCA reconstruction (1-based, per beach 1..16)
REMOVE = {
    1: _r(92, 112), 2: _r(1, 60) + _r(233, 293), 3: _r(1, 38), 4: _r(36, 43), 5: _r(14, 16),
    6: _r(165, 166), 7: _r(190, 196), 8: _r(1, 5) + _r(98, 106), 9: _r(145, 150),
    10: _r(110, 141), 11: _r(1, 136) + _r(207, 213) + _r(254, 264), 12: [], 13: _r(46, 49),
    14: [], 15: _r(26, 60), 16: _r(27, 43),
}

# NZGD2000 geographic -> NZTM2000 (as MATLAB's projfwd(projcrs(2193), lat, lon))
_TO_NZTM = Transformer.from_crs("EPSG:4167", "EPSG:2193", always_xy=True)


def projfwd(lat, lon):
    x, y = _TO_NZTM.transform(np.asarray(lon, dtype=float), np.asarray(lat, dtype=float))
    return np.asarray(x), np.asarray(y)


def _parse_latlon(cells):
    """Parse 'lat,lon' / 'lat lon' text cells; anything else gives NaN."""
    s = pd.Series(np.asarray(cells, dtype=object).ravel())
    txt = s.where(s.map(lambda v: isinstance(v, str)), None)
    parts = txt.str.replace(" ", ",", regex=False).str.split(",", n=2, expand=True)
    lat = np.full(s.size, np.nan)
    lon = np.full(s.size, np.nan)
    if parts is not None and parts.shape[1] >= 2:
        ok = parts[1].notna()
        lat[ok.values] = pd.to_numeric(parts[0][ok], errors="coerce").values
        lon[ok.values] = pd.to_numeric(parts[1][ok], errors="coerce").values
    return lat.reshape(np.shape(cells)), lon.reshape(np.shape(cells))


def initial_grid(beach, dy=None, n_d=None, data_dir=REPO, verbose=True):
    A = scipy.io.loadmat(os.path.join(data_dir, "A_site2.mat"))["A"]
    mint = scipy.io.loadmat(os.path.join(data_dir, "mint_site2.mat"))["mint"]
    x = y = None
    x_initial = y_initial = None
    for b in beach:
        if verbose:
            print(f"\nProcessing beach {b}...")
        fname = os.path.join(data_dir, FILE_NAMES[b - 1])
        T = pd.read_excel(fname, sheet_name="Transects")
        land_x = T.iloc[:, 17].to_numpy(float)          # longitude
        land_y = T.iloc[:, 18].to_numpy(float)          # latitude
        phi = T.iloc[:, 2].to_numpy(float)              # transect orientation
        nT = land_x.size
        idx_remove = np.asarray(REMOVE[b], dtype=int) - 1
        valid = np.ones(nT, dtype=bool)
        valid[idx_remove] = False

        land_x_utm, land_y_utm = projfwd(land_y, land_x)
        theta = np.deg2rad(phi)
        ux, uy = np.sin(theta), np.cos(theta)

        A_abs = A[0, b - 1] + mint[0, b - 1]            # time x retained transects
        if A_abs.shape[1] != valid.sum():
            raise ValueError(f"Beach {b}: number of PCA transects ({A_abs.shape[1]}) does not match "
                             f"number of valid transects ({valid.sum()}).")
        x_valid = land_x_utm[valid][None, :] + A_abs * ux[valid][None, :]
        y_valid = land_y_utm[valid][None, :] + A_abs * uy[valid][None, :]
        nTime = x_valid.shape[0]

        # Satellite observations ('Intersect points', one 'lat,lon' text per transect).
        # MATLAB's readtable(..., 'Range','B2:KH517') takes the first row of the
        # range (spreadsheet row 2, the first image date) as variable names, so the
        # observations start at row 3. This is reproduced here (it matches the
        # initial coastline stored in 30_sept_int_w/output.mat to 1e-8 m).
        obs = pd.read_excel(fname, sheet_name="Intersect points", header=None, skiprows=2,
                            usecols=range(1, 294), nrows=515, dtype=object)   # columns B:KH
        cols = [c for c in idx_remove if c < obs.shape[1]]
        x_obs = np.full(obs.shape, np.nan)
        y_obs = np.full(obs.shape, np.nan)
        if cols:
            lat, lon = _parse_latlon(obs.iloc[:, cols].to_numpy())
            x_obs[:, cols], y_obs[:, cols] = projfwd(lat, lon)

        x_b = np.full((nTime, nT), np.nan)
        y_b = np.full((nTime, nT), np.nan)
        x_b[:, valid] = x_valid
        y_b[:, valid] = y_valid
        for tr in idx_remove:
            if tr >= x_obs.shape[1]:
                print(f"Warning: beach {b}, transect {tr + 1}: no corresponding satellite column found.")
                continue
            ok = ~np.isnan(x_obs[:, tr]) & ~np.isnan(y_obs[:, tr])
            if ok.any():
                first = np.flatnonzero(ok)[0]
                x_b[:, tr] = x_obs[first, tr]
                y_b[:, tr] = y_obs[first, tr]
            else:
                print(f"Warning: beach {b}, transect {tr + 1}: no valid satellite observation found.")

        if x is None:
            x, y = x_b, y_b
            x_initial, y_initial = x_b[0], y_b[0]
        else:
            x = np.hstack([x, np.full((x.shape[0], 1), np.nan), x_b])
            y = np.hstack([y, np.full((y.shape[0], 1), np.nan), y_b])
            x_initial = np.concatenate([x_initial, [np.nan], x_b[0]])
            y_initial = np.concatenate([y_initial, [np.nan], y_b[0]])
        if verbose:
            print(f"  Total transects: {nT}\n  PCA transects:   {valid.sum()}\n  Removed:         {idx_remove.size}")
    g = np.stack([x_initial, y_initial], axis=1)
    if verbose:
        print(f"\nINITIAL GRID COMPLETE: x {x.shape}, separators at {np.flatnonzero(np.isnan(x_initial)) + 1}")
    return x, y, x_initial, y_initial, g


if __name__ == "__main__":
    initial_grid([2, 3], 153, 9130)
