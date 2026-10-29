
"""Python port of initial_grid.m.

Reconstructs CoastSat shoreline positions from the PCA model, excluding
transects removed from the PCA reconstruction. Beaches are joined with
one NaN separator between them.

    x, y, x_initial, y_initial, g = initial_grid([2, 3], dy, n_d)

x, y         : time x retained transects, with NaN separators between beaches
x_initial    : initial x-coordinate for each retained transect
y_initial    : initial y-coordinate for each retained transect
g            : [x_initial, y_initial], shape (n_retained, 2)

Coordinates are in NZTM2000 / EPSG:2193.
Time starts on 1 August 1999.
"""

import os

import numpy as np
import pandas as pd
import scipy.io
from pyproj import Transformer

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

FILE_NAMES = [
    "nzd0204.xlsx", "nzd0207.xlsx", "nzd0217.xlsx", "nzd0222.xlsx",
    "nzd0220.xlsx", "nzd0226.xlsx", "nzd0229.xlsx", "nzd0231.xlsx",
    "nzd0234.xlsx", "nzd0239.xlsx", "nzd0240.xlsx", "nzd0238.xlsx",
    "nzd0236.xlsx", "nzd0233.xlsx", "nzd0230.xlsx", "nzd0227.xlsx"
]


def _r(a, b):
    """MATLAB a:b (1-based, inclusive)."""
    return list(range(a, b + 1))


# Transects excluded from PCA, numbered from 1 for each beach
REMOVE = {
    1: _r(92, 112),
    2: _r(1, 60) + _r(233, 293),
    3: _r(1, 38),
    4: _r(36, 43),
    5: _r(14, 16),
    6: _r(165, 166),
    7: _r(190, 196),
    8: _r(1, 5) + _r(98, 106),
    9: _r(145, 150),
    10: _r(110, 141),
    11: _r(1, 136) + _r(207, 213) + _r(254, 264),
    12: [],
    13: _r(46, 49),
    14: [],
    15: _r(26, 60),
    16: _r(27, 43),
}

# NZGD2000 geographic -> NZTM2000
_TO_NZTM = Transformer.from_crs(
    "EPSG:4167", "EPSG:2193", always_xy=True
)


def projfwd(lat, lon):
    """Convert latitude/longitude to NZTM2000 coordinates."""
    x, y = _TO_NZTM.transform(
        np.asarray(lon, dtype=float),
        np.asarray(lat, dtype=float)
    )
    return np.asarray(x), np.asarray(y)


def initial_grid(beach, dy=None, n_d=None, data_dir=REPO, verbose=True):

    # PCA reconstruction and time-mean coordinates
    A = scipy.io.loadmat(
        os.path.join(data_dir, "PC_lst_site2.mat")
    )["PC_lst"]

    mint = scipy.io.loadmat(
        os.path.join(data_dir, "mean_PC_lst_site2.mat")
    )["mean_PC_lst"]

    x = y = None
    x_initial = y_initial = None

    for b in beach:

        if verbose:
            print(f"\nProcessing beach {b}...")

        fname = os.path.join(data_dir, FILE_NAMES[b - 1])
        T = pd.read_excel(fname, sheet_name="Transects")

        # Reference coordinates: longitude and latitude
        land_x = T.iloc[:, 17].to_numpy(dtype=float)
        land_y = T.iloc[:, 18].to_numpy(dtype=float)

        # Transect orientation
        phi = T.iloc[:, 2].to_numpy(dtype=float)

        nT = land_x.size

        # Identify transects retained in the PCA reconstruction
        idx_remove = np.asarray(REMOVE[b], dtype=int) - 1

        valid = np.ones(nT, dtype=bool)
        valid[idx_remove] = False

        # Convert reference coordinates to NZTM2000
        land_x_utm, land_y_utm = projfwd(land_y, land_x)

        theta = np.deg2rad(phi)
        ux = np.sin(theta)
        uy = np.cos(theta)

        # PCA reconstruction: time x retained transects
        A_abs = A[0, b - 1] + mint[0, b - 1]

        if A_abs.shape[1] != valid.sum():
            raise ValueError(
                f"Beach {b}: PCA contains {A_abs.shape[1]} transects, "
                f"but {valid.sum()} transects remain after exclusion."
            )

        # Reconstruct only retained transects
        x_b = (
            land_x_utm[valid][None, :]
            + A_abs * ux[valid][None, :]
        )

        y_b = (
            land_y_utm[valid][None, :]
            + A_abs * uy[valid][None, :]
        )

        # Initial coastline at the requested time index
        i0 = int(dy or 0)

        if i0 < 0 or i0 >= x_b.shape[0]:
            raise IndexError(
                f"dy={i0} is outside the available time range "
                f"for beach {b}."
            )

        x_b_initial = x_b[i0, :]
        y_b_initial = y_b[i0, :]

        # Combine beaches, with one NaN separator between them
        if x is None:
            x, y = x_b, y_b
            x_initial = x_b_initial.copy()
            y_initial = y_b_initial.copy()

        else:
            x = np.hstack([
                x,
                np.full((x.shape[0], 1), np.nan),
                x_b
            ])

            y = np.hstack([
                y,
                np.full((y.shape[0], 1), np.nan),
                y_b
            ])

            x_initial = np.concatenate([
                x_initial, [np.nan], x_b_initial
            ])

            y_initial = np.concatenate([
                y_initial, [np.nan], y_b_initial
            ])

        if verbose:
            print(f"  Original transects: {nT}")
            print(f"  Retained in PCA:    {valid.sum()}")
            print(f"  Excluded:           {idx_remove.size}")

    # Initial coastline coordinates: n transects x 2
    g = np.column_stack([x_initial, y_initial])

    if verbose:
        print(f"\nINITIAL GRID COMPLETE")
        print(f"x shape: {x.shape}")
        print(f"y shape: {y.shape}")
        print(f"g shape: {g.shape}")
        print(
            "NaN separators at transect positions:",
            np.flatnonzero(np.isnan(x_initial)) + 1
        )

    return x, y, x_initial, y_initial, g


if __name__ == "__main__":
    x, y, x_initial, y_initial, g = initial_grid(
        [2, 3], dy=0, n_d=9283
    )
import scipy.io

# Run initial_grid
x, y, x_initial, y_initial, g = initial_grid(
    [2, 3], dy=0, n_d=9283
)

# Save g and the initial coordinates to a MAT file
scipy.io.savemat(
    os.path.join(REPO, "g_initial.mat"),
    {
        "g": g,
        "x_initial": x_initial,
        "y_initial": y_initial
    }
)

print("Saved:", os.path.join(REPO, "g_initial.mat"))