"""Python port of project_output.m.

Projects the modelled and the observed (CoastSat/PCA) shorelines onto the
initial coastline (normals of +-Lcrit) and plots the cross-shore change of the
model, of the daily observations and of their 30-day moving mean.
"""
import os

import numpy as np
import scipy.io

from shorelines.geometry import get_polydistance

HERE = os.path.dirname(os.path.abspath(__file__))


def movmean_omitnan(a, window):
    """MATLAB movmean(a, window, 2, 'omitnan') along the rows (centred window)."""
    a = np.asarray(a, dtype=float)
    n = a.shape[1]
    before = window // 2
    after = window - 1 - before if window % 2 == 0 else window // 2
    if window % 2 == 1:
        before = after = (window - 1) // 2
    ok = ~np.isnan(a)
    vals = np.where(ok, a, 0.0)
    cs = np.concatenate([np.zeros((a.shape[0], 1)), np.cumsum(vals, axis=1)], axis=1)
    cn = np.concatenate([np.zeros((a.shape[0], 1)), np.cumsum(ok, axis=1)], axis=1)
    i = np.arange(n)
    lo = np.clip(i - before, 0, n)
    hi = np.clip(i + after + 1, 0, n)
    s = cs[:, hi] - cs[:, lo]
    c = cn[:, hi] - cn[:, lo]
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(c > 0, s / c, np.nan)


def project_output(x, y, file_out, dy, interv, x_initial, y_initial, n_d, O=None, plot=True,
                   Lcrit=500.0, figfile=None):
    """Returns dict with zg_model, zg_obs (every `interv` days), zg_obs_d (daily) and sm."""
    if O is None:
        raw = scipy.io.loadmat(os.path.join(file_out, "output.mat"), squeeze_me=True, struct_as_record=False)
        O = raw["O"]
        x_model, y_model = np.atleast_2d(O.x), np.atleast_2d(O.y)
    else:
        x_model, y_model = O["x"], O["y"]
    if x_model.shape[0] == 1:
        x_model, y_model = x_model.T, y_model.T

    # observed data (MATLAB rows dy:dy+n_d, 1-based)
    x_obs = x[dy - 1:dy + n_d, :]
    y_obs = y[dy - 1:dy + n_d, :]
    if interv > 1:
        idx = np.arange(0, x_obs.shape[0], interv)
        x_obst, y_obst = x_obs[idx].T, y_obs[idx].T
    else:
        x_obst, y_obst = x_obs.T, y_obs.T
    nGrid = x_obs.shape[1]
    ntime = min(x_model.shape[1], x_obst.shape[1])
    zg_model = np.full((nGrid, ntime), np.nan)
    zg_obs = np.full((nGrid, ntime), np.nan)
    for t in range(ntime):
        zg_model[:, t] = get_polydistance(x_initial, y_initial, x_model[:, t], y_model[:, t], Lcrit)[0]
        zg_obs[:, t] = get_polydistance(x_initial, y_initial, x_obst[:, t], y_obst[:, t], Lcrit)[0]
    # daily observed distance (independent of interv)
    x_obst_d, y_obst_d = x_obs.T, y_obs.T
    ntime_d = x_obst_d.shape[1]
    zg_obs_d = np.full((nGrid, ntime_d), np.nan)
    for t in range(ntime_d):
        zg_obs_d[:, t] = get_polydistance(x_initial, y_initial, x_obst_d[:, t], y_obst_d[:, t], Lcrit)[0]
    sm = movmean_omitnan(zg_obs_d, 30)
    if plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
        for ax, data, title in zip(axs, (zg_model, zg_obs_d, sm),
                                   ("Model", "CoastSat", "CoastSat (30-day moving mean)")):
            im = ax.imshow(data, aspect="auto", cmap="viridis", interpolation="nearest")
            fig.colorbar(im, ax=ax)
            ax.set_title(title)
            ax.set_xlabel("Time step")
        axs[0].set_ylabel("Alongshore position")
        fig.savefig(figfile or os.path.join(file_out, "project_output.png"), dpi=120)
        plt.close(fig)
    return dict(zg_model=zg_model, zg_obs=zg_obs, zg_obs_d=zg_obs_d, sm=sm)
