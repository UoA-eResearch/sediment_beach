"""Figures for the README (written to docs/figures/).

    python make_figures.py

Inputs (committed, Git LFS):
  results/python_hindcast_2000_2024/output.mat      full Python hindcast (O structure)
  results/python_hindcast_2000_2024/projection.mat  project_output results
  results/benchmarks.csv                            measured / estimated runtimes
  30_sept_int_w/output.mat                          the MATLAB run stored in the repo
"""
import csv
import datetime as dt
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import scipy.io  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

from shorelines.geometry import get_polydistance  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
RES = os.path.join(REPO, "results", "python_hindcast_2000_2024")
FIG = os.path.join(REPO, "docs", "figures")

# --- reference palette (light mode) ---------------------------------------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
NEUTRAL = "#f0efec"
BLUE, ORANGE = "#2a78d6", "#eb6834"
BLUE_LIGHT = "#9ec5f4"
# diverging: blue (seaward / accretion) <-> gray <-> red (landward / erosion)
DIVERGING = LinearSegmentedColormap.from_list(
    "accretion_erosion",
    ["#8a1f1f", "#c0392f", "#e34948", "#f4b9b4", NEUTRAL, "#b7d3f6", "#3987e5", "#1c5cab", "#0d366b"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.family": ["DejaVu Sans"], "font.size": 10,
    "text.color": INK, "axes.labelcolor": INK2, "axes.titlecolor": INK, "axes.titlesize": 11,
    "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "lines.linewidth": 1.6,
})

BEACH2 = slice(0, 293)          # transects of beach 2 (nzd0207)
BEACH3 = slice(294, 609)        # transects of beach 3 (nzd0217); 293 is the NaN separator
NOT_PCA = [(0, 60), (232, 293), (294, 332)]   # transects held fixed (removed from the PCA)


def datenum_to_date(d):
    return dt.date.fromordinal(int(d) - 366)


def load():
    O = scipy.io.loadmat(os.path.join(RES, "output.mat"), squeeze_me=True, struct_as_record=False)["O"]
    P = scipy.io.loadmat(os.path.join(RES, "projection.mat"))
    return O, P


def shade_not_pca(ax, orient="x"):
    for a, b in NOT_PCA:
        if orient == "x":
            ax.axvspan(a - 0.5, b - 0.5, color=NEUTRAL, zorder=0, lw=0)
        else:
            ax.axhspan(a - 0.5, b - 0.5, color=NEUTRAL, zorder=0, lw=0)


# ---------------------------------------------------------------------------
def fig_runtime():
    rows = list(csv.DictReader(open(os.path.join(REPO, "results", "benchmarks.csv"))))

    def get(task, impl):
        r = [r for r in rows if r["task"] == task and r["implementation"] == impl][0]
        return float(r["seconds"]), r["kind"]

    impls = ["Octave, original MATLAB code", "Octave, optimised MATLAB code", "Python port"]
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.3), constrained_layout=True)
    for ax, task, unit, scale, title in (
            (axs[0], "model_full_hindcast", "hours", 3600, "Full 2000–2024 hindcast (73,041 time steps)"),
            (axs[1], "project_output_full", "seconds", 1, "project_output, 9,131 daily coastlines")):
        vals, kinds = zip(*[get(task, i) for i in impls])
        vals = np.array(vals) / scale
        y = np.arange(len(impls))[::-1]
        for yi, v, k in zip(y, vals, kinds):
            ax.barh(yi, v, height=0.56, color=BLUE if k == "measured" else BLUE_LIGHT, zorder=2)
            txt = (f"{v:,.1f} h" if unit == "hours" and v >= 1 else
                   f"{v * 60:,.0f} min" if unit == "hours" else
                   f"{v / 3600:,.1f} h" if v >= 3600 else f"{v:,.1f} s")
            ax.text(v * 1.12, yi, txt + ("  (estimated)" if k != "measured" else ""), va="center",
                    color=INK, fontsize=9.5)
        ax.set_xscale("log")
        ax.set_yticks(y, impls)
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="y", visible=False)
        ax.set_xlim(vals.min() / 3, vals.max() * 30)
        ax.set_xlabel(f"wall-clock time, {unit} (log scale)")
        ax.set_title(title)
    axs[1].set_yticklabels([])
    fig.suptitle("Runtime, same machine, one run at a time (GNU Octave 8.4 vs Python 3.11)",
                 x=0.01, ha="left", color=INK2, fontsize=9.5)
    fig.savefig(os.path.join(FIG, "runtime.png"), dpi=150)
    plt.close(fig)


def fig_validation():
    M = scipy.io.loadmat(os.path.join(REPO, "30_sept_int_w", "output.mat"), squeeze_me=True,
                         struct_as_record=False)
    OM, S = M["O"], M["S"]
    O, _ = load()
    tm, tp = np.atleast_1d(OM.timenum), np.atleast_1d(O.timenum)
    dates, dmax, dmed = [], [], []
    for k, t in enumerate(tm):
        j = np.flatnonzero(np.abs(tp - t) < 1e-6)
        if not j.size:
            continue
        j = j[0]
        nm = np.max(np.flatnonzero(~np.isnan(OM.x[:, k]))) + 1
        npy = np.max(np.flatnonzero(~np.isnan(O.x[:, j]))) + 1
        assert nm == npy, "different number of grid points"
        dM = get_polydistance(S.xmc, S.ymc, OM.x[:nm, k], OM.y[:nm, k], 500)[0]
        dP = get_polydistance(S.xmc, S.ymc, O.x[:npy, j], O.y[:npy, j], 500)[0]
        d = np.abs(dM - dP)
        dates.append(datenum_to_date(t))
        dmax.append(np.nanmax(d))
        dmed.append(np.nanmedian(d))
    fig, ax = plt.subplots(figsize=(9, 3.6), constrained_layout=True)
    ax.plot(dates, dmax, color=BLUE, label="largest difference (any transect)")
    ax.plot(dates, dmed, color=ORANGE, label="median difference")
    ax.set_yscale("log")
    ax.set_ylim(1e-11, 1e-3)
    ax.axhline(1e-3, color=AXIS, lw=0.8)
    ax.text(dates[-1], dmax[-1] * 1.8, "largest", color=INK2, ha="right", fontsize=9)
    ax.text(dates[-1], dmed[-1] / 3.2, "median", color=INK2, ha="right", fontsize=9)
    ax.set_ylabel("|Python − MATLAB| cross-shore position (m)")
    ax.legend(loc="upper left", ncols=2, fontsize=9)
    ax.set_title("Python port vs. the MATLAB run stored in 30_sept_int_w (2000-01 to 2004-03)")
    ax.text(0.995, 0.03, "1 mm would be at the top of the axis; differences stay below 1e-7 m",
            transform=ax.transAxes, ha="right", color=MUTED, fontsize=8.5)
    fig.savefig(os.path.join(FIG, "validation_python_vs_matlab.png"), dpi=150)
    plt.close(fig)


def fig_crossshore():
    O, P = load()
    zm, sm = P["zg_model"], P["sm"]
    t_model = [datenum_to_date(t) for t in np.atleast_1d(O.timenum)[:zm.shape[1]]]
    d0 = dt.date(1999, 12, 31)          # first observation row used by project_output
    t_obs = [d0 + dt.timedelta(days=i) for i in range(sm.shape[1])]
    lim = 50
    fig, axs = plt.subplots(2, 1, figsize=(11, 7.4), sharex=True, constrained_layout=True)
    import matplotlib.dates as mdates
    for ax, data, tt, title in ((axs[0], zm, t_model, "ShorelineS (Python port), every 30 days"),
                                (axs[1], sm, t_obs, "CoastSat / PCA reconstruction, 30-day moving mean")):
        x = mdates.date2num(tt)
        im = ax.imshow(data, aspect="auto", cmap=DIVERGING, vmin=-lim, vmax=lim, interpolation="nearest",
                       extent=[x[0], x[-1], data.shape[0] - 0.5, -0.5])
        ax.axhline(293, color=INK, lw=0.8)
        for a, b in NOT_PCA:
            ax.plot([x[0], x[0]], [a, b - 1], color=MUTED, lw=4, solid_capstyle="butt", clip_on=False)
        ax.set_title(title)
        ax.set_ylabel("transect (west → east)")
        ax.grid(False)
        ax.text(x[-1], 140, " beach 2\n nzd0207", va="center", color=INK2, fontsize=9)
        ax.text(x[-1], 450, " beach 3\n nzd0217", va="center", color=INK2, fontsize=9)
    axs[1].xaxis_date()
    cb = fig.colorbar(im, ax=axs, extend="both", shrink=0.8, pad=0.07)
    cb.set_label(f"cross-shore change since start (m)\n+ seaward (accretion), − landward (erosion); clipped at ±{lim} m",
                 color=INK2)
    cb.outline.set_visible(False)
    fig.supxlabel("Grey bars on the left: transects not in the PCA. Their 'observations' are held at one "
                  "CoastSat position, so they show no change by construction.", color=MUTED, fontsize=8.5)
    fig.savefig(os.path.join(FIG, "crossshore_change_model_vs_coastsat.png"), dpi=150)
    plt.close(fig)


def trends(P, O):
    zm, zd = P["zg_model"], P["zg_obs_d"]
    tm = (np.atleast_1d(O.timenum)[:zm.shape[1]] - O.timenum[0]) / 365.25
    td = np.arange(zd.shape[1]) / 365.25

    def fit(z, t):
        out = np.full(z.shape[0], np.nan)
        for i in range(z.shape[0]):
            ok = ~np.isnan(z[i])
            if ok.sum() > 10:
                out[i] = np.polyfit(t[ok], z[i, ok], 1)[0]
        return out
    return fit(zm, tm), fit(zd, td)


def fig_trends(tr_model, tr_obs):
    n = tr_model.size
    i = np.arange(n)
    lim = 8
    fig, ax = plt.subplots(figsize=(11, 3.8), constrained_layout=True)
    shade_not_pca(ax)
    ax.axhline(0, color=AXIS, lw=0.8)
    ax.axvline(293, color=INK, lw=0.8)
    ax.plot(i, np.clip(tr_model, -lim, lim), color=BLUE, label="ShorelineS (Python port)")
    ax.plot(i, np.clip(tr_obs, -lim, lim), color=ORANGE, label="CoastSat / PCA")
    ax.set_ylim(-lim - 0.5, lim + 0.5)
    ax.set_xlim(-1, n)
    ax.set_xlabel("transect (west → east)")
    ax.set_ylabel("trend 2000–2024 (m/yr)\n+ seaward, − landward")
    ax.legend(loc="upper right", ncols=2, fontsize=9)
    ax.set_title("Linear shoreline trend per transect")
    ax.text(146, -lim + 0.2, "beach 2 (nzd0207)", ha="center", color=INK2, fontsize=9)
    ax.text(451, -lim + 0.2, "beach 3 (nzd0217)", ha="center", color=INK2, fontsize=9)
    ax.text(30, lim - 0.3, "not in PCA", ha="center", va="top", color=MUTED, fontsize=8.5)
    ax.text(0.995, -0.2, f"values beyond ±{lim} m/yr are drawn at the limit", transform=ax.transAxes,
            ha="right", color=MUTED, fontsize=8.5)
    fig.savefig(os.path.join(FIG, "shoreline_trends.png"), dpi=150)
    plt.close(fig)


def main():
    os.makedirs(FIG, exist_ok=True)
    O, P = load()
    from plot_waves import plot_waves
    plot_waves(os.path.join(REPO, "wave_data.mat"), os.path.join(FIG, "waves_along_contour.png"))
    fig_runtime()
    fig_validation()
    fig_crossshore()
    tm, to = trends(P, O)
    fig_trends(tm, to)
    with open(os.path.join(REPO, "results", "transect_trends.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["transect", "beach", "in_pca", "model_trend_m_per_yr", "coastsat_trend_m_per_yr",
                    "model_change_2024_m"])
        for k in range(tm.size):
            if k == 293:
                continue
            in_pca = not any(a <= k < b for a, b in NOT_PCA)
            beach = 2 if k < 293 else 3
            local = k + 1 if k < 293 else k - 293
            w.writerow([local, beach, int(in_pca), f"{tm[k]:.3f}", f"{to[k]:.3f}",
                        f"{P['zg_model'][k, -1]:.2f}"])
    print("figures written to", FIG)


if __name__ == "__main__":
    main()
