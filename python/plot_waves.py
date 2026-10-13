"""Python port of plot_waves.m: wave direction, Hs and Tp along the 25 m contour.

    python plot_waves.py [wave_data.mat] [out.png]
"""
import os
import sys
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)


def load_wave_data(fname):
    try:
        import scipy.io
        d = scipy.io.loadmat(fname, squeeze_me=True, struct_as_record=False)["data"]
        return (np.asarray(d.dates).astype("datetime64[s]"), np.ravel(d.ysave), d.hs_int_m, d.dpm_int_m,
                d.tps_int_m)
    except (NotImplementedError, ValueError):
        import h5py                       # MATLAB v7.3 file written by create_nc.m
        f = h5py.File(fname, "r")
        g = f["data"]
        ysave = g["ysave"][()].ravel()
        # MATLAB datetimes are stored as opaque objects; rebuild the 3-hourly axis
        n = g["hs_int_m"].shape[1]
        dates = np.datetime64("1999-08-01T00:00:00") + np.arange(n) * np.timedelta64(3, "h")
        return dates, ysave, g["hs_int_m"][()].T, g["dpm_int_m"][()].T, g["tps_int_m"][()].T


def plot_waves(fname=os.path.join(REPO, "wave_data.mat"), out=os.path.join(HERE, "output", "waves.png")):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t, ll, h, w, p = load_wave_data(fname)
    fig, axs = plt.subplots(3, 1, figsize=(12, 9), sharex=True, constrained_layout=True)
    from matplotlib.colors import LinearSegmentedColormap
    blues = LinearSegmentedColormap.from_list("blue", ["#f4f8fd", "#9ec5f4", "#3987e5", "#1c5cab", "#0d366b"])
    oranges = LinearSegmentedColormap.from_list("orange", ["#fdf3ee", "#f6b99f", "#eb6834", "#a9441c", "#5c2109"])
    # direction is cyclic (0 = 360 degrees), so it gets a cyclic colour map
    for ax, data, cmap, label in zip(axs, (w.T, h.T, p.T), ("twilight", blues, oranges),
                                     ("Wave direction (°)", "H$_s$ (m)", "T$_p$ (s)")):
        with warnings.catch_warnings():
            # the contour's longitude is not monotonic (MATLAB's pcolor accepts that too)
            warnings.simplefilter("ignore", UserWarning)
            pc = ax.pcolormesh(t, ll, data, shading="auto", cmap=cmap)
        ax.grid(False)
        ax.invert_yaxis()
        ax.set_ylabel("Longitude")
        fig.colorbar(pc, ax=ax, label=label)
    axs[-1].set_xlabel("Timestep = 3 hours")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


if __name__ == "__main__":
    plot_waves(*sys.argv[1:])
