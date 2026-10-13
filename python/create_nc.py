"""Python port of create_nc.m.

Interpolates the Oceanum wave hindcast (oceanum_paper2.nc, 48 nodes) onto the
25 m depth contour and writes the ShorelineS wave-station file
(input_25_contour.nc format) and wave_data.mat.

By default the outputs go to python/output/ so the MATLAB-generated files in
the repository are not overwritten:

    python create_nc.py [--out-dir DIR] [--plot]
"""
import argparse
import datetime as dt
import os
import time

import contourpy
import netCDF4
import numpy as np
import scipy.io
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator, griddata
from scipy.spatial import Delaunay

from initial_grid import projfwd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)


def colon(a, step, b):
    """Octave a:step:b (used for the 0.01 degree grid)."""
    import sys
    sys.path.insert(0, HERE)
    from shorelines.coastline import colon as _colon
    return _colon(a, step, b)


def contour_level(x, y, Z, level):
    """Contour lines of Z (rows ~ y, cols ~ x) at one level, concatenated in
    the order returned by the contour generator. Returns xs (along x), ys."""
    gen = contourpy.contour_generator(x=x, y=y, z=np.ma.masked_invalid(Z), name="serial",
                                      line_type=contourpy.LineType.Separate, corner_mask=True)
    lines = gen.lines(level)
    # MATLAB's contourc only puts points on grid-cell edges. contourpy's corner
    # masking (cells with one NaN corner) can add points on the cell diagonal;
    # drop those so the stations are the same as in MATLAB.
    def on_grid(v, g):
        return np.min(np.abs(v[:, None] - g[None, :]), axis=1) < 1e-9
    lines = [ln[on_grid(ln[:, 0], np.asarray(x)) | on_grid(ln[:, 1], np.asarray(y))] for ln in lines]
    lines = [ln for ln in lines if len(ln)]
    # longest line first (reproduces the order of MATLAB's contourc for this data set)
    lines = sorted(lines, key=len, reverse=True)
    xs = np.concatenate([ln[:, 0] for ln in lines]) if lines else np.zeros(0)
    ys = np.concatenate([ln[:, 1] for ln in lines]) if lines else np.zeros(0)
    return xs, ys, lines


def interpolation_weights(lon, lat, qlon, qlat):
    """Weights W (n_query x n_nodes) of linear interpolation with nearest-neighbour
    extrapolation (MATLAB scatteredInterpolant(...,'linear','nearest'))."""
    pts = np.column_stack([lon, lat])
    tri = Delaunay(pts)
    q = np.column_stack([qlon, qlat])
    nn = len(lon)
    W = np.zeros((len(qlon), nn))
    eye = np.eye(nn)
    for k in range(nn):
        f = LinearNDInterpolator(tri, eye[:, k])(q)
        W[:, k] = f
    outside = np.isnan(W).any(axis=1)
    if outside.any():
        idx = NearestNDInterpolator(pts, np.arange(nn))(q[outside]).astype(int)
        W[outside] = 0.0
        W[np.flatnonzero(outside), idx] = 1.0
    return W, tri


def interp_fields(W, F, lon, lat, qlon, qlat):
    """W @ F, with per-time-step interpolation where nodes have NaN values."""
    out = W @ np.nan_to_num(F)
    bad = np.flatnonzero(np.isnan(F).any(axis=0))
    if bad.size:
        pts = np.column_stack([lon, lat])
        q = np.column_stack([qlon, qlat])
        for t in bad:
            ok = ~np.isnan(F[:, t])
            lin = LinearNDInterpolator(pts[ok], F[ok, t])(q)
            m = np.isnan(lin)
            if m.any():
                lin[m] = NearestNDInterpolator(pts[ok], F[ok, t])(q[m])
            out[:, t] = lin
    return out


def create_nc(out_dir=os.path.join(HERE, "output"), ncfile=os.path.join(REPO, "oceanum_paper2.nc"),
              level=25.0, dy=0, plot=False, verbose=True):
    t0 = time.perf_counter()
    os.makedirs(out_dir, exist_ok=True)
    with netCDF4.Dataset(ncfile) as d:
        lat = np.asarray(d["latitude"][:], dtype=float)
        lon = np.asarray(d["longitude"][:], dtype=float)
        time_raw = np.asarray(d["time"][:], dtype=float)
        hs = np.asarray(d["hs"][:].filled(np.nan), dtype=float).T      # node x time
        dpm = np.asarray(d["dpm"][:].filled(np.nan), dtype=float).T
        depth = np.asarray(d["botl"][:].filled(np.nan), dtype=float).T
        tps = np.asarray(d["tps"][:].filled(np.nan), dtype=float).T
    t_read = time.perf_counter() - t0

    # select the dates that match the shoreline data
    origin = np.datetime64("1979-02-01T00:00:00")
    times = origin + (time_raw * 3600).astype("timedelta64[s]")
    start = max(np.datetime64("1999-08-01T00:00:00") + np.timedelta64(dy, "D"), times.min())
    end = min(np.datetime64("2024-12-30T00:00:00"), times.max())
    sel = np.flatnonzero((times >= start) & (times <= end))
    dates = times[sel]
    hs_idx, dpm_idx, tps_idx = hs[:, sel], dpm[:, sel], tps[:, sel]

    # 25 m depth contour on a 0.01 degree grid
    llat = colon(lat.min(), 0.01, lat.max())
    llon = colon(lon.min(), 0.01, lon.max())
    LAT, LON = np.meshgrid(llat, llon)
    DEP = griddata(np.column_stack([lat, lon]), depth.mean(axis=1), (LAT, LON), method="linear")
    xsave, ysave, lines = contour_level(llat, llon, DEP, level)   # xsave: latitude, ysave: longitude

    # time in seconds since 1999-08-01
    time_sec = (dates - np.datetime64("1999-08-01T00:00:00")).astype("timedelta64[s]").astype(float)
    station_x, station_y = projfwd(xsave, ysave)

    # interpolation (linear, nearest-neighbour extrapolation) via a weight matrix
    W, _ = interpolation_weights(lon, lat, ysave, xsave)
    hs_int = interp_fields(W, hs_idx, lon, lat, ysave, xsave)
    u_int = interp_fields(W, np.cos(np.deg2rad(dpm_idx)), lon, lat, ysave, xsave)
    v_int = interp_fields(W, np.sin(np.deg2rad(dpm_idx)), lon, lat, ysave, xsave)
    dpm_int = np.mod(np.rad2deg(np.arctan2(v_int, u_int)), 360)
    tps_int = interp_fields(W, tps_idx, lon, lat, ysave, xsave)

    # write the ShorelineS wave-station file
    fn = os.path.join(out_dir, "input_25_contour.nc")
    if os.path.exists(fn):
        os.remove(fn)
    with netCDF4.Dataset(fn, "w") as d:
        d.createDimension("time", time_sec.size)
        d.createDimension("station", xsave.size)
        v = d.createVariable("time", "f8", ("time",))
        v[:] = time_sec
        v.units = "seconds since 1999-08-01 00:00:00"
        v.long_name = "time"
        d.createVariable("station_x", "f8", ("station",))[:] = station_x
        d.createVariable("station_y", "f8", ("station",))[:] = station_y
        for name, data, unit in (("point_hm0", hs_int, "m"), ("point_wavdir", dpm_int, "degrees"),
                                 ("point_tp", tps_int, "s")):
            v = d.createVariable(name, "f8", ("station", "time"))
            v[:] = data
            v.units = unit
    scipy.io.savemat(os.path.join(out_dir, "wave_data.mat"),
                     {"data": {"dates": dates.astype("datetime64[s]").astype(str), "xsave": xsave,
                               "ysave": ysave, "dpm_int_m": dpm_int.T, "hs_int_m": hs_int.T,
                               "tps_int_m": tps_int.T}}, do_compression=False)
    if plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 6))
        pc = ax.pcolormesh(llon, llat, DEP.T if DEP.shape != (llon.size, llat.size) else DEP.T, shading="auto",
                           cmap="viridis_r")
        fig.colorbar(pc, label="depth (m)")
        ax.plot(ysave, xsave, "k.", ms=2)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        fig.savefig(os.path.join(out_dir, "contour.png"), dpi=150)
    elapsed = time.perf_counter() - t0
    if verbose:
        print(f"create_nc: {xsave.size} stations, {time_sec.size} times, read {t_read:.1f}s, total {elapsed:.1f}s")
    return dict(xsave=xsave, ysave=ysave, station_x=station_x, station_y=station_y, time_sec=time_sec,
                hs_int=hs_int, dpm_int=dpm_int, tps_int=tps_int, DEP=DEP, llat=llat, llon=llon,
                lines=lines, elapsed=elapsed)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default=os.path.join(HERE, "output"))
    ap.add_argument("--plot", action="store_true")
    a = ap.parse_args()
    create_nc(a.out_dir, plot=a.plot)
