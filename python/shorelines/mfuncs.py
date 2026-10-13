"""MATLAB/Octave-compatible numerical helpers.

The ShorelineS model relies on details of MATLAB/Octave numerics (``mod``,
``sind``/``cosd`` with exact multiples of 90 degrees, linear ``interp1`` with
NaN outside the data range, ``get_one_polygon`` sections separated by NaN, ...).
These helpers reproduce the GNU Octave 8 implementations so that the Python
port gives (nearly) the same floating-point results as the MATLAB/Octave code.
"""
import ctypes as _ctypes
import ctypes.util as _ctypes_util
import datetime as _dt
import os as _os

import numpy as np

__all__ = [
    "EXACT_LIBM", "m_sin", "m_exp", "m_sinh", "m_asin", "m_acos", "m_atan2",
    "m_hypot", "m_pow", "mmod", "sind", "cosd", "asind", "acosd", "atan2d", "interp1",
    "interp1_nearest_extrap", "interpNANs", "datenum", "datestr",
    "get_one_polygon", "insert_section", "get_clockpoly", "get_cyclic",
    "get_disper", "get_mod", "get_cumdist", "get_nansremoved",
]

_PI = np.pi

# ---------------------------------------------------------------------------
# Elementary functions.
# NumPy's SIMD implementations of pow/exp/sinh/asin/acos/atan2/hypot can differ
# from glibc's libm (used by Octave) in the last bit. By default the fast NumPy
# versions are used. Set the environment variable SHORELINES_EXACT_LIBM=1
# (before importing this package) to call glibc's libm element by element
# instead; this is slow but reproduces Octave's floating-point results, which
# is used to verify the port.
# ---------------------------------------------------------------------------
EXACT_LIBM = _os.environ.get("SHORELINES_EXACT_LIBM", "0") == "1"


def _libm_ufunc(name, nin):
    libm = _ctypes.CDLL(_ctypes_util.find_library("m"))
    f = getattr(libm, name)
    f.restype = _ctypes.c_double
    f.argtypes = [_ctypes.c_double] * nin
    uf = np.frompyfunc(f, nin, 1)

    def call(*args):
        out = uf(*[np.asarray(a, dtype=float) for a in args])
        return np.asarray(out, dtype=float)
    return call


if EXACT_LIBM:
    m_sin = _libm_ufunc("sin", 1)
    m_exp = _libm_ufunc("exp", 1)
    m_sinh = _libm_ufunc("sinh", 1)
    m_asin = _libm_ufunc("asin", 1)
    m_acos = _libm_ufunc("acos", 1)
    m_atan2 = _libm_ufunc("atan2", 2)
    m_hypot = _libm_ufunc("hypot", 2)
    _libm_pow = _libm_ufunc("pow", 2)

    def m_pow(x, p):
        """Octave ``x.^p`` (x*x for p == 2, libm pow otherwise)."""
        if np.isscalar(p) and p == 2:
            return np.asarray(x, dtype=float) * np.asarray(x, dtype=float)
        return _libm_pow(x, p)
else:
    m_sin = np.sin
    m_exp = np.exp
    m_sinh = np.sinh
    m_asin = np.arcsin
    m_acos = np.arccos
    m_atan2 = np.arctan2
    m_hypot = np.hypot

    def m_pow(x, p):
        return np.power(np.asarray(x, dtype=float), p)


def mmod(x, y):
    """Octave ``mod(x, y)`` for scalar ``y`` (x - floor(x/y)*y, exact zero for
    integer quotients, result carries the sign of y)."""
    x = np.asarray(x, dtype=float)
    if y == 0:
        return x.copy()
    with np.errstate(invalid="ignore"):
        q = x / y
        n = np.floor(q)
        r = x - y * n
        r = np.where(np.rint(q) == q, 0.0, r)
        if y > 0:
            r = np.where(r < 0, -r, r)
        else:
            r = np.where(r > 0, -r, r)
    return r


def sind(x):
    """Octave ``sind``: sine of an angle in degrees (exact zeros at multiples of 180)."""
    x = mmod(np.asarray(x, dtype=float) - 180.0, 360.0) - 180.0
    y = m_sin(x / 180.0 * _PI)
    return np.where(x == -180.0, 0.0, y)


def cosd(x):
    """Octave ``cosd`` = ``sind(x + 90)``."""
    return sind(np.asarray(x, dtype=float) + 90.0)


def asind(x):
    with np.errstate(invalid="ignore"):
        return m_asin(x) * 180.0 / _PI


def acosd(x):
    with np.errstate(invalid="ignore"):
        return m_acos(x) * 180.0 / _PI


def atan2d(y, x):
    return 180.0 / _PI * m_atan2(y, x)


def interp1(x, y, xi):
    """Linear ``interp1(x, y, xi)`` as in Octave: NaN outside [min(x), max(x)].

    ``y`` may be 1-D (len(x)) or 2-D with interpolation along axis 0.
    Unsorted ``x`` is sorted first (as Octave does)."""
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float)
    xi_arr = np.asarray(xi, dtype=float)
    scalar = xi_arr.ndim == 0
    xi_flat = xi_arr.ravel()
    if y.ndim == 1:
        y = y.reshape(-1)
    if x.size >= 2 and np.any(np.diff(x) < 0):
        order = np.argsort(x, kind="stable")
        x = x[order]
        y = y[order]
    n = x.size
    if n < 2:
        raise ValueError("interp1: minimum of 2 points required")
    dx = np.diff(x)
    dy = np.diff(y, axis=0)
    if y.ndim == 1:
        slope = dy / dx
    else:
        slope = dy / dx[:, None]
    # lookup(x, xi, "lr"): last index with x(idx) <= xi, clamped to [0, n-2]
    idx = np.searchsorted(x, xi_flat, side="right") - 1
    idx = np.clip(idx, 0, n - 2)
    d = xi_flat - x[idx]
    if y.ndim == 1:
        yi = slope[idx] * d + y[idx]
    else:
        yi = slope[idx] * d[:, None] + y[idx]
    out = (xi_flat < x[0]) | ~(xi_flat <= x[-1])
    if np.any(out):
        yi = yi.copy()
        yi[out] = np.nan
    if scalar:
        return yi[0]
    if y.ndim == 1:
        return yi.reshape(xi_arr.shape)
    return yi


def interp1_nearest_extrap(x, y, xi):
    """``interp1(x, y, xi, 'nearest', 'extrap')`` as in Octave."""
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    xi = np.asarray(xi, dtype=float)
    if x.size >= 2 and np.any(np.diff(x) < 0):
        order = np.argsort(x, kind="stable")
        x = x[order]
        y = y[order]
    mid = 0.5 * (x[:-1] + x[1:])
    idx = np.searchsorted(mid, xi.ravel(), side="right")
    return y[idx].reshape(xi.shape)


def interpNANs(data):
    """ShorelineS ``interpNANs`` for a vector: fill NaNs by linear
    interpolation and nearest-neighbour extrapolation."""
    data = np.asarray(data, dtype=float)
    shape = data.shape
    y = data.ravel()
    x = np.arange(1, y.size + 1, dtype=float)
    ok = ~np.isnan(y)
    if ok.all():
        return data.copy()
    x2 = x[ok]
    y2 = y[ok]
    if x2.size > 1:
        y3 = interp1(x2, y2, x)
        ok3 = ~np.isnan(y3)
        ynew = interp1_nearest_extrap(x[ok3], y3[ok3], x)
    elif x2.size == 1:
        ynew = np.full(x.shape, y2[0])
    else:
        ynew = np.full(x.shape, np.nan)
    return ynew.reshape(shape)


_DATENUM_OFFSET = 366  # datenum(0000-01-00) vs. Python proleptic ordinal


def datenum(s):
    """MATLAB datenum of an ISO date string 'yyyy-mm-dd[ HH:MM:SS]'."""
    s = s.strip()
    if " " in s or "T" in s:
        t = _dt.datetime.fromisoformat(s.replace("T", " "))
    else:
        t = _dt.datetime.fromisoformat(s + " 00:00:00")
    day = t.toordinal() + _DATENUM_OFFSET
    frac = (t.hour * 3600 + t.minute * 60 + t.second) / 86400.0
    return day + frac


def datestr(dn, fmt="%Y-%m-%d %H:%M"):
    day = int(np.floor(dn))
    frac = dn - day
    t = _dt.datetime.fromordinal(day - _DATENUM_OFFSET) + _dt.timedelta(seconds=round(frac * 86400))
    return t.strftime(fmt)


def get_one_polygon(x_mc, y_mc, i_mc):
    """Section ``i_mc`` (1-based) of NaN-separated polylines.

    Returns x, y, n_mc, i1, i2 (i1, i2 1-based like MATLAB)."""
    x_mc = np.asarray(x_mc, dtype=float).ravel()
    y_mc = np.asarray(y_mc, dtype=float).ravel()
    nans = np.flatnonzero(np.isnan(x_mc))
    if nans.size > 1:
        nansremove = nans[1:][np.diff(nans) == 1]
        if nansremove.size:
            keep = np.ones(x_mc.size, bool)
            keep[nansremove] = False
            x_mc = x_mc[keep]
            y_mc = y_mc[keep]
            nans = np.flatnonzero(np.isnan(x_mc))
    nans1 = nans + 1  # 1-based
    n_mc = nans.size + 1
    if nans.size == 0:
        i1, i2 = 1, x_mc.size
    elif i_mc == 1:
        i1, i2 = 1, nans1[0] - 1
    elif i_mc == n_mc:
        i1, i2 = nans1[i_mc - 2] + 1, x_mc.size
    elif i_mc > n_mc:
        i1, i2 = nans1[n_mc - 2] + 1, x_mc.size
    else:
        i1, i2 = nans1[i_mc - 2] + 1, nans1[i_mc - 1] - 1
    x = x_mc[i1 - 1:i2].copy()
    y = y_mc[i1 - 1:i2].copy()
    if not np.any(~np.isnan(x_mc)):
        n_mc = 0
    return x, y, n_mc, int(i1), int(i2)


def get_one_section(v_mc, i_mc):
    """``get_one_polygon(v_mc, i_mc)`` with 2 arguments: one NaN-separated vector."""
    x, _, _, _, _ = get_one_polygon(v_mc, v_mc, i_mc)
    return x


def insert_section(xnew, ynew, x_mc, y_mc, i_mc):
    """Replace section ``i_mc`` of (x_mc, y_mc) by (xnew, ynew)."""
    _, _, _, i1, i2 = get_one_polygon(x_mc, y_mc, i_mc)
    x_mc = np.asarray(x_mc, dtype=float).ravel()
    y_mc = np.asarray(y_mc, dtype=float).ravel()
    xnew = np.atleast_1d(np.asarray(xnew, dtype=float)).ravel()
    ynew = np.atleast_1d(np.asarray(ynew, dtype=float)).ravel()
    xo = np.concatenate([x_mc[:i1 - 1], xnew, x_mc[i2:]])
    yo = np.concatenate([y_mc[:i1 - 1], ynew, y_mc[i2:]])
    return xo, yo


def insert_section1(vnew, v_mc, i_mc):
    """``insert_section(vnew, v_mc, i_mc)`` with 3 arguments (one vector)."""
    xo, _ = insert_section(vnew, vnew, v_mc, v_mc, i_mc)
    return xo


def get_clockpoly(x, y):
    """Orientation of a polygon: 1 clockwise, -1 counter-clockwise, 0 undetermined."""
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    if x.size == 0 or y.size == 0:
        return 0
    xc = np.mean(x)
    xd = np.diff(np.concatenate([x, x[:1]]))
    yd = np.diff(np.concatenate([y, y[:1]]))
    mark = xd == 0
    xd = xd.copy()
    xd[mark] = np.nan
    with np.errstate(invalid="ignore", divide="ignore"):
        a = (xc - x) / xd
        mark2 = (a < 0) | (a > 1)
        yc = y + a * yd
    yc[mark & (xc != x)] = -np.inf
    yc[mark2] = -np.inf
    mark = mark & (xc == x)
    yc[mark] = y[mark] + np.maximum(0, yd[mark])
    ycm = np.nanmax(yc) if np.any(~np.isnan(yc)) else np.nan
    i = np.flatnonzero(yc == ycm)
    xdi = xd[i]
    if np.all(xdi >= 0) and np.any(xdi > 0):
        return 1
    if np.all(xdi <= 0) and np.any(xdi < 0):
        return -1
    return 0


def get_cyclic(x, y, ds0):
    frac = 0.2
    return bool(m_hypot(x[-1] - x[0], y[-1] - y[0]) < frac * ds0)


def get_disper(h, T):
    """Wave dispersion (Hunt approximation): k, C, Cg, n."""
    g = 9.81
    sigma = 2 * _PI / np.asarray(T, dtype=float)
    h = np.asarray(h, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        k = sigma * sigma / g * m_pow(1 - m_exp(-m_pow(sigma * np.sqrt(h / g), 2.5)), -0.4)
        C = sigma / k
        n = 0.5 + k * h / m_sinh(2 * k * h)
        Cg = n * C
    return k, C, Cg, n


def get_mod(x, y):
    """ShorelineS ``get_mod``: mod with 0 replaced by y (1-based cyclic index)."""
    z = int(x % y)
    return y if z == 0 else z


def get_cumdist(x, y):
    return np.concatenate([[0.0], np.cumsum(m_hypot(np.diff(x), np.diff(y)))])


def get_nansremoved(x_mc, y_mc):
    x_mc = np.asarray(x_mc, dtype=float).copy()
    y_mc = np.asarray(y_mc, dtype=float).copy()
    nans = np.flatnonzero(np.isnan(x_mc))
    if nans.size > 1:
        rem = nans[:-1][np.diff(nans) == 1]
        x_mc[rem] = -2e10
        y_mc[rem] = -2e10
    x_mc = x_mc[x_mc != -2e10]
    y_mc = y_mc[y_mc != -2e10]
    if x_mc.size and np.isnan(x_mc[0]):
        x_mc, y_mc = x_mc[1:], y_mc[1:]
    if x_mc.size and np.isnan(x_mc[-1]):
        x_mc, y_mc = x_mc[:-1], y_mc[:-1]
    return x_mc, y_mc
