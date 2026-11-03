"""Settings of the model runs that the tests compare with ShorelineS in GNU Octave.

Shared by test_shorelines.py and make_references.py (which runs these cases in
Octave with the ShorelineS submodule and writes tests/data/octave_<case>_output.mat)."""
import os

import numpy as np
import scipy.io

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DATA = os.path.join(HERE, "data")

# [x, y, d] for S.d (case F12): first, middle and last point of the test coastline
_D_XYD = [[1865023.6030439204, 5848190.137671227, 8.0], [1881127.2152420229, 5830048.109426575, 12.0],
          [1900633.8465522733, 5816429.303559587, 10.0]]


def _grid(name):
    c = scipy.io.loadmat(os.path.join(DATA, name))
    return c["xmc"].ravel(), c["ymc"].ravel()


def _synthetic_closed(**kw):
    """2 days on the test coastline (two beaches), closed boundaries, wave stations."""
    x, y = _grid("synthetic_closed_input.mat")
    S = dict(reftime="1999-08-01", endofsimulation="1999-08-03", xmc=x, ymc=y,
             wvcfile=os.path.join(REPO, "input_25_contour.nc"), interpolationmethod="alongshoremapping",
             ddeep=25, dnearshore=25, ds0=75, trform="KAMP", d50=1e-4, boundaryconditionstart="closed",
             boundaryconditionend="closed", d=10, dt=1 / (365 * 8), tc=0, storageinterval=1)
    S.update(kw)
    return S


def _hindcast_periodic():
    r = scipy.io.loadmat(os.path.join(REPO, "30_sept_int_w", "output.mat"), squeeze_me=True,
                         struct_as_record=False)["S"]
    return dict(reftime="2000-01-01", endofsimulation="2000-01-03", xmc=r.xmc, ymc=r.ymc,
                wvcfile=os.path.join(REPO, "input_25_contour.nc"), interpolationmethod="alongshoremapping",
                ddeep=25, dnearshore=25, ds0=75, trform="KAMP", boundaryconditionstart="periodic",
                boundaryconditionend="periodic", d=10, dt=1 / (365 * 8), tc=0, storageinterval=30)


def _cyclic_island():
    x, y = _grid("cyclic_coastline_input.mat")
    return dict(reftime="2000-01-01", endofsimulation="2000-01-07", xmc=x, ymc=y,
                wvcfile=os.path.join(REPO, "input_25_contour.nc"), interpolationmethod="alongshoremapping",
                ddeep=25, dnearshore=25, ds0=75, trform="KAMP", boundaryconditionstart="periodic",
                boundaryconditionend="periodic", d=10, dt=1 / (365 * 8), tc=0, storageinterval=1)


# name -> function returning the settings S (without S.outputdir)
MODEL_CASES = {
    "hindcast_periodic_2day": _hindcast_periodic,
    "cyclic_island_6day": _cyclic_island,
    "synthetic_waves_closed_2day": lambda: _synthetic_closed(wvcfile=os.path.join(REPO, "sinthe_waves.nc"),
                                                             gamma=0.5),
    # smoothrefrac and boundary conditions
    "bc_A": lambda: _synthetic_closed(smoothrefrac=0.5),
    "bc_B": lambda: _synthetic_closed(smoothrefrac=1),
    "bc_C": lambda: _synthetic_closed(boundaryconditionstart="Angleconstant", boundaryconditionend="Angleconstant"),
    "bc_D": lambda: _synthetic_closed(boundaryconditionstart="Fixed", boundaryconditionend="Fixed"),
    "bc_E": lambda: _synthetic_closed(boundaryconditionstart=["Closed", 5000],
                                      boundaryconditionend=["Gradient", 200], smoothrefrac=0.3),
    # options for sensitivity analyses
    "F01": lambda: _synthetic_closed(trform="CERC"),
    "F02": lambda: _synthetic_closed(trform="CERC2"),
    "F03": lambda: _synthetic_closed(trform="CERC3"),
    "F04": lambda: _synthetic_closed(trform="MILH"),
    "F05": lambda: _synthetic_closed(trform="VR14"),
    "F06": lambda: _synthetic_closed(suppresshighangle=1),
    "F07": lambda: _synthetic_closed(relaxationlength=500),
    "F08": lambda: _synthetic_closed(phif=20),
    "F09": lambda: _synthetic_closed(tanbetasetup=0.5),
    "F10": lambda: _synthetic_closed(sphimax="estimate"),
    "F11": lambda: _synthetic_closed(trform="CERC3", qscal=0.3, suppresshighangle=1, relaxationlength=300),
    "F12": lambda: _synthetic_closed(d=np.array(_D_XYD)),
    "F13": lambda: _synthetic_closed(preserveorientation=1),
    # time steps: sub-steps when the stable time step is smaller than S.dt (large
    # transport), and output at given dates (S.storagedate)
    "T01": lambda: _synthetic_closed(qscal=2000, endofsimulation="1999-08-02"),
    "T02": lambda: _synthetic_closed(storagedate=["1999-08-01 06:00:00", "1999-08-02 00:00:00", "1999-08-02 12:00:00"]),
}


def model_settings(name, outputdir):
    S = MODEL_CASES[name]()
    S["outputdir"] = str(outputdir)
    return S
