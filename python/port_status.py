"""Which ShorelineS functions (ShorelineS_functions/*.m) the Python port covers.

    python port_status.py          # rewrites ../docs/shorelines_port_status.md

STATUS maps every MATLAB file to (status, Python location, note):

  full        every code path is ported
  partial     the code paths of the supported options are ported; branches for
              features the port does not have (structures, dunes, ...) are not
  folded      the function's job is done inside other Python code (named)
  not-needed  plotting, video, interactive input or MATLAB utilities: not part
              of the model computation
  no          not ported (a feature the port does not have)

tests/test_shorelines.py checks that every .m file is listed here.
"""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
MATLAB_DIR = os.path.join(REPO, "ShorelineS_functions")

STATUS = {
    # --- main program and set-up
    "ShorelineS": ("partial", "model.ShorelineS.run, model.ShorelineS.step", "time loop of the supported processes"),
    "initialize_defaultvalues": ("partial", "model.DEFAULTS", "defaults of the supported settings"),
    "initialize_time": ("partial", "model.ShorelineS.__init__", "fixed time step only"),
    "initialize_output": ("folded", "model._O_FIELDS", ""),
    "prepare_coastline": ("partial", "model.ShorelineS.__init__", "coastline from S.xmc/S.ymc"),
    "prepare_waveconditions": ("partial", "model.ShorelineS.__init__", "wave time series at stations"),
    "prepare_transport": ("partial", "model.ShorelineS.__init__", "no sediment limitation"),
    "prepare_spit": ("full", "model.ShorelineS.__init__", ""),
    "collect_variables": ("partial", "model.ShorelineS.step, model._collect", "no dune/mud fields"),
    "get_segmentdata": ("partial", "model.ShorelineS._segment", "no dune/mud fields"),
    "save_shorelines": ("partial", "model.ShorelineS._store, model.ShorelineS.save", "no output grids or profiles"),
    "get_timestep": ("folded", "model.ShorelineS.step", "fixed time step only"),
    "get_coastlineupdate_timestep": ("folded", "model.ShorelineS.step, waves.introduce_wave",
                                     "fixed time step only"),
    # --- grid and coastline
    "make_sgrid_mc": ("partial", "coastline.make_sgrid_mc", "griddingmethod 2 only"),
    "get_gridsize": ("no", "", "used by griddingmethod 1"),
    "get_transportpoints": ("full", "coastline.get_transportpoints", ""),
    "get_coastline_orientation": ("full", "coastline.get_coastline_orientation", ""),
    "get_foreshore_orientation": ("partial", "coastline.get_foreshore_orientation",
                                  "S.phif as one angle or [x y angle]; not a file or {method, window}"),
    "get_active_profile": ("partial", "coastline.get_active_profile", "S.d as one value or [x y d]; not a file"),
    "coastline_change": ("partial", "coastline.coastline_change",
                         "no nourishment, sea-level rise, dune or structure terms"),
    "find_overwash_mc": ("partial", "coastline.find_overwash_mc", "spitmethod 'default' only"),
    "merge_coastlines": ("full", "coastline.merge_coastlines", ""),
    "merge_coastlines_mc": ("full", "coastline.merge_coastlines_mc, merge.merge_crossing_sections", ""),
    "cleanup_nans": ("full", "coastline.cleanup_nans", ""),
    "get_smoothdata": ("partial", "coastline.get_smoothdata_angle", "the 'angle' method"),
    "get_one_polygon": ("full", "mfuncs.get_one_polygon", ""),
    "get_one_polyvar": ("full", "mfuncs.get_one_section", ""),
    "insert_section": ("full", "mfuncs.insert_section", ""),
    "get_nansremoved": ("full", "mfuncs.get_nansremoved", ""),
    "get_cumdist": ("full", "mfuncs.get_cumdist", ""),
    "get_cyclic": ("full", "mfuncs.get_cyclic", ""),
    "get_clockpoly": ("full", "mfuncs.get_clockpoly", ""),
    "get_mod": ("full", "mfuncs.get_mod", ""),
    "interpNANs": ("full", "mfuncs.interpNANs", ""),
    "get_intersections": ("full", "geometry.get_intersections", ""),
    "get_polydistance": ("full", "geometry.get_polydistance", ""),
    "find_shadows_mc": ("full", "geometry.find_shadows_mc", ""),
    # --- waves
    "introduce_wave": ("partial", "waves.introduce_wave", "wave time series (S.wvcfile); no wave climates"),
    "get_interpolation_on_grid": ("partial", "waves.interpolate_on_grid",
                                  "'alongshore_mapping' and 'weighted_distance'"),
    "wave_refraction": ("partial", "waves.wave_refraction", "no wave fields (S.waveinteraction)"),
    "wave_breakingheight": ("full", "waves.wave_breakingheight", ""),
    "wave_shoalref": ("folded", "waves._shoalref_vec", "vectorised"),
    "wave_angles": ("partial", "waves.wave_angles", "all transport formulas except RAY"),
    "get_Sphimax": ("full", "waves.get_Sphimax", ""),
    "get_disper": ("full", "mfuncs.get_disper", ""),
    # --- transport
    "transport": ("partial", "transport.transport",
                  "CERC, CERC2, CERC3, KAMP, MILH, VR14; not RAY or TIDEPROF"),
    "transport_shadow_treat": ("partial", "transport.transport_shadow_treat", "no shadow of structures"),
    "get_smoothangles": ("full", "transport.get_smoothangles", ""),
    "get_upwindcorrection": ("full", "transport.get_upwindcorrection", ""),
    "transport_boundary_condition": ("partial", "transport.transport_boundary_condition",
                                     "no groyne boundaries"),
    # --- not needed: plotting, video, input files, utilities
    "initialize_plot": ("not-needed", "", "plotting"),
    "plot_coast": ("not-needed", "", "plotting"),
    "plot_debug": ("not-needed", "", "plotting"),
    "plot_figureproperties": ("not-needed", "", "plotting"),
    "plot_fill_sections": ("not-needed", "", "plotting"),
    "plot_insert_land_fill": ("not-needed", "", "plotting"),
    "plot_profiles": ("not-needed", "", "plotting"),
    "make_video": ("not-needed", "", "video"),
    "extract_shoreline": ("not-needed", "", "figure export at the end of a run"),
    "get_landboundary": ("not-needed", "", "land boundaries for plots"),
    "select_multi_polygon": ("not-needed", "", "interactive input"),
    "readkeys": ("not-needed", "", "keyword input files"),
    "runShorelineS": ("not-needed", "", "starts ShorelineS from an input file"),
    "struct2log": ("not-needed", "", "writes logfile.txt"),
    "isoctave": ("not-needed", "", "MATLAB/Octave check"),
    "get_fields": ("not-needed", "", "MATLAB struct utility"),
    "get_spikesremoved": ("not-needed", "", "not called anywhere in ShorelineS"),
    "get_fallvelocity": ("no", "", "shoreface nourishments"),
    # --- not ported: hard structures and groynes
    "prepare_structures": ("no", "", "hard structures"),
    "prepare_grid_groyne": ("no", "", "groynes (does nothing without structures)"),
    "get_reconnectedgroynes": ("no", "", "groynes"),
    "transport_bypass": ("no", "", "groyne bypassing"),
    "transport_groynesubmerged": ("no", "", "submerged groynes"),
    "transport_revetment": ("no", "", "revetments"),
    "introduce_perm_structures": ("no", "", "permeable structures"),
    "wave_transmission": ("no", "", "transmission over breakwaters"),
    "wave_diffraction": ("no", "", "diffraction"),
    "wave_diffraction_coeff": ("no", "", "diffraction"),
    # --- not ported: other processes
    "prepare_dunes": ("no", "", "dunes"),
    "dune_erosion": ("no", "", "dunes"),
    "dune_flux": ("no", "", "dunes"),
    "dune_growth": ("no", "", "dunes"),
    "extract_berm": ("no", "", "dunes"),
    "prepare_nourishment": ("no", "", "nourishments"),
    "get_nourishments": ("no", "", "nourishments"),
    "prepare_fnourishment": ("no", "", "shoreface nourishments"),
    "get_fnourishment": ("no", "", "shoreface nourishments"),
    "get_fnourishment_diffusion": ("no", "", "shoreface nourishments"),
    "get_fnourishment_fraction": ("no", "", "shoreface nourishments"),
    "get_fnourishment_rate": ("no", "", "shoreface nourishments"),
    "prepare_mudcoast": ("no", "", "mud coasts"),
    "transport_mud": ("no", "", "mud coasts"),
    "get_riverdischarges": ("no", "", "river discharges (mud)"),
    "prepare_channel": ("no", "", "migrating inlet channels"),
    "move_channel": ("no", "", "migrating inlet channels"),
    "find_neighbours": ("no", "", "migrating inlet channels"),
    "get_sectionnr": ("no", "", "migrating inlet channels"),
    "initialize_grid": ("no", "", "migrating inlet channels"),
    "prepare_delta": ("no", "", "deltas"),
    "prepare_tide": ("no", "", "tides"),
    "interpolate_tide": ("no", "", "tides"),
    "tide_1d_ana_anycomp": ("no", "", "tides"),
    "transport_tidewave": ("no", "", "TIDEPROF transport"),
    "transport_soulsbyvanrijn": ("no", "", "TIDEPROF transport"),
    "wave_cur_1D": ("no", "", "TIDEPROF transport"),
    "prepare_climatechange": ("no", "", "climate change (sea-level rise, wave change)"),
    "introduce_climatechange": ("no", "", "climate change (sea-level rise, wave change)"),
    "interpNANsDIR": ("no", "", "climate change and data files"),
    "get_inputfiledata": ("no", "", "time series in data files (climate change)"),
    "prepare_windconditions": ("no", "", "wind (dunes)"),
    "introduce_wind": ("no", "", "wind (dunes)"),
    "prepare_runupconditions": ("no", "", "wave run-up (dunes)"),
    "introduce_runup": ("no", "", "wave run-up (dunes)"),
    "initialize_bathyupdate": ("no", "", "bathymetry update"),
    "update_bathy": ("no", "", "bathymetry update"),
    "update_shoreline": ("no", "", "coastline from a bathymetry grid"),
    "get_disttopolyline": ("no", "", "bathymetry update"),
    "insert_props": ("no", "", "dune/mud properties along the coast"),
    "interpolate_props": ("no", "", "dune/mud properties along the coast"),
    "update_props": ("no", "", "dune/mud properties along the coast"),
    "get_gridprojection": ("no", "", "output grids (S.xyout)"),
    "get_interpolated_wavefield_Hs_dir_Tp": ("no", "", "wave fields (S.waveinteraction)"),
    "get_interpolated_wavefield_dir_Tp": ("no", "", "wave fields (S.waveinteraction)"),
    "get_interpolated_waves_from_grid": ("no", "", "wave fields (S.waveinteraction)"),
    "get_wave_fields_from_mat": ("no", "", "wave fields (S.waveinteraction)"),
    "find_in_grid": ("no", "", "wave fields (S.waveinteraction)"),
    "get_aggregatedtimeseries": ("no", "", "aggregating wave time series (time step > 2x wave interval)"),
    "get_dirspr_prob": ("no", "", "wave climates (directional spreading)"),
    "get_randsample": ("no", "", "wave climates (random sampling)"),
    "initialize_randomgenerator": ("no", "", "wave climates (random sampling)"),
    "get_twisterseed": ("no", "", "wave climates (random sampling)"),
}

LABELS = {"full": "ported", "partial": "ported in part", "folded": "folded into other Python code",
          "not-needed": "not needed (plotting, I/O, utilities)", "no": "not ported"}


def matlab_functions():
    """{name: number of code lines (no comments or blank lines)} of ShorelineS_functions/*.m."""
    out = {}
    for f in sorted(os.listdir(MATLAB_DIR)):
        if f.endswith(".m"):
            with open(os.path.join(MATLAB_DIR, f), encoding="latin-1") as fh:
                out[f[:-2]] = sum(1 for ln in fh if ln.strip() and not ln.lstrip().startswith("%"))
    return out


def summary():
    """Counts and code lines per status, and the share of the model code that is ported."""
    funcs = matlab_functions()
    rows = {}
    for name, lines in funcs.items():
        st = STATUS[name][0]
        n, ln = rows.get(st, (0, 0))
        rows[st] = (n + 1, ln + lines)
    total_n = len(funcs)
    total_l = sum(funcs.values())
    model = {k: v for k, v in rows.items() if k != "not-needed"}
    model_n = sum(n for n, _ in model.values())
    model_l = sum(ln for _, ln in model.values())
    ported = [rows.get(k, (0, 0)) for k in ("full", "partial", "folded")]
    return dict(rows=rows, total_n=total_n, total_l=total_l, model_n=model_n, model_l=model_l,
                ported_n=sum(n for n, _ in ported), ported_l=sum(ln for _, ln in ported))


def markdown():
    funcs = matlab_functions()
    s = summary()
    pct = lambda a, b: f"{100 * a / b:.0f} %"
    out = ["# ShorelineS functions in the Python port", "",
           "Generated by `python/port_status.py` from its `STATUS` table; lines are code lines of the",
           "MATLAB file (no comments or blank lines). See the [Python README](../python/README.md).", "",
           "| Status | Files | Code lines |", "|---|---:|---:|"]
    for st in ("full", "partial", "folded", "not-needed", "no"):
        n, ln = s["rows"].get(st, (0, 0))
        out.append(f"| {LABELS[st]} | {n} | {ln} |")
    out.append(f"| **all** | **{s['total_n']}** | **{s['total_l']}** |")
    out += ["", f"Ported (fully, in part or folded in): **{s['ported_n']} of {s['total_n']} files "
            f"({pct(s['ported_n'], s['total_n'])})**, {s['ported_l']} of {s['total_l']} code lines "
            f"({pct(s['ported_l'], s['total_l'])}). Without plotting, I/O and utilities: "
            f"{s['ported_n']} of {s['model_n']} files ({pct(s['ported_n'], s['model_n'])}), "
            f"{s['ported_l']} of {s['model_l']} code lines ({pct(s['ported_l'], s['model_l'])}).", ""]
    for st in ("full", "partial", "folded", "not-needed", "no"):
        names = [n for n in funcs if STATUS[n][0] == st]
        out += [f"## {LABELS[st].capitalize()} ({len(names)})", "",
                "| MATLAB function | Lines | Python | Note |", "|---|---:|---|---|"]
        for n in names:
            _, py, note = STATUS[n]
            py = ", ".join(f"`{p.strip()}`" for p in py.split(",")) if py else ""
            out.append(f"| `{n}` | {funcs[n]} | {py} | {note} |")
        out.append("")
    return "\n".join(out)


if __name__ == "__main__":
    path = os.path.join(REPO, "docs", "shorelines_port_status.md")
    with open(path, "w") as f:
        f.write(markdown())
    s = summary()
    print(f"{s['ported_n']}/{s['total_n']} files, {s['ported_l']}/{s['total_l']} code lines -> {path}")
