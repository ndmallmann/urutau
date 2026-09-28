"""
script_import.py — Loads an existing run_*.py script (written against the
Urutau scripting API) into a GUI configuration dict.

This is the inverse of script_export.py. Rather than parsing the script's
text, it runs the script with a stand-in Urutau class that records every
add_module() / read_csv() / execute() call instead of actually building and
running a pipeline. Python resolves all of the script's own variables and
expressions for us, so whatever dict ends up passed to add_module() is
captured exactly as the script computed it — no need to re-implement
variable tracking.

The recorded module sequence is then pattern-matched against the pipeline
shape urutau_wrapper.build_module_chain() knows how to build (an optional
paired Spatial stage, an optional paired Resolution stage, an optional
paired Spectral Binning stage, an optional single Dereddening stage, an
optional single S/N Mask stage, an optional single Starlight stage — in
that order). Modules that don't fit this shape are reported back as
warnings rather than silently dropped or causing a crash.

Since this executes the script's own code, only load scripts you trust —
this carries the same trust level as running `python your_script.py`
yourself.
"""

import os
import re

import pandas as pd

import urutau as _urutau_pkg
from urutau.modules import (
    SpatialResampler,
    SpectralResampler,
    DegradeDataFlex,
    StarlightOnUrutau,
)

from .urutau_wrapper import default_config, DEREDDENING_LAWS, SN_MASK_METHODS


class ScriptImportError(Exception):
    """Raised when a script can't be run, or builds no usable pipeline."""


class _RecordingUrutau:
    """
    Stands in for urutau.Urutau while the target script executes: same
    public surface (__init__/add_module/add_target/read_csv/execute), but
    only records what was called instead of doing any real work.
    """

    _last_instance = None

    def __init__(self, num_threads: int = 1) -> None:
        self.num_threads = num_threads
        self.modules = []
        self.targets_dir = None
        self.csv_file = None
        self.execute_kwargs = None
        _RecordingUrutau._last_instance = self

    def add_module(self, module, config=None) -> None:
        self.modules.append((module, dict(config) if config else {}))

    def add_target(self, target, config=None) -> None:
        pass  # individually-added targets aren't reflected in the GUI's targets table

    def read_csv(self, targets_dir, csv_file) -> None:
        self.targets_dir = targets_dir
        self.csv_file = csv_file

    def execute(self, save_path_root="./", save_config=True, debug=False, overwrite=True) -> None:
        self.execute_kwargs = {
            "save_path_root": save_path_root,
            "save_config": save_config,
            "debug": debug,
            "overwrite": overwrite,
        }


def _run_script(path: str) -> "_RecordingUrutau":
    with open(path, "r", encoding="utf-8") as f:
        source = f.read()

    real_urutau_cls = _urutau_pkg.Urutau
    _RecordingUrutau._last_instance = None
    _urutau_pkg.Urutau = _RecordingUrutau
    try:
        namespace = {"__name__": "__main__", "__file__": path}
        try:
            code = compile(source, path, "exec")
            exec(code, namespace)
        except Exception as e:
            raise ScriptImportError(
                f"Error while running '{os.path.basename(path)}': {e}"
            ) from e
    finally:
        _urutau_pkg.Urutau = real_urutau_cls

    recorder = _RecordingUrutau._last_instance
    if recorder is None:
        raise ScriptImportError(
            "The script ran but never created an Urutau instance — nothing to import."
        )
    return recorder


_DEREDDENING_CLASSES = {cls: name for name, cls in DEREDDENING_LAWS.items()}
_SN_MASK_CLASSES = {cls: name for name, (cls, _stat_key) in SN_MASK_METHODS.items()}


def _pop_pair(modules, cls):
    """If the next two recorded calls are both `cls`, consume and return them."""
    if len(modules) >= 2 and modules[0][0] is cls and modules[1][0] is cls:
        pair = (modules[0], modules[1])
        del modules[0:2]
        return pair
    return None


def _pop_single(modules, classes):
    if modules and modules[0][0] in classes:
        item = modules[0]
        del modules[0]
        return item
    return None


def import_script(path: str) -> tuple:
    """
    Runs the script at `path` and returns (cfg, warnings): a GUI
    configuration dict built from what it recorded, and a list of
    human-readable warnings about anything that couldn't be mapped.
    Raises ScriptImportError if the script itself can't be run.
    """
    recorder = _run_script(path)
    cfg = default_config()
    warnings = []

    cfg["num_threads"] = recorder.num_threads

    modules = list(recorder.modules)
    data_hdu = None
    stat_hdu = None

    # --- Spatial resampling -------------------------------------------------
    pair = _pop_pair(modules, SpatialResampler)
    cfg["spatial"]["enabled"] = bool(pair)
    if pair:
        (_, dcfg), (_, scfg) = pair
        cfg["spatial"]["size"] = dcfg.get("resample size", cfg["spatial"]["size"])
        data_hdu = dcfg.get("hdu target")
        stat_hdu = scfg.get("hdu target")

    # --- Resolution degradation ----------------------------------------------
    pair = _pop_pair(modules, DegradeDataFlex)
    cfg["resolution"]["enabled"] = bool(pair)
    if pair:
        (_, dcfg), (_, scfg) = pair
        res = cfg["resolution"]
        res["input_type"] = dcfg.get("input type", res["input_type"])
        res["input_value"] = dcfg.get("input value", res["input_value"])
        res["input_ref_wave"] = dcfg.get("input ref wave")
        res["output_type"] = dcfg.get("output type", res["output_type"])
        res["output_value"] = dcfg.get("output value", res["output_value"])
        res["output_ref_wave"] = dcfg.get("output ref wave")
        if data_hdu is None:
            data_hdu = dcfg.get("hdu target")
            stat_hdu = scfg.get("hdu target")

    # --- Spectral resampling / binning -----------------------------------------
    pair = _pop_pair(modules, SpectralResampler)
    cfg["spectral_bin"]["enabled"] = bool(pair)
    if pair:
        (_, dcfg), (_, scfg) = pair
        cfg["spectral_bin"]["size"] = dcfg.get("resample size", cfg["spectral_bin"]["size"])
        if data_hdu is None:
            data_hdu = dcfg.get("hdu target")
            stat_hdu = scfg.get("hdu target")

    # --- Galactic dereddening -----------------------------------------------
    item = _pop_single(modules, _DEREDDENING_CLASSES)
    cfg["dereddening"]["enabled"] = bool(item)
    if item:
        cls, dcfg = item
        cfg["dereddening"]["law"] = _DEREDDENING_CLASSES[cls]
        cfg["dereddening"]["ebv"] = dcfg.get("ebv", cfg["dereddening"]["ebv"])
        cfg["dereddening"]["rv"] = dcfg.get("rv", cfg["dereddening"]["rv"])
        if data_hdu is None:
            data_hdu = dcfg.get("hdu flux")

    # --- Signal-to-noise mask -------------------------------------------------
    item = _pop_single(modules, _SN_MASK_CLASSES)
    sn_thresholds = list(cfg["sn_mask"]["thresholds"])
    cfg["sn_mask"]["enabled"] = bool(item)
    if item:
        cls, scfg = item
        cfg["sn_mask"]["method"] = _SN_MASK_CLASSES[cls]
        cfg["sn_mask"]["window"] = list(scfg.get("sn window", cfg["sn_mask"]["window"]))
        sn_thresholds = list(scfg.get("thresholds", sn_thresholds))
        cfg["sn_mask"]["thresholds"] = sn_thresholds
        cfg["sn_mask"]["redshift"] = scfg.get("redshift", cfg["sn_mask"]["redshift"])
        for key in ("hdu var", "hdu error", "hdu ivar"):
            if key in scfg:
                stat_hdu = stat_hdu or scfg[key]
        if data_hdu is None:
            data_hdu = scfg.get("hdu flux")

    # --- Starlight -------------------------------------------------------------
    item = _pop_single(modules, {StarlightOnUrutau})
    cfg["starlight"]["enabled"] = bool(item)
    if item:
        _, scfg = item
        sl = cfg["starlight"]
        sl["path"] = scfg.get("starlight path", sl["path"])
        sl["grid_file"] = scfg.get("default grid file", sl["grid_file"])
        if "number of threads" in scfg:
            sl["num_threads"] = scfg["number of threads"]
            sl["auto_threads"] = False  # the script set this explicitly; respect it
        else:
            sl["auto_threads"] = True  # not set in the script; let the GUI auto-balance it
        sl["galaxy_distance"] = scfg.get("galaxy distance", sl["galaxy_distance"])
        sl["redshift"] = scfg.get("redshift", sl["redshift"])
        sl["normalization_factor"] = scfg.get("normalization factor", sl["normalization_factor"])
        sl["flux_unit"] = scfg.get("flux unit", sl["flux_unit"])
        sl["keep_tmp"] = scfg.get("keep tmp", sl["keep_tmp"])
        sl["mask_file"] = scfg.get("mask file") or ""
        sl["timeout_mode"] = scfg.get("timeout mode", sl["timeout_mode"])
        sl["timeout_minutes"] = scfg.get("timeout minutes", sl["timeout_minutes"])
        sl["timeout_window"] = scfg.get("timeout window", sl["timeout_window"])
        sl["timeout_multiplier"] = scfg.get("timeout multiplier", sl["timeout_multiplier"])
        sl["population_ages"] = scfg.get("population ages", {"x": (0, 13e9)})
        sl["sfr_ages"] = scfg.get("sfr ages", {})
        sl["ret_mass_ages"] = scfg.get("ret mass ages", {})
        sl["fc_exps"] = scfg.get("fc exps", {})
        sl["bb_temps"] = scfg.get("bb temps", {})

        flag = scfg.get("hdu flag")
        match = re.match(r"SN_MASKS_(\d+)", str(flag)) if flag else None
        if match:
            sl["flag_threshold"] = int(match.group(1))
        elif sn_thresholds:
            sl["flag_threshold"] = sn_thresholds[0]

        for key in ("hdu var", "hdu error", "hdu ivar"):
            if key in scfg:
                stat_hdu = stat_hdu or scfg[key]
        if data_hdu is None:
            data_hdu = scfg.get("hdu flux")

    if data_hdu:
        cfg["input"]["data_hdu"] = data_hdu
    if stat_hdu:
        cfg["input"]["stat_hdu"] = stat_hdu

    # --- Targets -----------------------------------------------------------
    if recorder.targets_dir is not None:
        cfg["targets"]["dir"] = recorder.targets_dir
    if recorder.csv_file is not None:
        cfg["targets"]["csv"] = recorder.csv_file
        if os.path.isfile(recorder.csv_file):
            try:
                df = pd.read_csv(recorder.csv_file, dtype=str).fillna("")
                cfg["targets"]["columns"] = list(df.columns)
                cfg["targets"]["rows"] = df.values.tolist()
            except Exception as e:
                warnings.append(f"Could not read targets CSV '{recorder.csv_file}': {e}")
        else:
            warnings.append(
                f"Targets CSV '{recorder.csv_file}' was not found on this machine — the path "
                "was imported, but you'll need to add rows by hand or use 'Load CSV…'."
            )
    elif recorder.targets_dir is None and recorder.csv_file is None:
        warnings.append("The script never called urutau.read_csv() — no targets were imported.")

    # --- Output --------------------------------------------------------------
    if recorder.execute_kwargs:
        cfg["output"].update(recorder.execute_kwargs)
    else:
        warnings.append("The script never called urutau.execute() — output settings left at defaults.")

    if modules:
        leftover = ", ".join(sorted({m[0].__name__ for m in modules}))
        warnings.append(
            f"{len(modules)} module call(s) didn't match the pipeline shape the GUI understands "
            f"and were skipped: {leftover}. Configure them by hand if you still need them."
        )

    return cfg, warnings
