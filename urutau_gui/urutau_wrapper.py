"""
urutau_wrapper.py — Translates a GUI configuration dict into an Urutau
pipeline (module chain + targets) and runs it.

The configuration dict mirrors exactly what MainWindow._gather_config()
produces / MainWindow._apply_config() consumes, so the same dict can be
saved/loaded as JSON and reused from scripts.
"""

from urutau import Urutau
from urutau.modules import (
    SpatialResampler,
    SpectralResampler,
    DegradeDataFlex,
    CcmLaw,
    SeatonLaw,
    SNMaskWithVar,
    SNMaskWithError,
    SNMaskWithIVar,
    SNMaskMeanStd,
    StarlightOnUrutau,
)
from urutau.modules._dereddening import CalzettiLaw, Fitzpatrick

DEREDDENING_LAWS = {
    "CCM": CcmLaw,
    "Seaton": SeatonLaw,
    "Calzetti": CalzettiLaw,
    "Fitzpatrick": Fitzpatrick,
}

# label -> (module class, stat-hdu config key or None)
SN_MASK_METHODS = {
    "Variance": (SNMaskWithVar, "hdu var"),
    "Error": (SNMaskWithError, "hdu error"),
    "Inverse Variance": (SNMaskWithIVar, "hdu ivar"),
    "Mean / Std (flux only)": (SNMaskMeanStd, None),
}

# label used in the GUI -> "data type" string understood by resampling/
# resolution modules
STAT_TYPE_MAP = {
    "variance": "variance",
    "error": "error",
    "inverse variance": "inv_variance",
}


def default_config() -> dict:
    """A default configuration matching run_newurutau_bass27_upto8850.py."""
    return {
        "num_threads": 1,
        "input": {
            "data_hdu": "DATA",
            "stat_hdu": "STAT",
            "stat_type": "variance",
        },
        "spatial": {"enabled": False, "size": 4},
        "resolution": {
            "enabled": True,
            "input_type": "R", "input_value": 3027.0, "input_ref_wave": None,
            "output_type": "FWHM", "output_value": 2.51, "output_ref_wave": None,
        },
        "spectral_bin": {"enabled": True, "size": 1.0},
        "dereddening": {"enabled": True, "law": "CCM", "ebv": 0.0, "rv": 3.1},
        "sn_mask": {
            "enabled": True,
            "method": "Variance",
            "window": (5650.0, 5750.0),
            "thresholds": [10, 20],
            "redshift": 0.0,
        },
        "starlight": {
            "enabled": True,
            "path": "./starlight/StarlightChains_v04.amd64_g77-3.4.6-r1_static.exe",
            "grid_file": "./starlight/reference_grid_muse_newMiles.in",
            "num_threads": 52,
            "flag_threshold": 10,
            "galaxy_distance": 0.0,
            "redshift": 0.0,
            "normalization_factor": 1.0,
            "flux_unit": "",
            "keep_tmp": False,
            "mask_file": "",
            "timeout_mode": "none",
            "timeout_minutes": None,
            "timeout_window": 15,
            "timeout_multiplier": 2.0,
            "population_ages": {
                "xyy": (100, 10.1e6), "xyo": (10.1e6, 56.3e6),
                "xiy": (56.3e6, 502.0e6), "xy": (100, 56.3e6),
                "xii": (502.0e6, 795.0e6), "xio": (795.0e6, 2.01e9),
                "xi": (56.3e6, 2.01e9), "xo": (2.01e9, 13e9),
            },
            "sfr_ages": {},
            "ret_mass_ages": {"inter": (10.0, 2.0e9)},
            "fc_exps": {"FC_50": (0.49, 0.51)},
            "bb_temps": {},
        },
        "targets": {"dir": "", "csv": "", "columns": ["target"], "rows": []},
        "output": {
            "save_path_root": "./urutau_output/", "save_config": True, "debug": True,
            "overwrite": True,
        },
    }


def build_module_chain(cfg: dict) -> list:
    """
    Returns the enabled processing stages as [(ModuleClass, config_dict), ...],
    in execution order, with HDU names chained exactly like build_pipeline()
    would add them to a live Urutau instance — but without needing one (or
    any target files to exist on disk). Used both by build_pipeline() and by
    the standalone-script exporter.
    """
    chain = []

    data_hdu = cfg["input"]["data_hdu"]
    stat_hdu = cfg["input"]["stat_hdu"]
    stat_type = STAT_TYPE_MAP[cfg["input"]["stat_type"]]

    # --- Spatial resampling -------------------------------------------------
    spatial = cfg["spatial"]
    if spatial["enabled"]:
        size = spatial["size"]
        chain.append((SpatialResampler, {
            "hdu target": data_hdu, "data type": "flux", "resample size": size,
        }))
        chain.append((SpatialResampler, {
            "hdu target": stat_hdu, "data type": stat_type, "resample size": size,
        }))
        data_hdu = f"{data_hdu}_RSP"
        stat_hdu = f"{stat_hdu}_RSP"

    # --- Spectral resolution degradation ------------------------------------
    res = cfg["resolution"]
    if res["enabled"]:
        def _res_cfg(target):
            c = {
                "hdu target": target,
                "input type": res["input_type"], "input value": res["input_value"],
                "output type": res["output_type"], "output value": res["output_value"],
            }
            if res.get("input_ref_wave") is not None:
                c["input ref wave"] = res["input_ref_wave"]
            if res.get("output_ref_wave") is not None:
                c["output ref wave"] = res["output_ref_wave"]
            return c

        data_cfg = _res_cfg(data_hdu)
        data_cfg["data type"] = "flux"
        stat_cfg = _res_cfg(stat_hdu)
        stat_cfg["data type"] = stat_type

        chain.append((DegradeDataFlex, data_cfg))
        chain.append((DegradeDataFlex, stat_cfg))
        data_hdu = f"{data_hdu}_DEGR"
        stat_hdu = f"{stat_hdu}_DEGR"

    # --- Spectral resampling / binning --------------------------------------
    sbin = cfg["spectral_bin"]
    if sbin["enabled"]:
        size = sbin["size"]
        chain.append((SpectralResampler, {
            "hdu target": data_hdu, "data type": "flux", "resample size": size,
        }))
        chain.append((SpectralResampler, {
            "hdu target": stat_hdu, "data type": stat_type, "resample size": size,
        }))
        data_hdu = f"{data_hdu}_BIN"
        stat_hdu = f"{stat_hdu}_BIN"

    # --- Galactic dereddening ------------------------------------------------
    drd = cfg["dereddening"]
    if drd["enabled"]:
        law_cls = DEREDDENING_LAWS[drd["law"]]
        chain.append((law_cls, {
            "hdu flux": data_hdu, "ebv": drd["ebv"], "rv": drd["rv"],
        }))
        data_hdu = "FLUX_DRD"

    # --- Signal-to-noise mask -------------------------------------------------
    sn = cfg["sn_mask"]
    if sn["enabled"]:
        mask_cls, stat_key = SN_MASK_METHODS[sn["method"]]
        sn_cfg = {
            "hdu flux": data_hdu,
            "sn window": tuple(sn["window"]),
            "thresholds": list(sn["thresholds"]),
            "redshift": sn["redshift"],
        }
        if stat_key is not None:
            sn_cfg[stat_key] = stat_hdu
        chain.append((mask_cls, sn_cfg))

    # --- Starlight -------------------------------------------------------------
    sl = cfg["starlight"]
    if sl["enabled"]:
        starlight_cfg = {
            "starlight path": sl["path"],
            "default grid file": sl["grid_file"],
            "hdu flux": data_hdu,
            "number of threads": sl["num_threads"],
            "population ages": sl["population_ages"],
            "galaxy distance": sl["galaxy_distance"],
            "redshift": sl["redshift"],
            "normalization factor": sl["normalization_factor"],
            "flux unit": sl["flux_unit"],
            "keep tmp": sl["keep_tmp"],
            "timeout mode": sl.get("timeout_mode", "none"),
            "timeout minutes": sl.get("timeout_minutes"),
            "timeout window": sl.get("timeout_window", 15),
            "timeout multiplier": sl.get("timeout_multiplier", 2.0),
        }
        if sl["sfr_ages"]:
            starlight_cfg["sfr ages"] = sl["sfr_ages"]
        if sl["ret_mass_ages"]:
            starlight_cfg["ret mass ages"] = sl["ret_mass_ages"]
        if sl["fc_exps"]:
            starlight_cfg["fc exps"] = sl["fc_exps"]
        if sl["bb_temps"]:
            starlight_cfg["bb temps"] = sl["bb_temps"]

        if sn["enabled"]:
            _, stat_key = SN_MASK_METHODS[sn["method"]]
            if stat_key is not None:
                starlight_cfg[stat_key] = stat_hdu
            if sl.get("flag_threshold") is not None:
                starlight_cfg["hdu flag"] = f"SN_MASKS_{int(sl['flag_threshold'])}"

        if sl.get("mask_file"):
            starlight_cfg["mask file"] = sl["mask_file"]

        chain.append((StarlightOnUrutau, starlight_cfg))

    return chain


def build_pipeline(cfg: dict) -> Urutau:
    urutau = Urutau(num_threads=int(cfg.get("num_threads", 1)))

    for module_cls, module_cfg in build_module_chain(cfg):
        urutau.add_module(module_cls, module_cfg)

    # --- Targets -----------------------------------------------------------
    targets = cfg["targets"]
    urutau.read_csv(targets_dir=targets["dir"], csv_file=targets["csv"])

    return urutau


def run_pipeline(cfg: dict) -> None:
    """Builds and executes the pipeline described by cfg (blocking call)."""
    urutau = build_pipeline(cfg)
    out = cfg["output"]
    urutau.execute(
        save_path_root=out["save_path_root"],
        save_config=out["save_config"],
        debug=out["debug"],
        overwrite=out.get("overwrite", True),
    )
