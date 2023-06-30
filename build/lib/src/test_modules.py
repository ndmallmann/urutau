"""
    Test Urutau
"""

from ._urutau import Urutau
from .modules._filters import ButterworthFilter
from .modules._sn_map import SignalToNoiseMask
from .modules._dereddening import CcmLaw
from .modules._starlight import StarlightOnUrutau
from .modules._resampling import Resampler


def test():

    # Start urutau and set number of threads to be used
    urutau = Urutau()
    urutau.set_number_of_threads(1)

    # Configure modules
    modules = [ButterworthFilter, CcmLaw,
               SignalToNoiseMask, StarlightOnUrutau]
    urutau.set_modules(modules)

    butter_cfg = {"hdu flux": "FLUX"}
    urutau.config_module(ButterworthFilter, butter_cfg)

    ccm_cfg = {"hdu flux": "FLUX_BW"}
    urutau.config_module(CcmLaw, ccm_cfg)

    sn_mask_cfg = {
        "hdu flux": "FLUX_BW",
        "hdu ivar": "IVAR",
        "sn window": (5650, 5750),
        "thresholds": [1, 5, 10, 20]
    }
    urutau.config_module(SignalToNoiseMask, sn_mask_cfg)

    pop_age_par = {
        "xyy": [0, 10.1E6],
        "xyo": [10.1E6, 56.3E6],
        "xiy": [56.3E6, 502.0E6],
        "xii": [502.0E6, 795.0E6],
        "xio": [795.0E6, 2.01E9],
        "xo": [2.01E9, 13E9]
    }
    starlight_cfg = {
        "hdu flux": "FLUX_DRD",
        "hdu flag": "SN_MASKS_10",
        "starlight path": "./starlight/StarlightChains_v04.amd64_gfortran-4.1.1_static.exe",
        "default grid file": "./starlight/reference_grid.in",
        "population ages": pop_age_par,
    }
    urutau.config_module(StarlightOnUrutau, starlight_cfg)

    # Load targets
    urutau.read_csv(targets_dir="./test_data/", csv_file="targets_data.csv")

    # Execute urutau
    urutau.execute("./urutau_results/", save_config=True, debug=False)


def test_sifs():

    # Start urutau and set number of threads to be used
    urutau = Urutau()
    urutau.set_number_of_threads(1)

    # Configure modules
    modules = [ButterworthFilter, CcmLaw, StarlightOnUrutau]
    urutau.set_modules(modules)

    butter_cfg = {"hdu flux": "PRIMARY"}
    urutau.config_module(ButterworthFilter, butter_cfg)

    ccm_cfg = {"hdu flux": "FLUX_BW"}
    urutau.config_module(CcmLaw, ccm_cfg)

    pop_age_par = {
        "xyy": [0, 10.1E6],
        "xyo": [10.1E6, 56.3E6],
        "xiy": [56.3E6, 502.0E6],
        "xii": [502.0E6, 795.0E6],
        "xio": [795.0E6, 2.01E9],
        "xo": [2.01E9, 13E9]
    }
    starlight_cfg = {
        "hdu flux": "FLUX_BW",
        "hdu flag": "SN_MASKS_10",
        "starlight path": "./starlight/StarlightChains_v04.amd64_gfortran-4.1.1_static.exe",
        "default grid file": "./starlight/reference_grid_sifs.in",
        "population ages": pop_age_par,
    }
    urutau.config_module(StarlightOnUrutau, starlight_cfg)

    # Load targets
    urutau.read_csv(targets_dir="./test_data/",
                    csv_file="targets_data_sifs.csv")

    # Execute urutau
    urutau.execute("./sifs_results/", save_config=True, debug=False)


def test_muse():

    # Start urutau and set number of threads to be used
    urutau = Urutau()
    urutau.set_number_of_threads(4)

    # Configure modules
    modules = [Resampler, ButterworthFilter, CcmLaw,
               SignalToNoiseMask, StarlightOnUrutau]
    urutau.set_modules(modules)

    resampler_cfg = {
        "hdu list": ["DATA", "STAT"],
        "data types list": ["flux", "variance"],
        "resampler size": 4
    }
    urutau.config_module(Resampler, resampler_cfg)

    butter_cfg = {"hdu flux": "DATA_RSP"}
    urutau.config_module(ButterworthFilter, butter_cfg)

    ccm_cfg = {"hdu flux": "FLUX_BW"}
    urutau.config_module(CcmLaw, ccm_cfg)

    sn_mask_cfg = {
        "hdu flux": "FLUX_DRD",
        "hdu var": "STAT_RSP",
        "sn window": (5350, 5450),
        "thresholds": [1, 5, 10, 20]
    }
    urutau.config_module(SignalToNoiseMask, sn_mask_cfg)

    pop_age_par = {
        "xyy": [0., 10.1E6],
        "xyo": [10.1E6, 56.3E6],
        "xiy": [56.3E6, 502.0E6],
        "xii": [502.0E6, 795.0E6],
        "xio": [795.0E6, 2.01E9],
        "xo": [2.01E9, 13E9]
    }
    starlight_cfg = {
        "hdu flux": "FLUX_DRD",
        "hdu flag": "SN_MASKS_1",
        "starlight path": "./starlight/StarlightChains_v04.amd64_gfortran-4.1.1_static.exe",
        "default grid file": "./starlight/reference_grid_muse.in",
        "population ages": pop_age_par,
    }
    urutau.config_module(StarlightOnUrutau, starlight_cfg)

    # Load targets
    urutau.read_csv(targets_dir="./test_data/",
                    csv_file="targets_data_muse.csv")

    # Execute urutau
    urutau.execute("./urutau_results_muse/", save_config=True, debug=True)


def quick_test():

    urutau = Urutau()
    urutau.set_number_of_threads(1)

    # Configure modules
    modules = [ButterworthFilter, CcmLaw, SignalToNoiseMask]
    urutau.set_modules(modules)

    butter_cfg = {"hdu flux": "FLUX"}
    urutau.config_module(ButterworthFilter, butter_cfg)

    ccm_cfg = {"hdu flux": "FLUX_BW"}
    urutau.config_module(CcmLaw, ccm_cfg)

    sn_mask_cfg = {
        "hdu flux": "FLUX_DRD",
        "hdu ivar": "IVAR",
        "sn window": (5650, 5750),
        "thresholds": [1, 5, 10, 20]
    }
    urutau.config_module(SignalToNoiseMask, sn_mask_cfg)

    # Load targets
    urutau.read_csv(targets_dir="./test_data/", csv_file="targets_data.csv")

    # Execute urutau
    urutau.execute("./urutau_results/", save_config=True, debug=True)


def quick_test_muse():

    urutau = Urutau()
    urutau.set_number_of_threads(1)

    # Configure modules
    modules = [Resampler,
               ButterworthFilter, CcmLaw, SignalToNoiseMask]
    urutau.set_modules(modules)

    resampler_cfg = {
        "hdu list": ["DATA", "STAT"],
        "data types list": ["flux", "variance"],
        "resampler size": 4
    }
    urutau.config_module(Resampler, resampler_cfg)

    butter_cfg = {"hdu flux": "DATA_RSP"}
    urutau.config_module(ButterworthFilter, butter_cfg)

    ccm_cfg = {"hdu flux": "FLUX_BW"}
    urutau.config_module(CcmLaw, ccm_cfg)

    sn_mask_cfg = {
        "hdu flux": "FLUX_DRD",
        "hdu var": "STAT_RSP",
        "sn window": (5350, 5450),
        "thresholds": [1, 5, 10, 20]
    }
    urutau.config_module(SignalToNoiseMask, sn_mask_cfg)

    # Load targets
    urutau.read_csv(targets_dir="./test_data/",
                    csv_file="targets_data_muse.csv")

    # Execute urutau
    urutau.execute("./urutau_results/", save_config=True, debug=True)
