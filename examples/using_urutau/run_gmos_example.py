"""
    Test Urutau
"""

from urutau import Urutau

from urutau.modules import (
    SpatialResampler,
    SpectralResampler,
    DegradeDataFlex,
    SNMaskMeanStd,
    StarlightOnUrutau,
)



def quick_resampling_gmos():
    
    # Start urutau with 3 threads
    urutau = Urutau(num_threads = 12)

    # Spatial resampler — commented out (not needed for this datacube).
    # Uncomment if spatial binning is required.
    #
    # data_spat_resample_cfg = {
    #     "hdu target": "SCI",
    #     "data type": "flux",
    #     "resample size": 4
    # }
    # urutau.add_module(SpatialResampler, data_spat_resample_cfg)

    # Spectral resampler for the data hdu
    data_spec_resample_cfg = {
        "hdu target": "SCI",
        "data type": "flux",
        "resample size": 1.
    }
    urutau.add_module(SpectralResampler, data_spec_resample_cfg)

    # Degrade from GMOS (R=2260 at λ=5620 Å) to MILES library resolution (FWHM=2.5 Å)
    # "input ref wave" pins sigma to the value at 5620 Å: sigma = 5620/(2260×2.35482) ≈ 1.06 Å
    data_degrade_cfg = {
        "hdu target": "SCI_RSP",
        "data type": "flux",
        "input type":     "R",    "input value":    2260, "input ref wave": 5620,
        "output type":    "FWHM", "output value":   2.5,  # MILES library FWHM in Angstrom
    }
    urutau.add_module(DegradeDataFlex, data_degrade_cfg)

    # No error/variance/ivar HDU available in this datacube — commented out.
    # Uncomment and adjust if a STAT extension becomes available.
    #
    # data_degrade_cfg = {
    #     "hdu target": "STAT_RSP_BIN",
    #     "data type": "variance",
    #     "input type": "R",    "input value": 2200,
    #     "output type": "FWHM", "output value": 2.5,
    # }
    # urutau.add_module(DegradeDataFlex, data_degrade_cfg)

    # Signal To Noise mask using mean/std of the flux window (no error HDU needed)
    sn_mask_cfg = {
        "hdu flux": "SCI_RSP_DEGR",
        "sn window": [5650, 5750],
        "thresholds": [1, 5, 10]
    }
    urutau.add_module(SNMaskMeanStd, sn_mask_cfg)

    # Run Starlight
    population_ages = {
        "xyy": (0., 1E7),
        "xyo": (1E7, 5.6E7),
        "xiy": (5.6E7, 5E8),
        "xii": (5E8, 8E8),
        "xio": (8E8, 2E9),
        "xo": (2E9, 13E9)
    }
    starlight_cfg = {
        "starlight path": "./starlight/StarlightChains_v04.amd64_g77-3.4.6-r1_static.exe",
        "default grid file": "./starlight/reference_grid_gmos.in",
        "hdu flux": "SCI_RSP_DEGR",
        "hdu flag": "SN_MASKS_1",
        "population ages": population_ages,
    }
    urutau.add_module(StarlightOnUrutau, starlight_cfg)

    # Load targets
    urutau.read_csv(targets_dir="./datacubes/",
                    csv_file="./datacubes/gmos_target.csv")

    # Execute urutau
    urutau.execute("./gmos_fit/", save_config=True, debug=True)

if __name__=="__main__":
    quick_resampling_gmos()
