"""
    Urutau to run starlight on JWST NIRSPEC datacubes
"""

from urutau import Urutau

from urutau.modules import (
    SpatialResampler,
    SpectralResampler,
    DegradeData,
    SNMaskWithVar,
    StarlightOnUrutau,
    SNMaskWithError,
)



def nirspec_JWST():
    
    # Start urutau with 3 threads
    urutau = Urutau(num_threads = 1)

    # First spatial resampler for the data hdu (not used, if you whant to use it change the HDU names acordingly on the next steps)
#    data_spat_resample_cfg = {
#        "hdu target": "SCI",
#        "data type": "flux",
#        "resample size": 2
#    }
#    urutau.add_module(SpatialResampler, data_spat_resample_cfg)

#    # Second spatial resampler for the stat hdu (not used, if you whant to use it change the HDU names acordingly on the next steps)
#    stat_spat_resample_cfg = {
#        "hdu target": "ERR",
#        "data type": "error",
#        "resample size": 2
#    }
#    urutau.add_module(SpatialResampler, stat_spat_resample_cfg)

    # First spectral resampler for the data hdu
    data_spec_resample_cfg = {
        "hdu target": "SCI",
        "data type": "flux",
        "resample size": 1.
    }
    urutau.add_module(SpectralResampler, data_spec_resample_cfg)

    # Second spectral resampler for the stat hdu
    stat_spec_resample_cfg = {
        "hdu target": "ERR",
        "data type": "error",
        "resample size": 1.
    }
    urutau.add_module(SpectralResampler, stat_spec_resample_cfg)

#    # Degrade resolving power for the data hdu (not used, if you whant to use it change the HDU names acordingly on the next steps)
#    data_degrade_cfg = {
#        "hdu target": "DATA_RSP_BIN",
#        "data type": "flux",
#        "r input": 3027,
#        "r output": 2000
#    }
#    urutau.add_module(DegradeData, data_degrade_cfg)

#    # Degrade resolving power for the stat hdu (not used, if you whant to use it change the HDU names acordingly on the next steps)
#    data_degrade_cfg = {
#        "hdu target": "STAT_RSP_BIN",
#        "data type": "variance",
#        "r input": 3027,
#        "r output": 2000
#    }
#    urutau.add_module(DegradeData, data_degrade_cfg)

    # Signal To Noise Mask With Error
    sn_mask_cfg = {
        "hdu flux": "SCI_BIN",
        "hdu error": "ERR_BIN",
        "sn window": [21910.00, 21966.00],
        "thresholds": [5, 10]
    }
    urutau.add_module(SNMaskWithError, sn_mask_cfg)

    # Run Starlight
    population_ages = {
        "xyy": (2., 1E7),
        "xyo": (1E7, 5.6E7),
        "xiy": (5.6E7, 5E8),
        "xii": (5E8, 8E8),
        "xio": (8E8, 2E9),
        "xo": (2E9, 13E9)
    }
    
    # define SFRs (see Riffel+2021 for details)
    sfr_age_par = {
    "SFR_100": (2.,1.E8),
    "SFR_200": (2.,2.E8)
    }
    # define AGN parameters (see RIffel+09 and Riffel+22)
    fc_par = {
    "FC_25":  (0.25,0.25),
    "FC_50":  (0.5,0.5),
    "FC_75":  (0.75,0.75),
    "FC_tot": (0.245,0.75)
    }
    
    bb_par = {
    "BB_Cool": (699,1000),
    "BB_hot": (1000,1400),
    "BB_Tot": (699,1400)
    }
    
   # you have first to install starlight and configure configure the grid_example accordingly.
    starlight_cfg = {
        "starlight path": "/home/riffel/WorkOn/nirspec_JWST/starlight/StarlightChains_v04RR_25klines_1000Base.exe",
        "default grid file": "/home/riffel/WorkOn/nirspec_JWST/starlight/grid_example.inp",
        "hdu flux": "SCI_BIN",
        "hdu error": "ERR_BIN",
        "hdu flag": "SN_MASKS_5",
        "number of threads": 48 ,
        "population ages": population_ages,
        "sfr ages" : sfr_age_par,
        "fc exps": fc_par,
        "bb temps": bb_par,
        "keep tmp": True,

    }
    urutau.add_module(StarlightOnUrutau, starlight_cfg)

    # Load targets
    urutau.read_csv(targets_dir="./cubes/",
                    csv_file="./cubes/cubes.csv")

    # Execute urutau
    urutau.execute("./runs/", save_config=True, debug=True)

if __name__=="__main__":
    nirspec_JWST()