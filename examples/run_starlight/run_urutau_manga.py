"""
    Urutau for MANGA
"""

from urutau import Urutau

from urutau.modules import (
    SpatialResampler,
    SpectralResampler,
    ButterworthFilter,
    DegradeData,
    SNMaskWithVar,
    SNMaskWithIVar,
    StarlightOnUrutau,
)



def quick_manga():
    
    # Start urutau with 3 threads
    urutau = Urutau(num_threads = 1)

    # First BW filter for the data hdu
    data_BW_cfg = {
        "hdu flux": "FLUX",
        "order": 3,
        "range": 0.3
    }
    urutau.add_module(ButterworthFilter, data_BW_cfg)


    # Run Starlight
    population_ages = {
        "xyy": (0., 1E7),
        "xyo": (1E7, 5.6E7),
        "xiy": (5.6E7, 5E8),
        "xii": (5E8, 8E8),
        "xio": (8E8, 2E9),
        "xo": (2E9, 13E9)
    }
    
    # Just an example
    sfr_age_par = {

    "SFR_1E6": (2,1.001E6), 
	"SFR_5E6": (2,5.621E6), 
	"SFR_10E6": (2,10.001E6),
	"SFR_14E6":(2,14.1001E6),
	"SFR_20E6": (2,20.001E6),
	"SFR_30E6": (2,31.6001E6),
	"SFR_56E6": (2,56.201E6),
	"SFR_100E6":(2,100.001E6),
	"SFR_200E6": (2,200.001E6)
    }



    # AGN featureless just an example
    fc_par = {
    "FC_150": (1.49,1.51)
    }
    
    
    
    starlight_cfg = {
        "starlight path": "./starlight/StarlightChains_v04.amd64_g77-3.4.6-r1_static.exe",
        "default grid file": "./starlight/reference_grid_manga.in",
        "hdu flux": "FLUX_BW",
        "hdu ivar": "IVAR",
        "hdu flag": "SN_MASKS_1",
        "number of threads": 16 ,
        "population ages": population_ages,
		"sfr ages" : sfr_age_par,
        "fc exps": fc_par,
		"keep tmp": True

    }
    urutau.add_module(StarlightOnUrutau, starlight_cfg)

    # Load targets
    urutau.read_csv(targets_dir="./mangaTest/",
                    csv_file="./mangaTest/targets.csv")

    # Execute urutau
    urutau.execute("./mangaTest/out/", save_config=True, debug=True)

if __name__=="__main__":
    quick_manga()
