# Urutau

Urutau is a modular pipeline tool developed to work with fits files and multithreading.

[![DOI](https://zenodo.org/badge/633615260.svg)](https://zenodo.org/badge/latestdoi/633615260)

## Requirements

Urutau requires Python 3.10 and the following python libraries:

- astropy
- pandas
- scipy

They are very popular libraries used in astronomy, but, in case they are not installed in your machine, urutau's installation should include them automatically.

## Installation

### Native installation

In order to install the package, just run the following command on your terminal:

```
pip install urutau@git+https://github.com/ndmallmann/urutau.git
```

Alternatively, you could clone the repository and install it using the following commands:

```
git clone https://github.com/ndmallmann/urutau.git
pip install urutau/
```

### Using CONDA enviroments

```
conda create -n urutau python=3.10 
source activate urutau
pip install urutau@git+https://github.com/ndmallmann/urutau.git
```



## How does it Work?

Urutau was developed for users to quickly assemble a pipeline for astronomy FITS files as well as process multiple datacubes in parallel.

The user feeds the software a series of modules to be executed in sequence (pipeline) as well as a list of target datacubes as input.

Each module can be configured based on default/general parameters (such as the name of the HDU extension to be used during the computations).

Each target can be loaded with specific parameters that will be automatically fed to each module in the pipeline (such as the redshift).

### Example code

As an example, here's a bit of code that generates data based on a FLUX HDU with resampled X and Y dimension:

```
from urutau import Urutau

from urutau.modules import SpatialResampler

ur = Urutau(num_threads = 1)

resampler_config = {
    "hdu target": "DATA",
    "data type": "flux",
    "resampler size": 4
}
ur.add_module(SpatialResampler, resampler_config)

ur.read_csv(targets_dir="./cubes/",
            csv_file="./target_parameters.csv")

ur.execute("./save_folder/")
```

This snippet of code uses one single module (SpatialResampler) to resize the spatial dimensions of the "DATA" HDU from all the cubes inside the directory "./cubes/" and listed in "./target_parameters.csv".

See another example code [here](/examples/using_urutau/using_urutau.py).

### Simple example to run URUTAU with Starlight (http://www.starlight.ufsc.br/) on NIRSPEC/JWST cubes

You can see the script [here](/examples/run_starlight/run_urutau_starlight_nirspec_jwst.py) 

```
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

    # First spatial resampler for the data hdu (NOT used here. To use it change the HDUs on the next steps acordingly)
    #    data_spat_resample_cfg = {
    #        "hdu target": "SCI",
    #        "data type": "flux",
    #        "resample size": 2
    #    }
    #    urutau.add_module(SpatialResampler, data_spat_resample_cfg)

    #    # Second spatial resampler for the stat hdu
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

    #    # Degrade resolving power for the data hdu
    #    data_degrade_cfg = {
    #        "hdu target": "DATA_RSP_BIN",
    #        "data type": "flux",
    #        "r input": 3027,
    #        "r output": 2000
    #    }
    #    urutau.add_module(DegradeData, data_degrade_cfg)

    #    # Degrade resolving power for the stat hdu
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

    
    # Starlight Module
    
    # Stellar populations ages parameter
    population_ages = {
        "xyy": (2., 1E7),
        "xyo": (1E7, 5.6E7),
        "xiy": (5.6E7, 5E8),
        "xii": (5E8, 8E8),
        "xio": (8E8, 2E9),
        "xo": (2E9, 13E9)
    }

    # SFR from stellar populations (see Riffel+2021 for details)
    sfr_age_par = {
        "SFR_100": (2.,1.E8),
        "SFR_200": (2.,2.E8)
    }

    # AGN featureless and hot dust components (see Riffel+2009 for details)
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
    
    # Setting up all the configurations 
    starlight_cfg = {
        "starlight path": "/home/riffel/WorkOn/nirspec_JWST/starlight/StarlightChains.exe",
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
```

just save this in a script to run it or download it [here](/examples/run_starilght/run_urutau_starlight_nirspec_jwst.py) 


## Creating Modules

Urutau wouldn't be very useful without the possibility of creating modules for the pipeline.

In order to create a module, the user needs to import the base module (an abstract class) created for this software as exemplified here:

```
from astropy.io import fits
from urutau.modules import AbstractModule

class MyModule(AbstractModule):
    """
        My module doc
    """

    def _set_init_default_parameters(self) -> None:
        ...

    def execute(self, input_hdu: fits.HDUList) -> fits.HDUList:
        ...
```

as a requirement, the module must define the methods "_set_init_default_parameters" and "execute".

Here is an [example](/examples/creating_modules/creating_simple.py) of how to create a very simple module.