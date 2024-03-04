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

### Under OS system

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
