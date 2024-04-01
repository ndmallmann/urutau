"""
    Starlight Urutau Module
"""

from astropy.io import fits

from ..starlight_utils import GridReader, StarlightGeneric

from ._module_base import AbstractModule


class StarlightOnUrutau(AbstractModule):
    """
        This module requires the stellar population fitting code
        STARLIGHT. For more details on STARLIGHT and installation instructions visit:
        http://www.starlight.ufsc.br/ 
        
        Starlight module for urutau.

        Default Urutau Parameters:
            - "number of threads" = number of starlight executables to run in parallel (default = 1)
            - "starlight path" = path for starlight's executable (default = "")
            - "default grid file" = reference grid file (default = "")
            - "hdu flux" = hdu name with flux data (default = "FLUX")
            - "galaxy distance" = galaxy distance in Mpc (default = 0.)
            - "population ages" = population age limits (default = {"x": [0, 13E9]})
            - "hdu ivar" = [Optional] hdu name with inverse variance data (default = None)
            - "hdu error" = [Optional] hdu name with error data (default = None)
            - "hdu flag" = [Optional] hdu name with mask data (default = None)
            - "sfr ages" = [Optional] star formation rate age limits (default = {})
            - "fc exps" = [Optional] featureless continuum exponent limits (default = {})
            - "bb temps" = [Optional] black body temperature limits (default = {})
            - "redshift" = [Optional] redshift (for correction) (default = 0.)
            - "normalization factor" = [Optional] flux's normalization value (default = 1.)
            - "flux unit" = [Optional] flux unit's unit (default = "")
            - "keep tmp" = [Optional] keep temp files (True) or not (False) (default = False)

        Resulting Extension Names = "BaseAgeMetal", "POPBINS", "PopVecsL",
            "PopVecsM", "FLXOBS", "FLXSYN", "WEIGHT"

        Obs:
            population ages is a dictionary containing tuples with min and max
            age limits, the keys are the population names. Ex:

            {
                "xy": [0, 1E3],
                "xi": [1E3, 10E6]
            }

            flag_hdu can be a 2D matrix (per spaxel mask) or 3D matrix (per 
            wavelength mask)

            starlight ignores flag if error is not available

        Datacubes must contain HDU with EXTNAME = flux_hdu (str or int) and the
        following header parameters:
            - "CD3_3" or "CDELT3"  =  DELTA LAMBDA
            - "CRPIX3" =  ARRAY POSITION OF CENTRAL WAVELENGTH
            - "CRVAL3" =  CENTRAL WAVELENGTH VALUE
    """

    name = "Starlight Wrapper"

    def _set_init_default_parameters(self) -> None:

        self.default_parameters["number of threads"] = 1
        self.default_parameters["starlight path"] = ""
        self.default_parameters["hdu flux"] = "FLUX"
        self.default_parameters["hdu ivar"] = None
        self.default_parameters["hdu var"] = None
        self.default_parameters["hdu error"] = None
        self.default_parameters["hdu flag"] = None

        self.default_parameters["population ages"] = {"x": [0, 13E9]}
        self.default_parameters["sfr ages"] = {}
        self.default_parameters["fc exps"] = {}
        self.default_parameters["bb temps"] = {}
        self.default_parameters["galaxy distance"] = 0.
        self.default_parameters["redshift"] = 0.
        self.default_parameters["normalization factor"] = 1.
        self.default_parameters["flux unit"] = ""

        self.default_parameters["default grid file"] = ""

        self.default_parameters["keep tmp"] = False

    def execute(self, input_hdu: fits.HDUList) -> fits.HDUList:

        num_threads = self["number of threads"]
        sl_path = self["starlight path"]
        flux = self["hdu flux"]
        ivar = self["hdu ivar"]
        var = self["hdu var"]
        error = self["hdu error"]
        flag = self["hdu flag"]

        extra_par = dict()
        if ivar is not None:
            extra_par["ivar_hdu"] = ivar
        if var is not None:
            extra_par["var_hdu"] = var
        if error is not None:
            extra_par["error_hdu"] = error
        if flag is not None:
            extra_par["flag_hdu"] = flag

        wrapper = StarlightGeneric(sl_path, flux, num_threads, **extra_par)

        pop_par = self["population ages"]
        sfr_par = self["sfr ages"]
        fc_par = self["fc exps"]
        bb_par = self["bb temps"]

        gal_dist = self["galaxy distance"]
        redshift = self["redshift"]
        norm_factor = self["normalization factor"]
        flux_unit = self["flux unit"]

        grid_reader = GridReader(self["default grid file"])

        first_entry = grid_reader.input_entries[0]
        grid = grid_reader.parameters

        grid.starlight_config = first_entry["config"]
        grid.master_base_file = first_entry["base"]
        grid.mask_file = first_entry["mask"]
        grid.reddening_law = first_entry["reddening"]
        grid.vel_recession = first_entry["v0"]
        grid.vel_dispersion = first_entry["vd"]

        synt_hdus = wrapper.run_starlight(input_hdu, grid, pop_par, sfr_par, fc_par,
                                          bb_par, gal_dist, norm_factor, flux_unit,
                                          redshift, keep_tmp=self["keep tmp"])

        return synt_hdus
