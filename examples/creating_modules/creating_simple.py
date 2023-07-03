"""
    Examples of how to create a module for Urutau.

    The function of this module is to multiply the data of an hdu by a number.
    It is not very useful, but it illustrates how to create one.
"""

# AbstractModule
# This class is the basis for any module.
from urutau import AbstractModule

# Astropy
# Urutau works with fits files and class structures from astropy
import astropy.io.fits as fits

class SimpleMultiplier(AbstractModule):
    """
        Simple module used to multiply hdu data by any scalar.

        Default Parameters:
            - "hdu" = HDU name/index (default = "FLUX")
            - "multiply" = value used to multiply the data (default = 2)

        Resulting HDU: "FLUX_MULT"
    """


    # This is one of the two methods required
    # It is used to set the default parameters
    def _set_init_default_parameters(self) -> None:
        self.default_parameters["hdu"] = "FLUX"
        self.default_parameters["multiply"] = 2.

    # This is one of the two methods required
    # It is the main body of the module
    # It must return a list of hdus from astropy
    def execute(self, input_hdu: fits.HDUList) -> fits.HDUList:
        # Each parameter can be loaded using self["parameter_name"]
        hdu_name = self["hdu"]
        multiplier = self["multiply"]

        # Generating a simple hdu
        cards_list = [fits.Card("EXTNAME", "FLUX_MULT"),
                      fits.Card("MULT", multiplier,  "multiplication value")]
        new_header = fits.Header(cards = cards_list)
        new_data = input_hdu[hdu_name].data * multiplier
        
        new_hdu = fits.ImageHDU(data=new_data, header=new_header)

        # Returning a list of hdus
        list_hdus = fits.HDUList(hdus=[new_hdu])
        return list_hdus