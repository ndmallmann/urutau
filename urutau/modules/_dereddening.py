"""
    Spectra dereddening for Urutau.
"""

from abc import ABC, abstractmethod

import astropy.io.fits as fits
import numpy as np
from scipy import interpolate

from ._module_base import AbstractModule


class GenericDereddening(AbstractModule, ABC):
    """
        Generic module to apply a dereddening law.

        Default Urutau Parameters:
            - "hdu flux" = Name of the hdu extension with flux data (3D matrix)
            - "ebv" = E(B-V) value for the target galaxy's angular position
            - "rv" = Galactic Rv parameter

        Resulting Extension Name = "FLUX_DRD"

        Datacubes must contain HDU with EXTNAME = flux_hdu (str or int) and the
        following header parameters:
            - "CD3_3" or "CDELT3"  =  DELTA LAMBDA
            - "CRPIX3" =  ARRAY POSITION OF CENTRAL WAVELENGTH
            - "CRVAL3" =  CENTRAL WAVELENGTH VALUE
    """

    name = "Generic Dereddening"

    def _set_init_default_parameters(self) -> None:
        self.default_parameters["hdu flux"] = "FLUX"
        self.default_parameters["ebv"] = 0.
        self.default_parameters["rv"] = 3.1

    def execute(self, input_hdu: fits.HDUList) -> fits.HDUList:
        hdu_name = self["hdu flux"]

        flux_header = input_hdu[hdu_name].header
        flux_data = input_hdu[hdu_name].data

        z_size, _, _ = flux_data.shape

        wavelength = self._wave_array(flux_header, z_size)
        exp_func = np.array([self._law(x) for x in wavelength])
        correction = 10. ** (0.4 * exp_func)

        deredded_flux = np.zeros_like(flux_data)
        for i in range(z_size):
            deredded_flux[i, :, :] = correction[i] * flux_data[i, :, :]

        hdu = fits.ImageHDU(data=deredded_flux, header=flux_header)
        hdu.header["EXTNAME"] = "FLUX_DRD"
        hdu.header["DRD_EBV"] = self["ebv"]
        hdu.header["DRD_RV"] = self["rv"]

        return fits.HDUList(hdus=[hdu])

    def _wave_array(self, flux_header: fits.Header, z_size: int) -> np.ndarray:
        delta_name = "CDELT3" if "CDELT3" in flux_header else "CD3_3"

        dt_wave = flux_header[delta_name]
        c_wave_position = flux_header["CRPIX3"] - 1
        c_wave_value = flux_header["CRVAL3"]

        ini_wave = c_wave_value - c_wave_position * dt_wave

        wave_array = ini_wave + np.array([x*dt_wave for x in range(0, z_size)])

        return wave_array

    @abstractmethod
    def _law(self, wave_lambda: float) -> float:
        """
            Dereddening law.
        """


class CcmLaw(GenericDereddening):
    """
        Module that applies CCM's dereddening (Galactic).

        Default Urutau Parameters:
            - "hdu flux" = Name of the hdu extension with flux data (3D matrix)
            - "ebv" = E(B-V) value for the target galaxy's angular position
            - "rv" = Galactic Rv parameter

        Resulting Extension Name = "FLUX_DRD"

        Datacubes must contain HDU with EXTNAME = flux_hdu (str or int) and the
        following header parameters:
            - "CD3_3" or "CDELT3"  =  DELTA LAMBDA
            - "CRPIX3" =  ARRAY POSITION OF CENTRAL WAVELENGTH
            - "CRVAL3" =  CENTRAL WAVELENGTH VALUE
    """

    name = "CCM Law Dereddening"

    def _law(self, wave_lambda: float) -> float:

        inv_microns = 10000. / wave_lambda

        if (inv_microns < 0.3) or (inv_microns > 10.):
            return np.nan

        if (inv_microns < 1.1):
            y_val = inv_microns ** 1.61
            a_val = 0.574 * y_val
            b_val = -0.527 * y_val
        elif (inv_microns < 3.3):
            y_val = inv_microns - 1.82

            a_aux = (0.01979 + y_val * (-0.77530 + y_val * 0.32999))
            a_aux = (-0.02427 + y_val * (0.72085 + y_val * a_aux))
            a_val = 1 + y_val * (0.17699 + y_val * (-0.50447 + y_val * a_aux))

            b_aux = (-0.62251 + y_val * (5.30260 + y_val * -2.09002))
            b_aux = (1.07233 + y_val * (-5.38434 + y_val * b_aux))
            b_val = y_val * (1.41338 + y_val * (2.28305 + y_val * b_aux))

        elif (inv_microns < 5.9):
            y_val = (inv_microns - 4.67) ** 2
            a_val = 1.752 - 0.316 * inv_microns - 0.104 / (y_val + 0.341)
            b_val = -3.090 + 1.825 * inv_microns + 1.206 / (y_val + 0.263)

        elif (inv_microns < 8.0):
            y_val = (inv_microns - 4.67) ** 2
            a_val = 1.752 - 0.316 * inv_microns - 0.104 / (y_val + 0.341)
            b_val = -3.090 + 1.825 * inv_microns + 1.206 / (y_val + 0.263)

            y_val = inv_microns - 5.9
            a_val = a_val - 0.04473 * y_val**2 - 0.009779 * y_val**3
            b_val = b_val + 0.2130 * y_val**2 + 0.1207 * y_val**3

        elif (inv_microns <= 10.0):
            y_val = inv_microns - 8
            a_val = -1.072 - 0.628 * y_val + 0.137 * y_val**2 - 0.070 * y_val**3
            b_val = 13.670 + 4.257 * y_val - 0.420 * y_val**2 + 0.374 * y_val**3

        y_val = a_val * self["ebv"] + b_val
        return y_val


class SeatonLaw(GenericDereddening):
    """
        Module that applies Seaton's dereddening (Galactic).

        Default Urutau Parameters:
            - "hdu flux" = Name of the hdu extension with flux data (3D matrix)
            - "ebv" = E(B-V) value for the target galaxy's angular position
            - "rv" = Galactic Rv parameter

        Resulting Extension Name = "FLUX_DRD"

        Datacubes must contain HDU with EXTNAME = flux_hdu (str or int) and the
        following header parameters:
            - "CD3_3" or "CDELT3"  =  DELTA LAMBDA
            - "CRPIX3" =  ARRAY POSITION OF CENTRAL WAVELENGTH
            - "CRVAL3" =  CENTRAL WAVELENGTH VALUE
    """

    name = "Seaton Law Dereddening"

    def _law(self, wave_lambda: float) -> float:

        inv_microns = 10000.0/wave_lambda

        func_x_val = np.arange(1., 2.71, 0.1)
        func_y_val = np.array([1.36, 1.44, 1.84, 2.04, 2.24, 2.44,
                               2.66, 2.88, 3.14, 3.36, 3.56, 3.77,
                               3.96, 4.15, 4.26, 4.40, 4.52, 4.64])

        if (inv_microns < 1.) or (inv_microns > 2.7):
            return np.nan

        seaton_func = interpolate.interp1d(func_x_val, func_y_val)

        cor = seaton_func(inv_microns) * self["ebv"]
        return cor
