"""
    Map filters for Urutau.
"""

import astropy.io.fits as fits
import numpy as np

from ._module_base import AbstractModule


class ButterworthFilter(AbstractModule):
    """
        Module to apply a 2D butterworth along the z axis of a datacube HDU.

        Default Urutau Parameters:
            - "hdu flux" = hdu name with flux data
            - "order" = butterworth order
            - "range" = butterworth range

        Resulting Extension Name = "FLUX_BW"
    """

    name = "Butterworth Filter"

    def _set_init_default_parameters(self) -> None:
        self.default_parameters["hdu flux"] = "FLUX"
        self.default_parameters["order"] = 3
        self.default_parameters["range"] = 0.3

    def execute(self, input_hdu: fits.HDUList) -> fits.HDUList:
        ext_name = self["hdu flux"]

        hdu_data = input_hdu[ext_name].data
        hdu_header = input_hdu[ext_name].header

        bw_data = self._filter(flux_data=hdu_data,
                               bw_order=self._config["order"],
                               bw_range=self._config["range"])

        hdu = fits.ImageHDU(data=bw_data, header=hdu_header)
        hdu.header["EXTNAME"] = "FLUX_BW"
        hdu.header["BW_ORDER"] = self["order"]
        hdu.header["BW_RANGE"] = self["range"]

        return fits.HDUList(hdus=[hdu])

    def _filter(self, flux_data: np.ndarray, bw_order: float, bw_range: float) -> np.ndarray:
        data_dimension = len(flux_data.shape)

        if data_dimension == 2:
            y_size, x_size = flux_data.shape
        elif data_dimension == 3:
            _, y_size, x_size = flux_data.shape
        else:
            raise IndexError("Dimension of hdu data must be 2 or 3.")

        radial_matrix = self._get_radial_matrix(y_size, x_size)

        h_matrix = np.zeros_like(radial_matrix)
        h_matrix = 1. - 1./(1. + (bw_range/radial_matrix)**(2.*bw_order))
        aux = np.roll(h_matrix, - int(y_size/2), axis=0)
        h_matrix = np.roll(aux, - int(x_size/2), axis=1)

        # Replace empty/inf values with 0 so the fourier doesn't return trash
        flux_data_clean = flux_data.copy()
        flux_data_clean[~np.isfinite(flux_data_clean)] = 0.

        temp = np.fft.fft2(flux_data_clean, axes=[2, 1])
        temp = temp*h_matrix

        results = np.fft.ifft2(temp).real

        results[~np.isfinite(flux_data)] = np.nan

        return results

    def _get_radial_matrix(self, y_size: int, x_size: int) -> np.ndarray:

        radial_matrix = np.zeros((y_size, x_size))

        x_correction = 0.5 if x_size % 2 == 0 else 0.0
        y_correction = 0.5 if y_size % 2 == 0 else 0.0

        for i in range(0, y_size):
            for j in range(0, x_size):
                x_element = (j - x_size/2 + x_correction) ** 2.
                y_element = (i - y_size/2 + y_correction) ** 2.
                radial_matrix[i, j] = np.sqrt(x_element + y_element)

        return radial_matrix/radial_matrix.max()
