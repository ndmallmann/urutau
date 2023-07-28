"""
    Map data resampling for Urutau.
"""

from math import floor

import astropy.io.fits as fits
import numpy as np

from ._module_base import AbstractModule


class Resampler(AbstractModule):
    """
        Abstract data resampler module for Urutau.

        Default Urutau Parameters:
            - "hdu target" = hdu extension to be resampled (default = "FLUX")
            - "data type" = data types of the hdu (default = "flux")
            - "resample size" = Number of pixels used per sample (default = 4)

        Obs:
            data types are "flux", "inv_flux", "error", "inv_error",
            "variance", and "inv_variance"

        Resulting Extension Name = "X_RSP", where X is the
            original extension name
    """

    name = "Resampler"

    def _set_init_default_parameters(self) -> None:
        self.default_parameters["hdu target"] = "FLUX"
        self.default_parameters["data type"] = "flux"
        self.default_parameters["resample size"] = 4

    def execute(self, input_hdu: fits.HDUList) -> fits.HDUList:
        hdu = self["hdu target"]
        type_data = self["data type"]

        header = input_hdu[hdu].header
        data = input_hdu[hdu].data

        new_hdu = self._resample_hdu_data(header, data, type_data)

        return fits.HDUList(hdus=[new_hdu])

    def _resample_hdu_data(self, header: fits.Header, data: np.ndarray, data_type: str) -> fits.FitsHDU:
        z_size, y_size, x_size = data.shape

        resample_size = self["resample size"]
        resample_x = floor(x_size/resample_size)
        resample_y = floor(y_size/resample_size)

        resampled_matrix = np.zeros((z_size, resample_y, resample_x))

        for i in range(resample_x):
            for j in range(resample_y):
                pos_x = i * resample_size
                pos_y = j * resample_size
                spectrum = self._sampled_spectrum(
                    data, pos_x, pos_y, data_type)
                resampled_matrix[:, j, i] = spectrum

        hdu = fits.ImageHDU(data=resampled_matrix, header=header)
        hdu.header["EXTNAME"] = f"{header['EXTNAME']}_RSP"
        hdu.header["RSP_SIZE"] = resample_size

        return hdu

    def _sampled_spectrum(self, matrix: np.ndarray, x_cent: int, y_cent: int, data_type: str) -> np.ndarray:
        """
            Generates sampled spectrum from matrix and position.
        """

        size = self["resample size"]

        matrix_slice = matrix[:, y_cent:y_cent+size, x_cent:x_cent+size]

        if "inv" in data_type:
            matrix_slice = 1. / matrix_slice

        if "variance" in data_type:
            matrix_slice = np.sqrt(matrix_slice)

        spec_result = self._sum_slice(matrix_slice)

        if "variance" in data_type:
            spec_result = spec_result ** 2.

        if "inv" in data_type:
            spec_result = 1. / spec_result

        return spec_result

    def _sum_slice(self, matrix_slice: np.ndarray) -> np.ndarray:
        axis_count = self["resample size"]
        matrix_aux = np.nanmean(matrix_slice, axis=1) * axis_count
        return np.nanmean(matrix_aux, axis=1) * axis_count
