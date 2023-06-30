"""
    Map data resampling for Urutau.
"""

from math import floor

import astropy.io.fits as fits
import numpy as np

from ._module_base import AbstractModule


class Resampler(AbstractModule):
    """
        Abstract data resampler module for Urutal.

        Default Urutau Parameters:
            - "hdu list" = List of hdu extensions to be resampled (3D matrix)
            - "data types list" = List of data types of each hdu (respectively)
            - "resample size" = Number of pixels used per sample

        Obs:
            data types list must be of same size as hdu list

            data types are "flux", "inv_flux", "error", "inv_error", "variance"
            "inv_variance"

        Resulting Extension Name = "X_RSP", where X is the
            original extension name
    """

    name = "Resampler"

    def _set_init_default_parameters(self) -> None:
        self.default_parameters["hdu list"] = ["FLUX"]
        self.default_parameters["data types list"] = ["flux"]
        self.default_parameters["resample size"] = 4

    def execute(self, input_hdu: fits.HDUList) -> fits.HDUList:
        hdu_list = self["hdu list"]
        type_list = self["data types list"]

        new_hdus = fits.HDUList()

        for index, res_hdu in enumerate(hdu_list):
            res_type = type_list[index]

            header = input_hdu[res_hdu].header
            data = input_hdu[res_hdu].data

            new_hdu = self._resample_hdu_data(header, data, res_type)
            new_hdus.append(new_hdu)

        return new_hdus

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
