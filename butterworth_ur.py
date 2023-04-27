import numpy as np
from urutau import UrutauModule
import astropy.io.fits as fits


class ButterworthFilter(UrutauModule):
    """
        Module to apply a 2D butterworth along the z axis of a datacube HDU.
    """

    name = "Butterworth Filter"

    def _set_init_default_parameters(self) -> None:
        self.default_parameters["hdu extension"] = "FLUX"
        self.default_parameters["order"] = 3
        self.default_parameters["range"] = 0.3

    def execute(self, target_file: str) -> fits.FitsHDU:
        hdu_data = fits.getdata(target_file, self["hdu extension"])
        hdu_header = fits.getheader(target_file, self["hdu extension"])

        bw_data = self._butterworth_filter(flux_data=hdu_data,
                                           bw_order=self._config["order"],
                                           bw_range=self._config["range"])

        hdu = fits.ImageHDU(data=bw_data, header=hdu_header)
        hdu.header["EXTNAME"] = "FLUX_BW"
        hdu.header["ORDER"] = self["order"]
        hdu.header["RANGE"] = self["range"]

        return hdu

    def _butterworth_filter(self, flux_data: np.ndarray, bw_order: float, bw_range: float) -> np.ndarray:
        data_dimension = len(flux_data.shape)

        if data_dimension == 2:
            y_size, x_size = flux_data.shape
            z_size = 1
        elif data_dimension == 3:
            z_size, y_size, x_size = flux_data.shape
        else:
            raise IndexError("Dimension of hdu data must be 2 or 3.")

        radial_matrix = self._get_radial_matrix(x_size, y_size)

        h_matrix = np.zeros_like(radial_matrix)
        h_matrix = 1. - 1./(1. + (bw_range/radial_matrix)**(2.*bw_order))
        aux = np.roll(h_matrix, - int(y_size/2), axis=0)
        h_matrix = np.roll(aux, - int(x_size/2), axis=1)

        results = np.zeros_like(flux_data)

        for k in range(0, z_size):
            temp = np.fft.fft2(flux_data[k, :, :])
            temp = temp*h_matrix
            results[k, :, :] = np.fft.ifft2(temp).real

        return results

    def _get_radial_matrix(self, x_size: int, y_size: int) -> np.ndarray:

        radial_matrix = np.zeros((x_size, y_size))

        x_correction = 0.5 if x_size % 2 == 0 else 0.0
        y_correction = 0.5 if y_size % 2 == 0 else 0.0

        for i in range(0, x_size):
            for j in range(0, y_size):
                x_element = (i - x_size/2 + x_correction) ** 2.
                y_element = (j - y_size/2 + y_correction) ** 2.
                radial_matrix[i, j] = np.sqrt(x_element + y_element)

        return radial_matrix/radial_matrix.max()
