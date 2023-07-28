"""
    Spectra resolution modifier for Urutau.
"""

import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d


class DegradeSpectrum():
    """
        Module to degrade data resolution.

        Default Urutau Parameters:
            - "hdu target" = hdu extension to be resampled (default = "FLUX")
            - "r input min" = input R value at the lowest wavelength (default = 2200)
            - "r input max" = input R value at the highest wavelength (default = None)
            - "r output min" = output R value at the lowest wavelength (default = 1000)
            - "r output max" = output R value at the highest wavelength (default = None)
            - "lambda step" = desired wavelength dispersion (default = 1.)

        Resulting Extension Name = "X_DEGR", where X is the
            original extension name
    """

    name = "Resampler"

    def _set_init_default_parameters(self) -> None:
        self.default_parameters["hdu target"] = "FLUX"
        self.default_parameters["r input min"] = 2200
        self.default_parameters["r input max"] = None
        self.default_parameters["r output min"] = 1000
        self.default_parameters["r output max"] = None
        self.default_parameters["lambda step"] = 1.

    def execute(self, input_hdu: fits.HDUList) -> fits.HDUList:
        hdu = self["hdu target"]

        header = input_hdu[hdu].header
        data = input_hdu[hdu].data

        wave = self._wave_array(header)

        new_hdu = self._degraded_hdu(wave, data)

        return fits.HDUList(hdus=[new_hdu])

    def _degraded_hdu(self, wave: np.ndarray, orig_data: np.ndarray) -> np.ndarray:
        FWHM_CONST = 2.35482

        lambda_step = self["lambda step"]

        r_in_min = self["r input min"]
        r_in_max = self["r input max"]

        r_out_min = self["r output min"]
        r_out_max = self["r output min"]

        if r_in_max is None:
            r_in = r_in_min
        else:
            r_in = np.linspace(r_in_min, r_in_max, wave.size)

        if r_out_max is None:
            r_out = r_out_min
        else:
            r_out = np.linspace(r_out_min, r_out_max, wave.size)

        sigma_in = wave / (r_in * FWHM_CONST)
        sigma_out = wave / (r_out * FWHM_CONST)

        conv_sigma = np.sqrt(sigma_out ** 2. - sigma_in ** 2.)

        NUM_SIGMA_TO_CALCULATE = 20.
        wave_min = wave - NUM_SIGMA_TO_CALCULATE * conv_sigma
        wave_max = wave + NUM_SIGMA_TO_CALCULATE * conv_sigma

        final_wave = np.arange(wave[0], wave[-1], lambda_step)
        _, y_size, x_size = orig_data.shape
        z_size = len(final_wave)

        degraded_flux = np.zeros((z_size, y_size, x_size), dtype=float)

        for i in range(x_size):
            for j in range(y_size):
                flux_conv = np.zeros_like(orig_data, dtype=float)

                for i in range(len(wave)):
                    index_cut = (wave >= wave_min[i]) * (wave <= wave_max[i])
                    psf = self._get_psf(
                        wave[i], wave[index_cut], conv_sigma[i])
                    flux_conv[i] = np.nansum(
                        psf * orig_data[index_cut]) / np.nansum(psf)

                interp_func = interp1d(wave, flux_conv)
                degraded_flux[:, j, i] = interp_func(final_wave)

        hdu = fits.ImageHDU(data=degraded_flux)
        hdu.header["CDELT3"] = lambda_step
        hdu.header["CRPIX3"] = 1
        hdu.header["CRVAL3"] = final_wave[0]

        return hdu

    def _get_psf(self, x0: float, x: np.ndarray, sigma: float) -> np.ndarray:
        exponent = - ((x - x0) ** 2.) / (2 * sigma)
        return np.exp(exponent)

    def _wave_array(self, flux_header: fits.Header) -> np.ndarray:
        z_size = flux_header["NAXIS3"]

        delta_name = "CDELT3" if "CDELT3" in flux_header else "CD3_3"

        dt_wave = flux_header[delta_name]
        c_wave_position = flux_header["CRPIX3"] - 1
        c_wave_value = flux_header["CRVAL3"]

        ini_wave = c_wave_value - c_wave_position * dt_wave

        wave_array = ini_wave + np.array([x*dt_wave for x in range(0, z_size)])

        return wave_array
