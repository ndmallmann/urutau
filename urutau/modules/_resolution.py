"""
    Spectra resolution modifier for Urutau.
"""

import numpy as np
from astropy.io import fits

from ._module_base import AbstractModule


class DegradeData(AbstractModule):
    """
        Module to degrade data resolution.

        Default Urutau Parameters:
            - "hdu target" = hdu extension to be resampled (default = "FLUX")
            - "data type" = data types of the hdu (default = "flux")
            - "r input" = input value/array with R values (default = 2200)
            - "r output" = output value/array with R values (default = 1000)

        Obs: r input and r output can be floats (for constant R) or arrays with
        size NAXIS3 (length of z axis on the HDU)

        Resulting Extension Name = "X_DEGR", where X is the
            original extension name

        Datacubes must contain HDU with EXTNAME = flux_hdu (str or int) and the
        following header parameters:
            - "CD3_3" or "CDELT3"  =  DELTA LAMBDA
            - "CRPIX3" =  ARRAY POSITION OF CENTRAL WAVELENGTH
            - "CRVAL3" =  CENTRAL WAVELENGTH VALUE
    """

    FWHM_CONST = 2.35482

    name = "Degrade Data"

    def _set_init_default_parameters(self) -> None:
        self.default_parameters["hdu target"] = "FLUX"
        self.default_parameters["data type"] = "flux"
        self.default_parameters["r input"] = 2200
        self.default_parameters["r output"] = 1000

    def execute(self, input_hdu: fits.HDUList) -> fits.HDUList:
        hdu = self["hdu target"]
        data_type = self["data type"]

        header = input_hdu[hdu].header
        data = input_hdu[hdu].data

        wave = self._wave_array(header)

        new_hdu = self._degraded_hdu(wave, data, data_type)
        orig_header_name = header["EXTNAME"]

        new_hdu.header["EXTNAME"] = f"{orig_header_name}_DEGR"
        delta_name = "CDELT3" if "CDELT3" in header else "CD3_3"
        new_hdu.header[delta_name] = header[delta_name]
        new_hdu.header["CRPIX3"] = header["CRPIX3"]
        new_hdu.header["CRVAL3"] = header["CRVAL3"]

        return fits.HDUList(hdus=[new_hdu])

    def _compute_conv_sigma(self, wave: np.ndarray) -> np.ndarray:
        r_in = self["r input"]
        r_out = self["r output"]

        sigma_in = wave / (r_in * self.FWHM_CONST)
        sigma_out = wave / (r_out * self.FWHM_CONST)

        diff_sigma = sigma_out ** 2. - sigma_in ** 2.
        diff_sigma[diff_sigma < 0] = 0.
        return np.sqrt(diff_sigma)

    def _degraded_hdu(self, wave: np.ndarray, orig_data: np.ndarray, data_type: str) -> fits.ImageHDU:
        conv_sigma = self._compute_conv_sigma(wave)

        NUM_SIGMA_TO_CALCULATE = 20.
        wave_min = wave - NUM_SIGMA_TO_CALCULATE * conv_sigma
        wave_max = wave + NUM_SIGMA_TO_CALCULATE * conv_sigma

        aux_matrix = orig_data

        if "inv_" in data_type:
            aux_matrix = 1. / aux_matrix

        if "variance" in data_type:
            aux_matrix = np.sqrt(aux_matrix)

        z_size, y_size, x_size = aux_matrix.shape
        degraded_matrix = np.zeros(aux_matrix.shape, dtype=float)
        for i in range(z_size):
            index_cut = (wave >= wave_min[i]) * (wave <= wave_max[i])
            psf = self._get_psf(
                wave[i], wave[index_cut], conv_sigma[i])

            sum_psf = np.nansum(psf)
            if sum_psf > 0.:
                conv_map = np.zeros((psf.size, y_size, x_size), dtype=float)
                orig_map = aux_matrix[index_cut, :, :]

                for k in range(psf.size):
                    conv_map[k, :, :] = psf[k] * orig_map[k, :, :]

                degraded_matrix[i, :, :] = np.nansum(
                    conv_map, axis=0) / sum_psf

        if "variance" in data_type:
            degraded_matrix = degraded_matrix ** 2.

        if "inv_" in data_type:
            degraded_matrix = 1. / degraded_matrix

        hdu = fits.ImageHDU(data=degraded_matrix)

        return hdu

    def _get_psf(self, x0: float, x: np.ndarray, sigma: float) -> np.ndarray:
        exponent = - ((x - x0) ** 2.) / (2 * sigma ** 2.)
        return np.exp(exponent)


class DegradeDataFlex(DegradeData):
    """
        Module to degrade data resolution with full flexibility to mix and match
        resolution units for the input and output spectra.

        Supported resolution units:
            - "R"     : Resolving power  R = λ / FWHM
                        → sigma(λ) = λ / (R × 2.35482)   [varies with λ]
            - "FWHM"  : Full Width at Half Maximum in Angstrom
                        → sigma(λ) = FWHM / 2.35482       [constant in λ]
            - "sigma" : Gaussian standard deviation in Angstrom
                        → sigma(λ) = sigma                [constant in λ]

        Physics reminder:
            Constant R  → FWHM and sigma grow proportionally with λ.
            Constant sigma/FWHM → R grows proportionally with λ.

        The convolution kernel applied at each wavelength is:
            conv_sigma(λ) = sqrt( sigma_out(λ)² − sigma_in(λ)² )

        Default Urutau Parameters:
            - "hdu target"   = hdu extension to be resampled (default = "FLUX")
            - "data type"    = data type of the hdu (default = "flux")
            - "input type"   = resolution unit of the input data:
                               "R", "sigma", or "FWHM" (default = "R")
            - "input value"  = input resolution value; scalar or array of size
                               NAXIS3 (default = 2200)
            - "output type"  = resolution unit of the output data:
                               "R", "sigma", or "FWHM" (default = "R")
            - "output value" = output resolution value; scalar or array of size
                               NAXIS3 (default = 1000)

        Examples:
            R=2200 (constant R, sigma varies with λ) → R=1000 :
                {"input type": "R", "input value": 2200,
                 "output type": "R", "output value": 1000}

            R=2200 at λ_ref=5500 Å (constant sigma/FWHM) → R=1000 at λ_ref=5500 Å :
                {"input type": "R",  "input value": 2200, "input ref wave": 5500,
                 "output type": "R", "output value": 1000, "output ref wave": 5500}

            R=2200 → FWHM=2.5 Å (MILES) :
                {"input type": "R",     "input value": 2200,
                 "output type": "FWHM", "output value": 2.5}

            FWHM=2 Å → sigma=3 Å :
                {"input type": "FWHM",  "input value": 2.0,
                 "output type": "sigma", "output value": 3.0}

        Resulting Extension Name = "X_DEGR", where X is the
            original extension name

        Datacubes must contain HDU with EXTNAME = hdu_target (str or int) and
        the following header parameters:
            - "CD3_3" or "CDELT3"  =  DELTA LAMBDA
            - "CRPIX3" =  ARRAY POSITION OF CENTRAL WAVELENGTH
            - "CRVAL3" =  CENTRAL WAVELENGTH VALUE
    """

    name = "Degrade Data Flex"

    def _set_init_default_parameters(self) -> None:
        super()._set_init_default_parameters()
        self.default_parameters["input type"] = "R"
        self.default_parameters["input value"] = 2200
        self.default_parameters["input ref wave"] = None   # optional reference wavelength for "R" type
        self.default_parameters["output type"] = "R"
        self.default_parameters["output value"] = 1000
        self.default_parameters["output ref wave"] = None  # optional reference wavelength for "R" type

    def _resolution_to_sigma(self, wave: np.ndarray, value, kind: str,
                             ref_wave=None) -> np.ndarray:
        """
        Convert a resolution specification to a per-wavelength sigma array (Å).

        For "R" type:
            - Without ref_wave : sigma(λ) = λ / (R × FWHM_CONST)  [varies with λ]
            - With    ref_wave : sigma     = λ_ref / (R × FWHM_CONST)  [constant,
                                 equivalent to fixing the FWHM at the reference λ]
        """
        kind = kind.lower()
        if kind == "r":
            if ref_wave is not None:
                # R given at a reference wavelength → constant sigma across all λ
                # (same as specifying FWHM = λ_ref / R)
                sigma_ref = float(ref_wave) / (float(value) * self.FWHM_CONST)
                return np.full_like(wave, sigma_ref)
            else:
                # Constant R → sigma grows proportionally with λ
                return wave / (np.asarray(value, dtype=float) * self.FWHM_CONST)
        elif kind == "fwhm":
            fwhm = np.full_like(wave, float(value)) if np.isscalar(value) else np.asarray(value, dtype=float)
            return fwhm / self.FWHM_CONST
        elif kind == "sigma":
            if np.isscalar(value):
                return np.full_like(wave, float(value))
            return np.asarray(value, dtype=float)
        else:
            raise ValueError(
                f"Unknown resolution type '{kind}'. "
                "Valid options are: 'R', 'FWHM', 'sigma'."
            )

    def _compute_conv_sigma(self, wave: np.ndarray) -> np.ndarray:
        sigma_in = self._resolution_to_sigma(
            wave, self["input value"], self["input type"], self["input ref wave"])
        sigma_out = self._resolution_to_sigma(
            wave, self["output value"], self["output type"], self["output ref wave"])

        diff_sigma = sigma_out ** 2. - sigma_in ** 2.
        diff_sigma[diff_sigma < 0] = 0.
        return np.sqrt(diff_sigma)


    def _wave_array(self, flux_header: fits.Header) -> np.ndarray:
        z_size = flux_header["NAXIS3"]

        delta_name = "CDELT3" if "CDELT3" in flux_header else "CD3_3"

        dt_wave = flux_header[delta_name]
        c_wave_position = flux_header["CRPIX3"] - 1
        c_wave_value = flux_header["CRVAL3"]

        ini_wave = c_wave_value - c_wave_position * dt_wave

        wave_array = ini_wave + np.array([x*dt_wave for x in range(0, z_size)])

        return wave_array
