"""
    Signal to noise map generator for Urutau.
"""

import astropy.io.fits as fits
import numpy as np

from ._module_base import AbstractModule


class GenericSNMask(AbstractModule):
    """
        Module to generate signal to noise maps based on threshold values. It
        utilizes average and standard deviation values of the flux to calculate
        the signal to noise map at the desired window range.

        Default Urutau Parameters:
            - "hdu flux" = hdu name with flux data (default = "FLUX")
            - "sn window" = signal to noise window (default = [4000, 6000])
            - "thresholds" = signal to noise thresholds list (default = [10])
            - "redshift" = redshift of the object (default = 0.)

        Resulting Extension Names = "SN_MASKS_X",
            where X is the threshold value

        Obs:
            each threshold value generates a different HDU extension.

        Datacubes must contain HDU with EXTNAME = flux_hdu (str or int) and the
        following header parameters:
            - "CD3_3" or "CDELT3"  =  DELTA LAMBDA
            - "CRPIX3" =  ARRAY POSITION OF CENTRAL WAVELENGTH
            - "CRVAL3" =  CENTRAL WAVELENGTH VALUE
    """

    name = "Generic SN Mask"

    def _set_init_default_parameters(self) -> None:
        self.default_parameters["hdu flux"] = "FLUX"
        self.default_parameters["sn window"] = [4000, 6000]
        self.default_parameters["thresholds"] = [10]
        self.default_parameters["redshift"] = 0.

    def execute(self, input_hdu: fits.HDUList) -> fits.HDUList:
        ext_name = self["hdu flux"]
        flux_data = input_hdu[ext_name].data
        flux_header = input_hdu[ext_name].header

        z_size, _, _ = flux_data.shape

        wavelength = self._redshifted_wave_array(flux_header, z_size)

        left_lambda, right_lambda = self["sn window"]
        left_index = self._min_index(left_lambda, wavelength)
        right_index = self._max_index(right_lambda, wavelength)

        sn_ratio_map = self._sn_map(input_hdu, left_index, right_index)

        hdus_list = fits.HDUList()

        for threshold in self["thresholds"]:
            hdu = self._sn_mask(left_lambda, right_lambda,
                                sn_ratio_map, threshold)
            hdus_list.append(hdu)

        return hdus_list

    def _sn_map(self, input_hdu: fits.HDUList, left_index: int, right_index: int) -> np.ndarray:
        flux_data = input_hdu[self["hdu flux"]].data

        flux_at_window = flux_data[left_index:right_index, :, :]
        mean_at_window = np.nanmean(flux_at_window, axis=0)
        std_deviation_map = np.nanstd(flux_at_window, axis=0)

        sn_ratio_map = mean_at_window / std_deviation_map

        return sn_ratio_map

    def _sn_mask(self, min_l: float, max_l: float, sn_ratio: np.ndarray, limit: float) -> fits.FitsHDU:
        sn_map = np.zeros_like(sn_ratio, dtype=int)
        sn_map[sn_ratio > limit] = 1

        hdu = fits.ImageHDU(data=sn_map)
        hdu.header["EXTNAME"] = f"SN_MASKS_{int(limit)}"
        hdu.header["SN_WIND"] = f"{min_l} - {max_l}"
        hdu.header["THRESH"] = limit
        return hdu

    def _min_index(self, l_value: float, wavelength: np.ndarray) -> int:
        return np.sum(wavelength < l_value)

    def _max_index(self, l_value: float, wavelength: np.ndarray) -> int:
        return np.sum(wavelength < l_value) - 1

    def _redshifted_wave_array(self, flux_header: fits.Header, z_size: int) -> np.ndarray:
        delta_name = "CDELT3" if "CDELT3" in flux_header else "CD3_3"

        dt_wave = flux_header[delta_name]
        c_wave_position = flux_header["CRPIX3"] - 1
        c_wave_value = flux_header["CRVAL3"]

        ini_wave = c_wave_value - c_wave_position * dt_wave

        wave_array = ini_wave + np.array([x*dt_wave for x in range(0, z_size)])

        return wave_array / (1. + self["redshift"])


class SNMaskWithError(GenericSNMask):
    """
        Module to generate signal to noise maps based on threshold values. It
        utilizes the error data from an hdu to calculate the desired signal
        to noise map at window range.

        Default Urutau Parameters:
            - "hdu flux" = hdu name with flux data (default = "FLUX")
            - "hdu error" = hdu name with error data (default = "ERROR")
            - "sn window" = signal to noise window (default = [4000, 6000])
            - "thresholds" = list of signal to noise thresholds (default = [10])
            - "redshift" = redshift of the object (default = 0.)

        Resulting Extension Names = "SN_MASKS_X",
            where X is the threshold value

        Obs:
            each threshold value generates a different HDU extension.

        Datacubes must contain HDU with EXTNAME = flux_hdu (str or int) and the
        following header parameters:
            - "CD3_3" or "CDELT3"  =  DELTA LAMBDA
            - "CRPIX3" =  ARRAY POSITION OF CENTRAL WAVELENGTH
            - "CRVAL3" =  CENTRAL WAVELENGTH VALUE
    """

    name = "SN Mask With Error"

    def _set_init_default_parameters(self) -> None:
        super()._set_init_default_parameters()
        self.default_parameters["hdu error"] = "ERROR"

    def _sn_map(self, input_hdu: fits.HDUList, left_index: int, right_index: int) -> np.ndarray:
        flux_data = input_hdu[self["hdu flux"]].data
        error_data = input_hdu[self["hdu error"]].data

        sn_ratio_map = np.zeros_like(flux_data[0, :, :])
        flux_cut = flux_data[left_index:right_index, :, :]
        mean_flux = np.mean(flux_cut, axis=0)

        error_cut = error_data[left_index:right_index, :, :]
        mean_error = np.mean(error_cut, axis=0)

        good_ind = mean_error > 0
        sn_ratio_map[good_ind] = mean_flux[good_ind] / mean_error[good_ind]

        return sn_ratio_map


class SNMaskWithIVar(GenericSNMask):
    """
        Module to generate signal to noise maps based on threshold values. It
        utilizes the inverse variance data from an hdu to calculate the desired
        signal to noise map at window range.

        Default Urutau Parameters:
            - "hdu flux" = hdu name with flux data (default = "FLUX")
            - "hdu ivar" = hdu name with inverse variance data (default = "IVAR")
            - "sn window" = signal to noise window (default = [4000, 6000])
            - "thresholds" = list of signal to noise thresholds (default = [10])
            - "redshift" = redshift of the object (default = 0.)

        Resulting Extension Names = "SN_MASKS_X",
            where X is the threshold value

        Obs:
            each threshold value generates a different HDU extension.

        Datacubes must contain HDU with EXTNAME = flux_hdu (str or int) and the
        following header parameters:
            - "CD3_3" or "CDELT3"  =  DELTA LAMBDA
            - "CRPIX3" =  ARRAY POSITION OF CENTRAL WAVELENGTH
            - "CRVAL3" =  CENTRAL WAVELENGTH VALUE
    """

    name = "SN Mask With Inverse Variance"

    def _set_init_default_parameters(self) -> None:
        super()._set_init_default_parameters()
        self.default_parameters["hdu ivar"] = "IVAR"

    def _sn_map(self, input_hdu: fits.HDUList, left_index: int, right_index: int) -> np.ndarray:
        flux_data = input_hdu[self["hdu flux"]].data
        ivar_data = input_hdu[self["hdu ivar"]].data

        flux_cut = flux_data[left_index:right_index, :, :]
        mean_flux = np.mean(flux_cut, axis=0)

        sqrt_ivar_data = np.sqrt(ivar_data)
        sqrt_ivar_cut = sqrt_ivar_data[left_index:right_index, :, :]
        mean_sqrt_ivar = np.mean(sqrt_ivar_cut, axis=0)

        sn_ratio_map = mean_flux * np.sqrt(mean_sqrt_ivar)

        return sn_ratio_map


class SNMaskWithVar(GenericSNMask):
    """
        Module to generate signal to noise maps based on threshold values. It
        utilizes the variance data from an hdu to calculate the desired signal
        to noise map at window range.

        Default Urutau Parameters:
            - "hdu flux" = hdu name with flux data (default = "FLUX")
            - "hdu var" = hdu name with inverse variance data (default = "VAR")
            - "sn window" = signal to noise window (default = [4000, 6000])
            - "thresholds" = list of signal to noise thresholds (default = [10])
            - "redshift" = redshift of the object (default = 0.)

        Resulting Extension Names = "SN_MASKS_X",
            where X is the threshold value

        Obs:
            each threshold value generates a different HDU extension.

        Datacubes must contain HDU with EXTNAME = flux_hdu (str or int) and the
        following header parameters:
            - "CD3_3" or "CDELT3"  =  DELTA LAMBDA
            - "CRPIX3" =  ARRAY POSITION OF CENTRAL WAVELENGTH
            - "CRVAL3" =  CENTRAL WAVELENGTH VALUE
    """

    name = "SN Mask With Variance"

    def _set_init_default_parameters(self) -> None:
        super()._set_init_default_parameters()
        self.default_parameters["hdu var"] = "VAR"

    def _sn_map(self, input_hdu: fits.HDUList, left_index: int, right_index: int) -> np.ndarray:
        flux_data = input_hdu[self["hdu flux"]].data
        var_data = input_hdu[self["hdu var"]].data

        flux_cut = flux_data[left_index:right_index, :, :]
        mean_flux = np.mean(flux_cut, axis=0)

        sqrt_var_data = np.sqrt(var_data)
        sqrt_var_cut = sqrt_var_data[left_index:right_index, :, :]
        mean_sqrt_var = np.mean(sqrt_var_cut, axis=0)

        good_ind = mean_sqrt_var > 0
        sn_ratio_map = np.zeros_like(flux_data[0, :, :])
        sn_ratio_map[good_ind] = mean_flux[good_ind] / mean_sqrt_var[good_ind]

        return sn_ratio_map
