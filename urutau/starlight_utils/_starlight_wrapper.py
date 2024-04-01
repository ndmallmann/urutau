"""
    A collection of classes to make it easier to use Starlight
    (from Cid Fernandes) with datacubes.
"""

import itertools as it
import math
import os
import re
import subprocess as sp
import threading as th
import uuid
from abc import ABC, abstractmethod
from typing import Any, Callable

import astropy.io.fits as fits
import numpy as np

from ._starlight_grid import GridGenerator, GridParameters
from ._starlight_reader import StarlightOutput, StarlightOutputReader


class StarlightWrapper(ABC):
    """
        Wrapper to execute starlight.
    """

    def __init__(self, sl_exec: str, num_threads: int = 1) -> None:

        # Starlight
        self._sl_exec_name = os.path.join(".", os.path.basename(sl_exec))
        self._sl_dir = os.path.join(os.path.dirname(sl_exec), "")

        # Execution parameters
        self._grid_generator = GridGenerator(GridParameters(), "")
        self._extracted_files = list()
        self._sl_output_files = list()
        self._num_threads = num_threads

        # Parameters specific to a cube when running starlight
        self._pop_age = dict()
        self._name_prefix = ""
        self._redshift = 0.
        self._galaxy_distance = 0.
        self._norm_factor = 1.
        self._flux_unit = ""

    @abstractmethod
    def _extract_cube(self, cube_data: fits.HDUList, out_dir: str) -> np.ndarray:
        """
            Extract spectra from datacube.

            Parameters:
                - cube_data =  datacube HDUList
                - out_dir   =  directory to extract the data

            Returns matrix with extracted spectra paths.

            Obs:
                    If a spaxel has no extracted spectrum, the value of the
                matrix element is the empty string "" or None.
        """

    @abstractmethod
    def _mount_data(self, sl_output_matrix: np.ndarray) -> fits.HDUList:
        """
            Compile starlight results into a datacube.

            Parameters:
                - sl_output_matrix = matrix containing starlight outputs

            Obs:
                    The output matrix contains empty parameters in points where
                the data wasn't extracted. It is also in the same shape (y, x)
                as the original data.

            Returns list of fits extensions.
        """

    def run_starlight(self, cube_data: fits.HDUList, grid_parameters: GridParameters, pop_age_par: dict, sfr_age_par: dict, fc_par: dict, bb_par: dict, galaxy_distance: float, norm_factor: float, flux_unit: str, redshift: float, keep_tmp: bool = False) -> fits.HDUList:
        """
            Run starlight for the listed spectra.

            Parameters:
                - cube_data        =  datacube path
                - grid_parameters  =  grid parameter object
                - pop_age_par      =  dictionary with population age ranges, ex: {"name_pop": [age_ini_exclusive, age_fin_inclusive]}
                - sfr_age_par      =  dictionary with star formation rate age ranges, ex: {"name_sfr": [age_ini_exclusive, age_fin_inclusive]}
                - fc_par           =  dictionary with featureless continuum exponent ranges, ex: {"name_fc": [exp_ini_exclusive, exp_fin_inclusive]}
                - bb_par           =  dictionary with black body temperature (kelvin) ranges, ex: {"name_bb": [temp_ini_exclusive, temp_fin_inclusive]}
                - galaxy_distance  =  galaxy distance in Mpc
                - norm_factor      =  flux normalization factor
                - flux_unit        =  flux unit name
                - redshift         =  galaxy redshift
                - keep_tmp         =  True/False to keep temporary files
                                      generated

            Return HDUList with created HDUs
        """

        self._clear_tmp_data()

        # Set specific cube parameters
        self._name_prefix = str(uuid.uuid4()).split("-", maxsplit=1)[0]
        self._pop_age = pop_age_par
        self._sfr_age = sfr_age_par
        self._fc_par = fc_par
        self._bb_par = bb_par
        self._redshift = redshift
        self._galaxy_distance = galaxy_distance
        self._norm_factor = norm_factor
        self._flux_unit = flux_unit

        # Create extraction and output dirs
        obs_dir = grid_parameters.obs_dir
        output_dir = grid_parameters.out_dir
        self._create_dir(obs_dir)
        self._create_dir(output_dir)

        # Setup grid generator
        self._grid_generator.clear_grids()
        self._grid_generator.change_grid_parameters(grid_parameters)
        self._grid_generator.change_grid_prefix(self._name_prefix)
        self._grid_generator.change_grid_dir(output_dir)

        # Extract data
        ext_matrix = self._extract_cube(cube_data, obs_dir)
        self._extracted_files = [x for x in ext_matrix.flat]

        # Generate grids, output path matrix, and run starlight
        res_matrix = self._create_all_grids(ext_matrix)
        self._sl_output_files = [x for x in res_matrix.flat]
        self._call_starlight(self._grid_generator.grids)

        # Read starlight output, generate param matrix and mount cube
        param_matrix = self._results_matrix(res_matrix)
        fits_hdus = self._mount_data(param_matrix)

        # Delete temporary files and clear data
        if not keep_tmp:
            self._delete_tmp_files()

        return fits_hdus

    def _create_all_grids(self, ext_matrix: np.ndarray) -> np.ndarray:

        dest_matrix = list()
        for extract in ext_matrix.flat:
            dest_file = ""

            if os.path.exists(extract):
                dest_file = self._grid_generator.create_grid(extract)[0]

            dest_matrix.append(dest_file)

        return np.array(dest_matrix, dtype=str).reshape(ext_matrix.shape)

    def _results_matrix(self, res_paths: np.ndarray) -> np.ndarray:
        output_reader = StarlightOutputReader()

        res_list = list()
        for result in res_paths.flat:
            parameters = output_reader.get_parameters(result)

            res_list.append(parameters)

        return np.array(res_list).reshape(res_paths.shape)

    def _create_dir(self, obs_dir):
        if not os.path.exists(obs_dir):
            os.makedirs(obs_dir)

    def _call_starlight(self, full_grids_list: list[str]) -> None:

        def _starlight_thread(grid_list: list[str], exec_name: str, exec_dir: str) -> None:
            for grid in grid_list:
                with open(grid, "rb") as grid_handler:
                    arguments = (exec_name,)
                    prog = sp.Popen(args=arguments, stdin=grid_handler,
                                    stdout=sp.DEVNULL, cwd=exec_dir)
                    prog.wait()

        workers: list[th.Thread] = list()

        for grid_list in np.array_split(full_grids_list, self._num_threads):
            arguments = (grid_list, self._sl_exec_name, self._sl_dir)
            worker = th.Thread(target=_starlight_thread, args=arguments)
            worker.start()
            workers.append(worker)

        for worker in workers:
            worker.join()


    def _delete_tmp_files(self) -> None:
        self._delete_extracted_files()
        self._delete_output_files()
        self._grid_generator.delete_grids()

    def _clear_tmp_data(self) -> None:
        self._extracted_files = list()
        self._sl_output_files = list()

    def _delete_output_files(self):
        for sl_out in self._sl_output_files:
            if os.path.exists(sl_out):
                os.remove(sl_out)
                continue

            base_name = os.path.basename(sl_out)
            error_path = os.path.join(self._sl_dir, base_name)

            if os.path.exists(error_path) and not os.path.isdir(error_path):
                os.remove(error_path)

        self._sl_output_files.clear()

    def _delete_extracted_files(self):
        for extracted in self._extracted_files:
            if os.path.exists(extracted):
                os.remove(extracted)
        self._extracted_files.clear()

class StarlightGeneric(StarlightWrapper):
    """
        Run starlight for generic datacubes.

        Datacubes must contain HDU with EXTNAME = flux_hdu (str or int) and the
        following header parameters:
            - "CD3_3" or "CDELT3"  =  DELTA LAMBDA
            - "CRPIX3" =  ARRAY POSITION OF CENTRAL WAVELENGTH
            - "CRVAL3" =  CENTRAL WAVELENGTH VALUE

        Optional:
            - ivar_hdu     =  Name of inverse variance hdu
            - var_hdu      =  Name of variance hdu (will override ivar)
            - error_hdu    =  Name of error hdu (will override ivar and var)
            - flag_hdu     =  Name of mask for bad spaxels (anything different
                              from 0 is a bad spaxel)

        Obs:
            error_hdu, var_hdu and ivar_hdu can only be a 3D matrix
            (per wavelength error)

            flag_hdu can be a 2D matrix (per spaxel mask) or 3D matrix (per 
            wavelength mask)

            starlight ignores flag if error is not available
    """

    def __init__(self, starlight_exec_path: str, flux_hdu: str | int, num_threads: int = 1, **kwargs) -> None:
        self._flux_loc = flux_hdu
        self._ivar_loc = self._optional_arg("ivar_hdu", **kwargs)
        self._var_loc = self._optional_arg("var_hdu", **kwargs)
        self._error_loc = self._optional_arg("error_hdu", **kwargs)
        self._flag_loc = self._optional_arg("flag_hdu", **kwargs)
        self._matrix = None
        self._card_index = 0
        self._x_size = 0
        self._z_size = 0
        self._y_size = 0
        self._new_cards = list()
        self._new_data = list()
        self._hdu_list = fits.HDUList()

        super().__init__(starlight_exec_path, num_threads)

    def _optional_arg(self, name: str, **kwargs) -> Any:
        if name in kwargs:
            return kwargs[name]
        return None

    def _extract_cube(self, cube_data: fits.HDUList, out_dir: str) -> np.ndarray:
        error_limit = 9999
        flag_limit = 100
        tolerance = 0.75

        # Read flux and optional [error/ivar/flag] data
        if self._flux_loc in cube_data:
            flux_data = cube_data[self._flux_loc].data
        else:
            return np.zeros((0, 0), dtype=str)

        has_var_val = not self._var_loc is None
        has_ivar_val = not self._ivar_loc is None
        has_error_val = not self._error_loc is None

        has_error = True
        if has_error_val and self._error_loc in cube_data:
            error_data = cube_data[self._error_loc].data
        elif has_ivar_val and self._ivar_loc in cube_data:
            ivar_data = cube_data[self._ivar_loc].data
            good_ind = ivar_data > 0

            error_data = np.full_like(ivar_data, error_limit)
            error_data[good_ind] = 1. / np.sqrt(ivar_data[good_ind])
        elif has_var_val and self._var_loc in cube_data:
            var_data = cube_data[self._var_loc].data
            good_ind = var_data > 0.

            error_data = np.full_like(var_data, error_limit)
            error_data[good_ind] = np.sqrt(var_data[good_ind])
        else:
            has_error = False

        has_3d_flag = False
        has_2d_flag = False
        if self._flag_loc in cube_data:
            flag_data = cube_data[self._flag_loc].data
            has_3d_flag = np.ndim(flag_data) == 3
            has_2d_flag = not has_3d_flag

        has_both_flags = has_error and has_3d_flag

        # Modify grid parameters to include or exclude table columns
        self._grid_generator.parameters.has_errors = 1 if has_error else 0
        self._grid_generator.parameters.has_flags = 1 if has_both_flags else 0

        # Generate wave array
        flux_header = cube_data[self._flux_loc].header
        self._z_size, self._y_size, self._x_size = flux_data.shape

        wave = self._redshifted_wave_array(flux_header, self._z_size)

        # Setup extraction file name and table format
        file_name_format = self._name_prefix + "_{0}_{1}.spec"
        file_path_format = os.path.join(out_dir, file_name_format)

        # Set table format
        fmt_text = ["%.2f", "%.3e"]
        if has_error:
            fmt_text.append("%.3e")
        if has_both_flags:
            fmt_text.append("%i")

        file_matrix = [[""] * self._x_size for x in range(self._y_size)]

        indices = it.product(range(self._y_size), range(self._x_size))
        for j, i in indices:
            # Read and verify flag data
            if has_2d_flag and not flag_data[j, i]:
                continue
            elif has_3d_flag:
                flag_spec = (flag_data[:, j, i] > 0) * flag_limit

            # Read flux data
            spec_flux = flux_data[:, j, i]
            if np.max(spec_flux) <= 0.:
                continue

            # Read error data
            if has_error:
                spec_error = error_data[:, j, i]

            # Verify bad flux data
            bad_indices = ~np.isfinite(spec_flux)
            spec_flux[bad_indices] = 0.

            # Ignore spectrum with low quantity of good points
            good_flux = spec_flux > 0.
            if has_error:
                good_error = (spec_flux / spec_error) > 3.
                good_flux = good_flux * good_error
            if has_3d_flag:
                good_flag = flag_spec < flag_limit
                good_flux = good_flux * good_flag

            good_count = np.sum(good_flux)
            good_percentage = good_count / spec_flux.size
            if good_percentage < tolerance:
                continue

            # Resample data
            wave_ini = math.ceil(wave[0])
            wave_fin = math.floor(wave[-1])
            wave_resample = np.arange(wave_ini, wave_fin, 1)

            flux_resample = np.interp(wave_resample, wave, spec_flux)

            # Format table data
            data_list = [wave_resample, flux_resample * self._norm_factor]
            if has_error:
                spec_error[bad_indices] = error_limit
                error_resample = np.interp(wave_resample, wave, spec_error)
                data_list.append(error_resample * self._norm_factor)

                if has_3d_flag:
                    flag_spec[bad_indices] = flag_limit
                    flag_resample = np.interp(wave_resample, wave, flag_spec)
                    data_list.append(flag_resample)
            file_data = np.column_stack(data_list)

            # Save table on file
            file_path = file_path_format.format(j, i)
            np.savetxt(file_path, file_data, fmt=fmt_text, delimiter=" ")
            file_matrix[j][i] = file_path

        return np.array(file_matrix, dtype=str)

    def _clear_tmp_data(self) -> None:
        self._matrix = None

        self._card_index = 0
        self._x_size = 0
        self._y_size = 0
        self._z_size = 0

        self._new_cards = list()
        self._new_data = list()

        self._hdu_list = fits.HDUList()
        super()._clear_tmp_data()

    def _mount_data(self, sl_output_matrix: np.ndarray) -> fits.HDUList:

        self._matrix = sl_output_matrix

        self._add_base_age_metal_hdu()
        self._add_popbins_hdu()
        self._add_popvecs_light_hdu()
        self._add_popvecs_mass_hdu()
        self._add_obs_flux_hdu()
        self._add_syn_flux_hdu()
        self._add_weight_hdu()
        self._add_failed_spaxels_hdu()

        return self._hdu_list.copy()

    def _add_base_age_metal_hdu(self):
        base_age_metal = self._get_ssp_base_parameters()

        name = "BaseAgeMetal"
        summary = "Ages & Metallicities in the Base"
        self._add_hdu(name, summary, list(), base_age_metal)

    def _get_ssp_base_parameters(self):
        base_name = self._grid_generator.parameters.master_base_file
        base_path = os.path.join(self._sl_dir, base_name)

        with open(base_path, mode="r", encoding="utf-8") as base_handler:
            base_data = base_handler.readlines()

        re_brackets = re.compile(r"\[.*?\]")
        re_spaces = re.compile(r"\s+")

        number_line = base_data[0]
        number_line = re.sub(re_spaces, "", number_line)
        number_line = re.sub(re_brackets, "", number_line)
        number_value = int(number_line)

        base_age_metal = list()
        for line in base_data[1:1+number_value]:
            line_data = re.split(re_spaces, line)
            age = float(line_data[1])
            metal = float(line_data[2])
            base_age_metal.append([age, metal])
        return base_age_metal

    def _add_popbins_hdu(self) -> None:
        if len(self._fc_par) > 0:
            self._add_fc_data_to_popbins_hdu()
            
        if len(self._bb_par) > 0:
            self._add_bb_data_to_popbins_hdu()

        self._add_population_data_to_popbins_hdu()
        self._add_star_formation_rate_data_to_popbins_hdu()
        self._add_other_data_to_popbins_hdu()

        name = "PopBins"
        summary = "Synthesis Parameters & Binned Population Vectors"
        self._add_hdu(name, summary, self._new_cards, self._new_data)

    def _add_popvecs_light_hdu(self) -> None:
        name = "PopVecsL"
        summary = "Population Vectors Not Binned in Light Fractions"
        hdu_data = self._array_matrix(lambda x: x.x_j)
        self._add_hdu(name, summary, list(), hdu_data)

    def _add_popvecs_mass_hdu(self) -> None:
        name = "PopVecsM"
        summary = "Population Vectors Not Binned in Mass Fractions"
        hdu_data = self._array_matrix(lambda x: x.m_cor_j)
        self._add_hdu(name, summary, list(), hdu_data)

    def _add_obs_flux_hdu(self) -> None:
        name = "FLXOBS"
        summary = "Observed Flux (STARLIGHT output) in input units"

        cards = self._wavelength_info_cards()

        unit_comment = "Unit"
        cards.append(fits.Card("BUNIT", 1, unit_comment))

        hdu_data = self._array_matrix(
            lambda x: x.f_obs * x.fobs_norm * self._norm_factor
        )
        self._add_hdu(name, summary, cards, hdu_data)

    def _add_syn_flux_hdu(self) -> None:
        name = "FLXSYN"
        summary = "Synthetic Flux (STARLIGHT output) in input units"

        cards = self._wavelength_info_cards()

        unit_comment = "Unit"
        cards.append(fits.Card("BUNIT", 1, unit_comment))

        hdu_data = self._array_matrix(
            lambda x: x.f_syn * x.fobs_norm * self._norm_factor
        )
        self._add_hdu(name, summary, cards, hdu_data)

    def _add_weight_hdu(self) -> None:
        name = "WEIGHT"
        summary = "Weight used in the fits (STARLIGHT output)"

        cards = self._wavelength_info_cards()

        hdu_data = self._array_matrix(lambda x: x.wei)
        self._add_hdu(name, summary, cards, hdu_data)


    def _add_failed_spaxels_hdu(self) -> None:
        name = "BADSPXL"
        summary = "Map of bad spaxels generated by starlight"

        hdu_data = np.ndarray((self._y_size, self._x_size), dtype=float)

        indices = it.product(range(self._y_size), range(self._x_size))
        for j, i in indices:
            hdu_data[j, i] = self._matrix[j, i].invalid_file

        self._add_hdu(name, summary, list(), hdu_data)

    def _add_hdu(self, extname: str, summary: str, cards: list[fits.Card], hdu_data: np.ndarray) -> None:
        header = fits.Header(cards=cards)
        header["EXTNAME"] = extname
        header["SUMMARY"] = summary
        header["Author_1"] = "Nicolas Dullius Mallmann"
        header["E-mail_1"] = "nicolas.mallmann@ufrgs.br"
        header["Author_2"] = "Rogerio Riffel"
        header["E-mail_2"] = "riffel@ufrgs.br"
        hdu = fits.PrimaryHDU(data=hdu_data, header=header)
        self._hdu_list.append(hdu)

    def _add_fc_data_to_popbins_hdu(self) -> None:
        for name, exponent in self._fc_par.items():
            card_n = f"{name}"
            card_comment = f"Featureless continuum for exponentes between {exponent[0]} and {exponent[1]}",
            self._add_card_and_data(card_name=card_n,
                                    card_comment=card_comment,
                                    data_matrix=self._property_matrix(
                                        lambda x: self._fc_component(x, exponent[0], exponent[1])))

    def _add_bb_data_to_popbins_hdu(self) -> None:
        for name, temp in self._bb_par.items():
            card_n = f"{name}"
            card_comment = f"Black body for temperatures between {temp[0]} and {temp[1]} Kelvin",
            self._add_card_and_data(card_name=card_n,
                                    card_comment=card_comment,
                                    data_matrix=self._property_matrix(
                                        lambda x: self._bb_component(x, temp[0], temp[1])))

    def _add_population_data_to_popbins_hdu(self) -> None:
        for name, age in self._pop_age.items():
            card_n = f"{name}_light"
            card_comment = f"Light binned pop: {age[0]:.1E} < age <= {age[1]:.1E}"
            self._add_card_and_data(
                card_name=card_n,
                card_comment=card_comment,
                data_matrix=self._property_matrix(
                    lambda x: self._pop_by_light(x, age[0], age[1])
                )
            )

        for name, age in self._pop_age.items():
            card_n = f"{name}_mass"
            card_comment = f"Mass binned pop: {age[0]:.1E} < age <= {age[1]:.1E}"
            self._add_card_and_data(
                card_name=card_n,
                card_comment=card_comment,
                data_matrix=self._property_matrix(
                    lambda x: self._pop_by_mass(x, age[0], age[1])
                )
            )

    def _add_star_formation_rate_data_to_popbins_hdu(self) -> None:

        for name, age in self._sfr_age.items():
            card_n = f"{name}"
            card_comment = f"for ages between {age[0]:.1E} and {age[1]:.1E} years"
            self._add_card_and_data(
                card_name=card_n,
                card_comment=card_comment,
                data_matrix=self._property_matrix(
                    lambda x: self._sfr_at_age(x, age[0], age[1])
                )
            )

    def _add_other_data_to_popbins_hdu(self) -> None:
        self._add_card_and_data(card_name="Av",
                                card_comment="Optical extinction",
                                data_matrix=self._property_matrix(lambda x: x.av_min))

        self._add_card_and_data(card_name="Mage_L",
                                card_comment="Mean age light weigthed",
                                data_matrix=self._property_matrix(lambda x: self._mean_by_light(x, np.log10(x.age_j))))

        self._add_card_and_data(card_name="Mage_M",
                                card_comment="Mean age mass weigthed",
                                data_matrix=self._property_matrix(lambda x: self._mean_by_mass(x, x.m_cor_j)))

        self._add_card_and_data(card_name="MZ_L",
                                card_comment="Mean metalicity light weigthed",
                                data_matrix=self._property_matrix(lambda x: self._mean_by_light(x, x.z_j)))

        self._add_card_and_data(card_name="MZ_M",
                                card_comment="Mean metalicity mass weigthed",
                                data_matrix=self._property_matrix(lambda x: self._mean_by_mass(x, x.z_j)))

        self._add_card_and_data(card_name="Mstar",
                                card_comment="Present mass in stars (Msun, starlight M* )",
                                data_matrix=self._property_matrix(self._m_cor_t))

        self._add_card_and_data(card_name="Mpcross",
                                card_comment="Mass processed in stars (Msun, starlight M*ini)",
                                data_matrix=self._property_matrix(self._m_ini_t))

        self._add_card_and_data(card_name="F_Norm",
                                card_comment="Normalization flux in input units",
                                data_matrix=self._property_matrix(lambda x: x.fobs_norm))

        self._add_card_and_data(card_name="Sigma_star",
                                card_comment="Stellar dispersion (see starlight manual)",
                                data_matrix=self._property_matrix(lambda x: x.vd_min))

        self._add_card_and_data(card_name="vrot_star",
                                card_comment="Stellar rotation (see starlight manual)",
                                data_matrix=self._property_matrix(lambda x: x.v0_min))

        self._add_card_and_data(card_name="Adev",
                                card_comment="Precentage mean deviation (see manual)",
                                data_matrix=self._property_matrix(lambda x: x.adev))

        self._add_card_and_data(card_name="ChiSqrt",
                                card_comment="ChiSqrt/Nl_eff (see manual)",
                                data_matrix=self._property_matrix(lambda x: x.fscale_chi2))

        self._add_card_and_data(card_name="SNR",
                                card_comment="SNR on normalization window",
                                data_matrix=self._property_matrix(lambda x: x.sn_in_norm_window))

    def _add_card_and_data(self, card_name: str, card_comment: str, data_matrix: np.ndarray) -> None:
        card_keyword = f"DATA{self._card_index}"
        card = fits.Card(keyword=card_keyword,
                         value=card_name,
                         comment=card_comment)

        self._new_cards.append(card)
        self._new_data.append(data_matrix)
        self._card_index += 1

    def _property_matrix(self, func: Callable[[StarlightOutput], float]) -> np.ndarray:
        matrix = np.ndarray((self._y_size, self._x_size), dtype=float)

        indices = it.product(range(self._y_size), range(self._x_size))
        for j, i in indices:
            if not self._matrix[j, i].invalid_file:
                matrix[j, i] = func(self._matrix[j, i])

        return matrix

    def _array_matrix(self, func: Callable[[StarlightOutput], np.ndarray]) -> np.ndarray:

        biggest_array_size = 0

        indices = it.product(range(self._y_size), range(self._x_size))
        for j, i in indices:
            if self._matrix[j, i].invalid_file:
                array_data = np.array([])
                array_size = 0
                continue

            array_data = func(self._matrix[j, i])
            array_size = np.size(array_data)
           
            if biggest_array_size < array_size:
                biggest_array_size = array_size

        matrix_shape = (biggest_array_size, self._y_size, self._x_size)
        matrix = np.zeros(matrix_shape, dtype=float)

        indices = it.product(range(self._y_size), range(self._x_size))
        for j, i in indices:
            if self._matrix[j, i].invalid_file:
                array_data = np.ones((biggest_array_size)) * np.nan
            else:
                array_data = np.array(func(self._matrix[j, i]))
                array_size = np.size(array_data)

                if array_size < biggest_array_size:
                    array_data = np.resize(array_data, (biggest_array_size))
                    array_data[array_size:] = np.nan

            matrix[:, j, i] = array_data

        return matrix

    def _mean_by_light(self, sl_out: StarlightOutput, attribute: np.ndarray) -> float:
        valid_age = sl_out.age_j > 0.

        sum_of_x = np.sum(sl_out.x_j)
        if sum_of_x <= 0.:
            return 0.

        light_contribution = sl_out.x_j[valid_age] / sum_of_x
        mean_by_light = np.sum(light_contribution * attribute[valid_age])
        return mean_by_light

    def _mean_by_mass(self, sl_out: StarlightOutput, attribute: np.ndarray) -> float:
        valid_age = sl_out.age_j > 0

        sum_of_x = np.sum(sl_out.m_cor_j)
        if sum_of_x <= 0.:
            return 0.

        mass_contribution = sl_out.m_cor_j[valid_age] / sum_of_x
        mean_by_mass = np.sum(mass_contribution * attribute[valid_age])
        return mean_by_mass

    def _m_ini_t(self, sl_out: StarlightOutput) -> float:
        gd_factor = self._gd_factor()
        m_ini_t = self._norm_factor * sl_out.m_ini_tot * gd_factor

        return m_ini_t

    def _m_cor_t(self, sl_out: StarlightOutput) -> float:
        gd_factor = self._gd_factor()
        m_cor_t = self._norm_factor * sl_out.m_cor_tot * gd_factor

        return m_cor_t

    def _gd_factor(self):
        unit_conv = 4 * np.pi / 3.826E33
        gd_factor = unit_conv * (3.08567758E24 * self._galaxy_distance) ** 2.
        return gd_factor

    def _fc_component(self, sl_out: StarlightOutput, exp_min: float, exp_max: float) -> float:

        only_fc = np.array([j.lower().startswith("agn_fc_") for j in sl_out.component_j])
        
        fc_contribution_j = sl_out.x_j[only_fc]
        fc_name_j = sl_out.component_j[only_fc]

        no_extension = re.compile(r"\.[A-Za-z0-9]+$")

        valid_exp = []
        for fc in fc_name_j:
            fc_no_extension = re.sub(no_extension, "", fc.lower())
            exp_val = float(fc_no_extension.lstrip("agn_fc_"))
            valid_exp.append(exp_val > exp_min and exp_val <= exp_max)
        valid_exp = np.array(valid_exp)

        total_x_j = np.sum(sl_out.x_j)
        if total_x_j <= 0.:
            return 0
        sum_factor = 100. / total_x_j

        return np.sum(fc_contribution_j[valid_exp]) * sum_factor

    def _bb_component(self, sl_out: StarlightOutput, temp_min: float, temp_max: float) -> float:

        only_bb = np.array([j.lower().startswith("agn_bb_") for j in sl_out.component_j])
        
        bb_contribution_j = sl_out.x_j[only_bb]
        bb_name_j = sl_out.component_j[only_bb]

        no_extension = re.compile(r"\.[A-Za-z0-9]+$")

        valid_temp = []
        for bb in bb_name_j:
            bb_no_extension = re.sub(no_extension, "", bb.lower())
            temp_val = float(bb_no_extension.lstrip("agn_bb_"))
            valid_temp.append(temp_val > temp_min and temp_val <= temp_max)
        valid_temp = np.array(valid_temp)

        total_x_j = np.sum(sl_out.x_j)
        if total_x_j <= 0.:
            return 0
        sum_factor = 100. / total_x_j

        return np.sum(bb_contribution_j[valid_temp]) * sum_factor

    def _pop_by_light(self, sl_out: StarlightOutput, age_min: float, age_max: float) -> float:
        age_index = (sl_out.age_j > age_min) * (sl_out.age_j <= age_max)

        exclude_bb = np.array([not j.lower().startswith("agn_bb_") for j in sl_out.component_j])
        exclude_fc = np.array([not j.lower().startswith("agn_fc_") for j in sl_out.component_j])

        total_x_j = np.sum(sl_out.x_j)
        if total_x_j <= 0.:
            return 0
        sum_factor = 100. / total_x_j

        return np.sum(sl_out.x_j[age_index * exclude_bb * exclude_fc]) * sum_factor

    def _pop_by_mass(self, sl_out: StarlightOutput, age_min: float, age_max: float) -> float:
        age_index = (sl_out.age_j > age_min) * (sl_out.age_j <= age_max)
        
        exclude_bb = np.array([not j.lower().startswith("agn_bb") for j in sl_out.component_j])
        exclude_fc = np.array([not j.lower().startswith("agn_fc") for j in sl_out.component_j])
        
        total_m_j = np.sum(sl_out.m_cor_j)
        if total_m_j <= 0.:
            return 0
        sum_factor = 100. / total_m_j

        return np.sum(sl_out.m_cor_j[age_index * exclude_bb * exclude_fc]) * sum_factor

    def _sfr_at_age(self, sl_out: StarlightOutput, age_min: float, age_max: float) -> float:
        m_cor_t = self._m_cor_t(sl_out)

        exclude_bb = np.array([not j.lower().startswith("agn_bb") for j in sl_out.component_j])
        exclude_fc = np.array([not j.lower().startswith("agn_fc") for j in sl_out.component_j])

        age_range = age_max - age_min
        age_index = (sl_out.age_j > age_min) * (sl_out.age_j <= age_max)

        m_total_factor = m_cor_t / (100. * age_range)
        sfr_value = np.sum(sl_out.m_ini_j[age_index * exclude_bb * exclude_fc]) * m_total_factor

        return sfr_value

    def _wavelength_info_cards(self):
        init_wave = self._grid_generator.parameters.olsyn_ini
        d_wave = self._grid_generator.parameters.odlsyn

        cards = list()
        cards.append(fits.Card("CRPIX3", 1))
        cards.append(fits.Card("CRVAL3", init_wave))
        cards.append(fits.Card("CD3_3", d_wave))
        return cards

    def _redshifted_wave_array(self, flux_header: fits.Header, z_size: int) -> np.ndarray:
        delta_name = "CDELT3" if "CDELT3" in flux_header else "CD3_3"

        dt_wave = flux_header[delta_name]
        c_wave_position = flux_header["CRPIX3"] - 1
        c_wave_value = flux_header["CRVAL3"]

        ini_wave = c_wave_value - c_wave_position * dt_wave

        wave_array = ini_wave + np.array([x*dt_wave for x in range(0, z_size)])
        return wave_array / (1. + self._redshift)
