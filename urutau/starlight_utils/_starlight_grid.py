"""
    A collection of classes to make it easier to use/create
    grids for Starlight (from Cid Fernandes).
"""

import os
import re
from dataclasses import dataclass
import numpy as np


@dataclass
class GridParameters:
    """Data containing starlight grid parameters (for starlight v04)."""

    bases_dir: str = ""
    obs_dir: str = ""
    mask_dir: str = ""
    starlight_config: str = ""
    master_base_file: str = ""
    mask_file: str = ""
    out_dir: str = ""

    seed: int = 0
    llow_sn: float = 0.0
    lupp_sn: float = 0.0
    olsyn_ini: float = 0.0
    olsyn_fin: float = 0.0
    odlsyn: float = 0.0
    fscale_chi2: float = 0.0
    fit_fxk: str = ""
    has_errors: int = 0
    has_flags: int = 0
    reddening_law: str = ""
    vel_recession: float = 0.0
    vel_dispersion: float = 0.0


class GridReader():
    """
        Starlight grid reader.
    """

    def __init__(self, sl_grid: str) -> None:
        self._re_brackets = re.compile(r"\[.*?\]")
        self._re_spaces = re.compile(r"\s+")
        self._re_trim = re.compile(r"^[\s]*|[\s]*$")

        self._input_entries = list()

        self._file_path = sl_grid
        self._sl_param = GridParameters()

        self._get_parameters()

    @property
    def parameters(self) -> GridParameters:
        """
            General parameters used by the grids.
        """
        return self._sl_param

    @property
    def input_entries(self) -> list:
        """
            Input entries and specific parameters used.
        """
        return self._input_entries.copy()

    def _get_parameters(self) -> GridParameters:

        if (not os.path.exists(self._file_path)) or os.path.isdir(self._file_path):
            return None

        with open(self._file_path, "r", encoding="utf8") as sl_file:
            lines = sl_file.readlines()

        num_entries = self._get_value(lines[0], int)

        self._read_parameters(lines)
        self._read_input_list(lines, num_entries)

    def _read_parameters(self, lines: list[str]) -> None:
        self._sl_param.bases_dir = self._get_value(lines[1], str)
        self._sl_param.obs_dir = self._get_value(lines[2], str)
        self._sl_param.mask_dir = self._get_value(lines[3], str)
        self._sl_param.out_dir = self._get_value(lines[4], str)

        self._sl_param.seed = self._get_value(lines[5], int)
        self._sl_param.llow_sn = self._get_value(lines[6], float)
        self._sl_param.lupp_sn = self._get_value(lines[7], float)
        self._sl_param.olsyn_ini = self._get_value(lines[8], float)
        self._sl_param.olsyn_fin = self._get_value(lines[9], float)
        self._sl_param.odlsyn = self._get_value(lines[10], float)
        self._sl_param.fscale_chi2 = self._get_value(lines[11], float)
        self._sl_param.fit_fxk = self._get_value(lines[12], str)
        self._sl_param.has_errors = self._get_value(lines[13], int)
        self._sl_param.has_flags = self._get_value(lines[14], int)

        self._sl_param.starlight_config = ""
        self._sl_param.master_base_file = ""
        self._sl_param.mask_file = ""

        self._sl_param.reddening_law = ""
        self._sl_param.vel_recession = 0.0
        self._sl_param.vel_dispersion = 0.0

    def _read_input_list(self, lines: list[str], num_entries: int) -> None:
        skip = 15
        input_lines = lines[skip:skip+num_entries]
        for line in input_lines:
            input_file = self._get_value(line, str, 0)

            starlight_config = self._get_value(line, str, 1)
            master_base_file = self._get_value(line, str, 2)
            mask_file = self._get_value(line, str, 3)

            reddening_law = self._get_value(line, str, 4)
            vel_recession = self._get_value(line, float, 5)
            vel_dispersion = self._get_value(line, float, 6)

            output_file = self._get_value(line, str, 7)

            input_data = dict()

            input_data["input"] = input_file
            input_data["config"] = starlight_config
            input_data["base"] = master_base_file
            input_data["mask"] = mask_file
            input_data["reddening"] = reddening_law
            input_data["v0"] = vel_recession
            input_data["vd"] = vel_dispersion
            input_data["output"] = output_file

            self._input_entries.append(input_data)

    def _get_value(self, line: str, dtype: str, index: int = 0) -> int:
        value_str = re.sub(self._re_brackets, "", line)
        value_str = re.sub(self._re_trim, "", value_str)
        values = re.split(self._re_spaces, value_str)
        return dtype(values[index])


class GridGenerator():
    """
        Generate grid files with selected parameters (for starlight v04).
    """

    def __init__(self, sl_par: GridParameters, grid_prefix: str, grid_dir: str = "./") -> None:
        self._sl_par = sl_par
        self._prepare_paths()

        self._grid_prefix = grid_prefix
        self._grid_dir = grid_dir
        self._generated_grids = list()

    @property
    def parameters(self) -> GridParameters:
        """
            Parameters used to create the grids.
        """
        return self._sl_par

    @property
    def grids(self) -> list[str]:
        """
            List with all generated grids.
        """
        return self._generated_grids.copy()

    def change_grid_parameters(self, sl_par: GridParameters) -> None:
        """
            Change current grid parameters.
        """
        self._sl_par = sl_par
        self._prepare_paths()

    def change_grid_prefix(self, grid_prefix: str) -> None:
        """
            Change current grid prefix.
        """
        self._grid_prefix = grid_prefix

    def change_grid_dir(self, grid_dir: str) -> None:
        """
            Change current grid dir.
        """
        self._grid_dir = grid_dir

    def clear_grids(self) -> None:
        """
            Clear current grid list without deleting the files.
        """
        self._generated_grids.clear()

    def delete_grids(self) -> None:
        """
            Delete all generated grid files.
        """
        for grid in self._generated_grids:
            os.remove(grid)

        self.clear_grids()

    def create_grid(self, specs: np.ndarray | str) -> np.ndarray:
        """
            Generate grid file with the input spectra.

            Returns list with possible starlight output destinations.
        """

        if isinstance(specs, str):
            specs = np.array([specs])
        else:
            specs = np.array(specs)

        number_of_grids = len(self._generated_grids)
        grid_name = f"{self._grid_prefix}_{number_of_grids}.inp"

        grid_path = os.path.join(self._grid_dir, grid_name)

        if not os.path.exists(self._grid_dir):
            os.makedirs(self._grid_dir)

        out_list = list()

        spec_list = specs.flat

        with open(grid_path, "w", encoding="utf8") as save:
            input_number = len(spec_list)
            base_parameters = self._base_lines(input_number)
            save.write(base_parameters)

            for spec in spec_list:
                if spec is None or spec == "":
                    out_list.append("")
                    continue

                base_name = os.path.basename(spec)
                output_name = f"out_{base_name}"
                input_line = self._spec_line(base_name, output_name)
                save.write(input_line)

                file_path = os.path.join(self.parameters.out_dir, output_name)
                out_list.append(file_path)

            self._generated_grids.append(grid_path)

            return np.array(out_list, dtype=str).reshape(np.shape(specs))

    def get_list_of_output_destinations(self) -> None:
        """
            Get the list of output result expected from the current grids.
        """
        list_output = list()
        grid_lines = list()

        re_spaces = re.compile(r"\s+")
        re_trim = re.compile(r"^[\s]*|[\s]*$")

        for grid in self._generated_grids:
            with open(grid, "r", encoding="utf8") as grid_file:
                current_lines = grid_file.readlines()[15:]
                grid_lines += current_lines

        for line in grid_lines:
            trimmed_line = re.sub(re_trim, "", line)
            output_name = re.split(re_spaces, trimmed_line)[-1]
            output_path = os.path.join(self._sl_par.out_dir, output_name)
            list_output.append(output_path)

        return list_output

    def _base_lines(self, input_number):
        base_format = "{0}\n{1}\n{2}\n{3}\n{4}\n{5}\n{6}\n{7}\n"
        base_format += "{8}\n{9}\n{10}\n{11}\n{12}\n{13}\n{14}\n"

        base_parameters = base_format.format(input_number, self._sl_par.bases_dir,
                                             self._sl_par.obs_dir, self._sl_par.mask_dir,
                                             self._sl_par.out_dir, self._sl_par.seed,
                                             self._sl_par.llow_sn, self._sl_par.lupp_sn,
                                             self._sl_par.olsyn_ini, self._sl_par.olsyn_fin,
                                             self._sl_par.odlsyn, self._sl_par.fscale_chi2,
                                             self._sl_par.fit_fxk, self._sl_par.has_errors,
                                             self._sl_par.has_flags)
        return base_parameters

    def _spec_line(self, spec_name: str, output_name: str):

        line_format = "{0} {1} {2} {3} {4} {5} {6} {7}\n"
        line = line_format.format(spec_name, self._sl_par.starlight_config,
                                  self._sl_par.master_base_file, self._sl_par.mask_file,
                                  self._sl_par.reddening_law, self._sl_par.vel_recession,
                                  self._sl_par.vel_dispersion, output_name)

        return line

    def _prepare_paths(self) -> None:
        self._sl_par.bases_dir = self._abs_path(self._sl_par.bases_dir) + "/"
        self._sl_par.obs_dir = self._abs_path(self._sl_par.obs_dir) + "/"
        self._sl_par.mask_dir = self._abs_path(self._sl_par.mask_dir) + "/"
        self._sl_par.out_dir = self._abs_path(self._sl_par.out_dir) + "/"

    def _abs_path(self, path: str) -> str:
        expanded = os.path.expanduser(path)
        absolute = os.path.abspath(expanded)
        return absolute
