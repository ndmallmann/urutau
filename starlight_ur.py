import os
import math
import uuid
import numpy as np
import shutil as sh
import subprocess as sp
from urutau import UrutauModule
import astropy.io.fits as fits


class StarlightWrap(UrutauModule):

    name = "Starlight Wrap"

    def _set_init_default_parameters(self) -> None:
        self.default_parameters["flux extension"] = "FLUX"
        self.default_parameters["central lambda position"] = "CRPIX3"
        self.default_parameters["central lambda value"] = "CRVAL3"
        self.default_parameters["delta lambda"] = "CD3_3"
        self.default_parameters["redshift"] = 0.0
        self.default_parameters["starlight directory path"] = "./starlight/"
        self.default_parameters["starlight executable"] = "StarlightChains_v04.amd64_gfortran-4.1.1_static.exe"
        self.default_parameters["bases directory"] = "./starlight/BasesDir/"
        self.default_parameters["mask path"] = "./starlight/Mask/"
        self.default_parameters["seed"] = 123456
        self.default_parameters["llow sn"] = 5650.0
        self.default_parameters["lupp sn"] = 5750.0
        self.default_parameters["olsyn ini"] = 3800.0
        self.default_parameters["olsyn fin"] = 7000.0
        self.default_parameters["odlsyn"] = 1.0
        self.default_parameters["fscale chi2"] = 1.0
        self.default_parameters["fit fxk"] = "FIT"
        self.default_parameters["iserrspecavailable"] = 0
        self.default_parameters["isflagspecavailable"] = 0
        self.default_parameters["starlight config"] = "StCv04.C11.config"
        self.default_parameters["master base file"] = "Base.BC03.N"
        self.default_parameters["mask file"] = "Masks.EmLines.SDSS.gm"
        self.default_parameters["reddening law"] = "CCM"
        self.default_parameters["velocity recession"] = 0.
        self.default_parameters["velocity dispersion"] = 150.

        self.default_parameters["galaxy distance"] = 0
        self.default_parameters["normalization factor"] = 1e-17

        self.default_parameters["signal to noise list"] = (1, 5, 10, 20)
        self.default_parameters["signal to noise window"] = (5650, 5750)
        self.default_parameters["ivar extension"] = None
        self.default_parameters["error extension"] = None

        self.default_parameters["BinFCVec"] = {"FC1.50": [1.49, 1.5]}
        self.default_parameters["BinHDVec"] = {"BB_h": [1001, 1500],
                                               "BB_c": [0, 1000]}
        self.default_parameters["BinPopVec"] = {"xyy": [0.9E6, 10.1E6], "xyo": [14.0E6, 56.3E6],
                                                "xiy": [99.9E6, 502.0E6], "xii": [630.0E6, 795.0E6],
                                                "xio": [890.0E6, 2.01E9], "xo": [5.0E9, 12.7E9]}
        self.default_parameters["BinSFR"] = {"SFR_1E6": [10, 1.001E6], "SFR_5E6": [10, 5.621E6],
                                             "SFR_10E6": [10, 10.001E6], "SFR_14E6": [10, 14.1001E6],
                                             "SFR_20E6": [10, 20.001E6], "SFR_30E6": [10, 31.6001E6],
                                             "SFR_56E6": [10, 56.201E6], "SFR_100E6": [10, 100.001E6],
                                             "SFR_200E6": [10, 200.001E6]}

    def execute(self, target_file: str) -> fits.FitsHDU:

        dir_files, list_files = self._extract_files(target_file)

        grid_file = self._generate_grid_files(dir_files, list_files)

        self._run_starlight(grid_file)

        # sh.rmtree(dir_files)

        return fits.ImageHDU()

    def _extract_files(self, target: str) -> tuple[str, list[str]]:

        tmp_location = self._generate_tmp_directory()

        target_data = fits.getdata(target, self._config["flux extension"])

        _, y_size, x_size = target_data.shape
        wavelength = self._wavelength_array(target)

        data_format = ["%.2f", "%.3e"]
        data_delimiter = " "

        list_files = []
        for i in range(x_size):
            for j in range(y_size):
                file_name = os.path.join(tmp_location, f"tmp_{i}_{j}.spec")

                spec_flux = target_data[:, j, i]

                if np.isnan(np.sum(spec_flux)) or math.isclose(np.mean(spec_flux), 0.):
                    continue

                file_data = np.column_stack((wavelength, spec_flux))
                np.savetxt(file_name, file_data, fmt=data_format,
                           delimiter=data_delimiter)

                list_files.append(file_name)

        return tmp_location, list_files

    def _generate_grid_files(self, dir_files: str, list_spec: list[str]) -> str:
        base_grid_format = "{0}\n{1}\n{2}\n{3}\n{4}\n{5}\n{6}\n{7}\n{8}\n{9}\n{10}\n{11}\n{12}\n{13}\n{14}\n"
        line_format = "{0} {1} {2} {3} {4} {5} {6} {7}\n"

        dir_files = self._abs_path(dir_files) + "/"
        bases_dir = self._abs_path(self["bases directory"]) + "/"
        mask_path = self._abs_path(self["mask path"]) + "/"

        # Create a grid file
        grid_file_name = "grid.inp"
        grid_file_path = os.path.join(dir_files, grid_file_name)

        with open(grid_file_path, "w", encoding="utf8") as save:

            base_parameters = base_grid_format.format(len(list_spec), bases_dir, dir_files, mask_path, dir_files,
                                                      self["seed"], self["llow sn"], self["lupp sn"], self["olsyn ini"],
                                                      self["olsyn fin"], self["odlsyn"], self["fscale chi2"], self["fit fxk"],
                                                      self["iserrspecavailable"], self["isflagspecavailable"])
            save.write(base_parameters)

            for spec in list_spec:
                spec_name = os.path.basename(spec)
                output_name = f"starlight_{spec_name}"
                save.write(line_format.format(spec_name, self["starlight config"],
                                              self["master base file"], self["mask file"],
                                              self["reddening law"], self["velocity recession"],
                                              self["velocity dispersion"], output_name))

        return grid_file_path

    def _run_starlight(self, grid_file_path: str) -> None:

        with open(grid_file_path, "r", encoding="utf8") as grid:
            directory = self._abs_path(self["starlight directory path"])
            executable = os.path.join("./", self["starlight executable"])
            arguments = (executable,)

            prog = sp.Popen(args=arguments, stdin=grid, cwd=directory)

            prog.wait()

    def _abs_path(self, path: str) -> str:
        expanded = os.path.expanduser(path)
        absolute = os.path.abspath(expanded)
        return absolute

    def _generate_tmp_directory(self) -> str:
        unique_id = uuid.uuid4().hex
        tmp_location = f"./{unique_id}/"

        while os.path.exists(tmp_location):
            unique_id = uuid.uuid4().hex
            tmp_location = f"./{unique_id}/"

        os.makedirs(tmp_location)

        return tmp_location

    def _wavelength_array(self, target: str) -> np.ndarray:
        hdu_header = fits.getheader(target, self["flux extension"])

        clp = hdu_header[self["central lambda position"]] - 1
        clv = hdu_header[self["central lambda value"]]
        delta_l = hdu_header[self["delta lambda"]]

        z_size, _, _ = fits.getdata(target, self["flux extension"]).shape

        start_wave = clv - clp * delta_l
        wave = [x * delta_l + start_wave for x in range(z_size)]
        return np.array(wave) / (1. + self["redshift"])

    def _popbins_header(self) -> fits.Header:

        popbins_header = []

        data_index = 0

        aux = self["BinFCVec"]
        for w in aux:
            comment = f"Featureless cont. F_l {aux[w][0]} & {str(aux[w][1])}"
            card = fits.Card(f"DATA{data_index}", f"FC{w}", comment)
            popbins_header.append(card)
            data_index += 1

        aux = self["BinHDVec"]
        for w in aux:
            comment = f"Planck function T: {aux[w][0]} & {aux[w][1]}"
            card = fits.Card(f"DATA{data_index}", f"BB{w}", comment)
            popbins_header.append(card)
            data_index += 1

        aux = self["BinPopVec"]
        for w in aux:
            comment = f"Light binned pop. vec {aux[w][0]:.1E} < age <= {aux[w][1]:.1E}"
            card = fits.Card(f"DATA{data_index}", f"{w}_light", comment)
            popbins_header.append(card)
            data_index += 1

        for w in aux:
            comment = f"Mass binned pop. vec {aux[w][0]:.1E} < age <= {aux[w][1]:.1E}"
            card = fits.Card(f"DATA{data_index}", f"{w}_mass", comment)
            popbins_header.append(card)
            data_index += 1

        aux = self["BinSFR"]
        for w in aux:
            comment = f"SFR for {aux[w][0]:.1E} < age <= {aux[w][1]:.1E}"
            card = fits.Card(f"DATA{data_index}", f"{w}", comment)
            popbins_header.append(card)
            data_index += 1

        ext_keys = [['Av', 'Optical extinction'],
                    ['Mage_L', 'Mean age light weigthed'],
                    ['Mage_M', 'Mean age mass weigthed'],
                    ['MZ_L', 'Mean metalicity light weigthed'],
                    ['MZ_M', 'Mean metalicity mass weigthed'],
                    ['Mstar', 'Mass in stars (M_sun, M* from starligh)'],
                    ['Mpross', 'Mass converted to stars (~2xMstar)'],
                    ['F_Norm', 'Normalization flux in input units'],
                    ['Sigma_star', 'Stellar dispersion'],
                    ['vrot_star', 'Stellar rotation'],
                    ['Adev', 'Precentage mean deviation (see manual)'],
                    ['ChiSqrt', 'ChiSqrt/Nl_eff (see manual)'],
                    ['SNR', 'SNR on normalization window']]

        for key, comment in ext_keys:
            card = fits.Card(f"DATA{data_index}", key, comment)
            popbins_header.append(card)
            data_index += 1

        return fits.Header(cards=popbins_header)
