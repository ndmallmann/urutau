"""
    Urutau, a module based pipeline.
"""

import os
import queue
import threading as th
from typing import Type

import astropy.io.fits as fits
import pandas as pd

from .modules._module_base import AbstractModule


class Urutau:
    """
        Urutau is a module based pipeline.
    """

    def __init__(self) -> None:
        self._modules = list()
        self._targets = list()

        self._modules_configs = dict()
        self._default_parameters = dict()
        self._specific_parameters = dict()

        self._num_procs = 1
        self._executing = False

        self._jobs = queue.Queue()

    def set_modules(self, modules: list[Type[AbstractModule]]) -> None:
        """
            Set modules list to be executed.

            Replaces old list and clears module configs.
        """

        if self._executing:
            return

        self._modules = modules.copy()
        self._modules_configs.clear()

    def config_module(self, module: Type[AbstractModule], configuration: dict) -> None:
        """
            Set initial configuration for module.
        """

        if self._executing:
            return

        self._modules_configs[module] = configuration

    def clear_targets(self) -> None:
        """
            Clear current target list.
        """

        if self._executing:
            return

        self._targets.clear()

    def clear_default_parameters(self) -> None:
        """
            Clear current default parameters.
        """

        if self._executing:
            return

        self._default_parameters.clear()

    def clear_target_specific_parameters(self) -> None:
        """
            Clear current specific parameters.
        """

        if self._executing:
            return

        self._specific_parameters.clear()

    def load_targets(self, targets: list[str]) -> None:
        """
            Load new target list to replace old ones.

            Warning: Duplicates will be removed due to racing conditions with multiple threads.
        """

        if self._executing:
            return

        for target in targets:
            if target not in self._targets:
                self._targets.append(target)

    def load_default_parameters(self, default_parameters: dict) -> None:
        """
            Load new set of default parameters.

            Default parameters are used by all targets.

            Ex:
               default_parameters = {
                "reddening law": "CCM",
                "initial velocity guess": 150,
                }
        """

        if self._executing:
            return

        self._default_parameters.update(default_parameters)

    def load_target_specific_parameters(self, specific_parameters: dict) -> None:
        """
            Load new set of specific parameters.

            Specific parameters are used to change default parameter values for specific targets.

            Ex:
               specific_parameters = {
                "./target1.fits": {"speed": 10},
                "./target2.fits": {"redshift": 0.1, "speed": 5},
                }

            Warning: Method does nothing if urutau is executing.
        """

        if self._executing:
            return

        self._specific_parameters.update(specific_parameters)

    def set_number_of_threads(self, num_threads: int) -> None:
        """
            Set number of working threads to process multiple targets at the same time.

            Warning: num_threads can't be less than 1
        """

        if self._executing:
            return

        self._num_procs = num_threads if num_threads > 1 else 1

    def execute(self, save_path_root: str = "./", save_config: bool = True, debug: bool = False) -> None:
        """
            Execute all modules from the list.

            Saves the final result on save_path.
        """

        if len(self._modules) == 0:
            print("No modules loaded!")
            return
        elif len(self._targets) == 0:
            print("No targets loaded!")
            return
        elif self._executing:
            print("Urutau already executing!")
            return
        elif os.path.exists(save_path_root) and not os.path.isdir(save_path_root):
            print(f"Invalid directory path '{save_path_root}'!")
            return

        if not os.path.exists(save_path_root):
            os.makedirs(save_path_root)

        self._executing = True

        for target in self._targets:
            self._jobs.put(target)

        for _ in range(self._num_procs):
            worker_args = (save_path_root, save_config, debug)
            worker = th.Thread(target=self._worker_task, args=worker_args)
            worker.start()

        self._jobs.join()
        self._targets.clear()
        self._executing = False

    def read_csv(self, targets_dir: str, csv_file: str) -> None:
        """
            Reads a csv file with the targets and their specific configurations.

            targets_dir: folder path containing the targets
            csv_file: column named "target" is used to identify the targets (otherwise, the first column is used)
        """

        if self._executing:
            return

        csv_df = pd.read_csv(csv_file)

        target_column = "target"

        if target_column in csv_df.columns:
            targets_names = list(csv_df[target_column].values)
            target_property = target_column
            properties = csv_df.columns.to_list()
            properties.remove(target_column)
        else:
            targets_names = list(csv_df.iloc[:, 0].values)
            columns_list = csv_df.columns.to_list()
            target_property = columns_list[0]
            properties = columns_list[1:]

        targets = [os.path.join(targets_dir, x) for x in targets_names]
        self.load_targets(targets)

        spec_config = {
            os.path.join(targets_dir, row[target_property]):
            row[properties].to_dict() for (_, row) in csv_df.iterrows()
        }
        self.load_target_specific_parameters(spec_config)

    def _worker_task(self, save_path_root: str, save_config: bool, debug: bool) -> None:

        while True:

            # Verify if the jobs queue is not empty
            try:
                target = self._jobs.get(False)
            except queue.Empty:
                break

            final_target_name = self._final_save_path(save_path_root, target)

            orig_file_data = fits.open(target)

            final_configuration = dict()

            for next_module in self._modules:

                config_parameters = dict()

                # Config priority order:
                #       module -> default targets -> specific target
                if next_module in self._modules_configs:
                    config_parameters.update(
                        self._modules_configs[next_module])

                config_parameters.update(self._default_parameters.copy())

                if target in self._specific_parameters:
                    config_parameters.update(self._specific_parameters[target])

                loaded_module = next_module(**config_parameters)

                if debug:
                    self._debug_message_parameters(loaded_module)

                # Execute module
                next_hdu_list = loaded_module.execute(orig_file_data)

                self._update_result_file(orig_file_data, next_hdu_list)

                final_configuration.update(loaded_module.config)

            if debug:
                self._debug_message_final_configuration(
                    target, final_configuration)

            if save_config:
                config_hdu = self._generate_config_hdu(final_configuration)
                self._update_result_file(orig_file_data, config_hdu)

            orig_file_data.writeto(final_target_name, overwrite=True)

            self._jobs.task_done()

    def _update_result_file(self, data: fits.HDUList, new_data: fits.HDUList) -> None:
        for hdu in new_data:
            data.append(hdu)

    def _final_save_path(self, save_path_root, target):
        target_base_name = os.path.basename(target)
        file_name = os.path.splitext(target_base_name)[0] + "_urutau.fits"
        save_path = os.path.join(save_path_root, file_name)
        return save_path

    def _generate_config_hdu(self, configuration: dict) -> fits.HDUList:
        hdu_list = fits.HDUList()

        config_hdu = fits.ImageHDU()
        header = config_hdu.header

        header["EXTNAME"] = "URUT_CFG"

        index_par = 0
        for key, value in configuration.items():
            card = fits.Card(
                keyword=f"CFGP{index_par}", value=str(value), comment=key)
            header.append(card=card)
            index_par += 1

        header.add_comment("Urutau Config Parameters")

        hdu_list.append(config_hdu)

        return hdu_list

    def _debug_message_final_configuration(self, target: str, final_configuration: dict) -> None:
        message = f"\n>>> Target {target}\n___ Final Configuration: {final_configuration}"
        print(message)

    def _debug_message_parameters(self, loaded_module: AbstractModule) -> None:
        print(f"\n>>> Module Loaded: {loaded_module.name}")
        print(f"___ Default Par: {loaded_module.default_parameters}")
        print(f"___ Received Par: {loaded_module.received_config}")
        print(f"___ Config Par: {loaded_module.config}")
