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


class ModuleContainer:

    def __init__(self, module: Type[AbstractModule], config: dict) -> None:
        self._module = module
        self._config = config

    def get_module(self) -> Type[AbstractModule]:
        return self._module

    def get_config(self) -> dict:
        return self._config.copy()

    def loaded_module(self, target_config: dict = None) -> AbstractModule:
        if target_config is None:
            target_config = dict()

        config_copy = self._config.copy()
        config_copy.update(target_config)

        loaded_module = self._module(**config_copy)

        return loaded_module


class TargetContainer:

    def __init__(self, target: str, config: dict = None) -> None:
        self._target = target
        self._config = dict() if config is None else config

    def get_target(self) -> str:
        return self._target

    def get_config(self) -> dict:
        return self._config.copy()


class Urutau:
    """
        Urutau is a module based pipeline.
    """

    def __init__(self, num_threads: int = 1) -> None:
        self._module_containers: list[ModuleContainer] = list()
        self._target_containers: list[TargetContainer] = list()

        self._num_threads = int(num_threads) if num_threads > 1 else 1
        self._jobs: queue.Queue[TargetContainer] = queue.Queue()

        self._executing = False

    def add_module(self, module: AbstractModule, config: dict = None) -> None:
        """
            Add module at the end of the execution list.
        """

        if self._executing:
            return

        module_cfg = dict() if config is None else config
        module_container = ModuleContainer(module, module_cfg)

        self._module_containers.append(module_container)

    def add_target(self, target: str, config: dict = None) -> None:
        """
            Add target to the list.
        """

        if self._executing:
            return

        target_cfg = dict() if config is None else config
        target_container = TargetContainer(target, target_cfg)

        self._target_containers.append(target_container)

    def read_csv(self, targets_dir: str, csv_file: str) -> None:
        """
            Reads a csv file with the targets and their specific configurations.

            targets_dir: folder path containing the targets
            csv_file: column named "target" is used to identify the targets (otherwise, the first column is used)
        """

        if self._executing:
            return

        csv_df: pd.DataFrame = pd.read_csv(csv_file)
        properties = csv_df.columns.to_list()

        target_column = "target"

        if target_column in properties:
            properties.remove(target_column)
        else:
            target_column = properties[0]
            properties.remove(target_column)

        for index, row in csv_df.iterrows():
            target_path = os.path.join(targets_dir, row[target_column])
            target_cfg = {prop: row[prop] for prop in properties}
            self.add_target(target_path, target_cfg)

    def execute(self, save_path_root: str = "./", save_config: bool = True, debug: bool = False) -> None:
        """
            Execute all modules from the list.

            Saves the final result on save_path.
        """

        if not self._can_run(save_path_root):
            return

        if not os.path.exists(save_path_root):
            os.makedirs(save_path_root)

        self._executing = True

        for target_container in self._target_containers:
            self._jobs.put(target_container)

        for _ in range(self._num_threads):
            worker_args = (save_path_root, save_config, debug)
            worker = th.Thread(target=self._worker_task, args=worker_args)
            worker.start()

        self._jobs.join()
        self._target_containers.clear()

        self._executing = False

    def _can_run(self, save_path_root: str) -> None:

        is_valid = True

        if len(self._module_containers) == 0:
            print("No modules loaded!")
            is_valid = False
        elif len(self._target_containers) == 0:
            print("No targets loaded!")
            is_valid = False
        elif self._executing:
            print("Urutau already executing!")
            is_valid = False
        elif os.path.exists(save_path_root) and not os.path.isdir(save_path_root):
            print(f"Invalid directory path '{save_path_root}'!")
            is_valid = False

        return is_valid

    def _worker_task(self, path_root: str, save_config: bool, debug: bool) -> None:

        while True:

            # Verify if the jobs queue is not empty
            try:
                target_container = self._jobs.get(False)
            except queue.Empty:
                return

            final_configuration = dict()

            target = target_container.get_target()
            target_cfg = target_container.get_config()

            if not os.path.exists(target):
                print(f"TARGET {target} NOT FOUND!")
                self._jobs.task_done()

            with fits.open(target) as opened_file:

                for module_capsule in self._module_containers:

                    loaded_module = module_capsule.loaded_module(target_cfg)

                    if debug:
                        self._debug_message_parameters(loaded_module)

                    result_hdus = loaded_module.execute(opened_file)
                    self._update_result_file(opened_file, result_hdus)

                    final_configuration.update(loaded_module.config)

                if debug:
                    self._debug_message_final_configuration(
                        target, final_configuration)

                if save_config:
                    cfg_hdus = self._generate_config_hdu(final_configuration)
                    self._update_result_file(opened_file, cfg_hdus)

                final_target_name = self._final_save_path(path_root, target)
                opened_file.writeto(final_target_name, overwrite=True)

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
                keyword=f"CFGP{index_par}", value=key, comment=str(value))
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
