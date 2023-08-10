"""
    Urutau, a module based pipeline.
"""

import os
import queue
import threading as th

import astropy.io.fits as fits
import pandas as pd

from ._containers import ModuleContainer, TargetContainer
from .modules._module_base import AbstractModule

_URUTAU_CONFIG_EXTNAME = "URUT_CFG"
_URUTAU_PATH_ENDING = "_urutau.fits"


class Urutau:
    """
        Urutau is a module based pipeline.
    """

    def __init__(self, num_threads: int = 1) -> None:
        self._module_containers: list[ModuleContainer] = list()
        self._target_containers: list[TargetContainer] = list()

        self._num_threads = int(num_threads) if num_threads > 1 else 1

        self._save_path_root = "./"
        self._executing = False

    def add_module(self, module: AbstractModule, config: dict = None) -> None:
        """
            Add module at the end of the execution list.
        """

        self._check_if_executing()

        module_cfg = dict() if config is None else config
        module_container = ModuleContainer(module, module_cfg)

        self._module_containers.append(module_container)

    def add_target(self, target: str, config: dict = None) -> None:
        """
            Add target to the list.
        """

        self._check_if_executing()

        if not os.path.exists(target) or os.path.isdir(target):
            raise RuntimeError(f"'{target}' is not a valid target path.")

        target_cfg = dict() if config is None else config
        target_container = TargetContainer(target, target_cfg)

        self._target_containers.append(target_container)

    def read_csv(self, targets_dir: str, csv_file: str) -> None:
        """
            Reads a csv file with the targets and their specific configurations.

            targets_dir: folder path containing the targets
            csv_file: column named "target" is used to identify the targets (otherwise, the first column is used)
        """

        self._check_if_executing()

        csv_df: pd.DataFrame = pd.read_csv(csv_file)
        properties = csv_df.columns.to_list()

        target_column = "target"

        if target_column in properties:
            properties.remove(target_column)
        else:
            target_column = properties[0]
            properties.remove(target_column)

        for _, row in csv_df.iterrows():
            target_path = os.path.join(targets_dir, row[target_column])
            target_cfg = {prop: row[prop] for prop in properties}
            self.add_target(target_path, target_cfg)

    def execute(self, save_path_root: str = "./", save_config: bool = True, debug: bool = False) -> None:
        """
            Execute all modules from the list (in sequential order).

            Saves the final result on save_path.
        """

        self._save_path_root = save_path_root

        self._runtime_checks()
        self._create_directory()

        self._executing = True

        jobs = self._create_jobs_list()
        
        for worker in self._create_workers(jobs, save_config, debug):
            worker.start()

        jobs.join()

        self._target_containers.clear()
        self._executing = False

    def _create_workers(self, jobs: queue.Queue[TargetContainer], save_config: bool, debug: bool) -> list[th.Thread]:
        workers: list[th.Thread] = list()

        for _ in range(self._num_threads):
            worker_args = (jobs, self._module_containers.copy(),
                           self._save_path_root, save_config, debug)
            worker = th.Thread(target=_worker_task, args=worker_args)
            workers.append(worker)

        return workers

    def _create_jobs_list(self) -> queue.Queue[TargetContainer]:
        jobs = queue.Queue()

        for target_container in self._target_containers:
            jobs.put(target_container)

        return jobs

    def _create_directory(self):
        if not os.path.exists(self._save_path_root):
            os.makedirs(self._save_path_root)

    def _runtime_checks(self) -> None:
        self._check_if_executing()
        self._check_modules()
        self._check_targets()
        self._check_if_valid_path()

    def _check_if_executing(self) -> None:
        if self._executing:
            raise RuntimeError("Urutau is running.")

    def _check_modules(self) -> None:
        if not len(self._module_containers) > 0:
            raise RuntimeError(f"No modules loaded.")

    def _check_targets(self) -> None:
        if not len(self._target_containers) > 0:
            raise RuntimeError(f"No targets loaded.")

    def _check_if_valid_path(self) -> None:
        if os.path.exists(self._save_path_root) and not os.path.isdir(self._save_path_root):
            print(f"Invalid directory path '{self._save_path_root}'!")


def _worker_task(jobs: queue.Queue[TargetContainer], module_containers: list[ModuleContainer], path_root: str, save_config: bool, debug: bool) -> None:

    while True:
        # Verify if the jobs queue is not empty
        try:
            target_container: TargetContainer = jobs.get(False)
        except queue.Empty:
            return

        final_configuration = dict()

        if not os.path.exists(target_container.target):
            print(f"TARGET {target_container.target} NOT FOUND!")
            jobs.task_done()

        with fits.open(target_container.target) as opened_file:

            for module_capsule in module_containers:

                loaded_module = module_capsule.get_loaded_module(
                    target_container.config)

                if debug:
                    _debug_message_parameters(loaded_module)

                result_hdus = loaded_module.execute(opened_file)
                _concatenate_hdus(opened_file, result_hdus)

                final_configuration.update(loaded_module.config)

            if debug:
                _debug_message_final_configuration(
                    target_container.target, final_configuration)

            if save_config:
                cfg_hdus = _generate_config_hdu(final_configuration)
                _concatenate_hdus(opened_file, cfg_hdus)

            final_target_name = _final_save_path(
                path_root, target_container.target)
            opened_file.writeto(final_target_name, overwrite=True)

        jobs.task_done()


def _concatenate_hdus(hdus: fits.HDUList, extra_hdus: fits.HDUList) -> None:
    for hdu in extra_hdus:
        hdus.append(hdu)


def _final_save_path(save_path_root: str, target: str):
    base_name = os.path.basename(target)
    base_name_no_ext, _ = os.path.splitext(base_name)

    final_name = f"{base_name_no_ext}{_URUTAU_PATH_ENDING}"
    save_path = os.path.join(save_path_root, final_name)

    return save_path


def _generate_config_hdu(configuration: dict) -> fits.HDUList:
    config_hdu = fits.ImageHDU()

    header = config_hdu.header
    header["EXTNAME"] = _URUTAU_CONFIG_EXTNAME

    index_par = 0
    for key, value in configuration.items():
        card = fits.Card(
            keyword=f"CFGP{index_par}",
            value=key,
            comment=str(value)
        )
        header.append(card=card)
        index_par += 1

    header.add_comment("Urutau Config Parameters")

    return fits.HDUList(hdus=[config_hdu])


def _debug_message_final_configuration(target: str, final_configuration: dict) -> None:
    message = f"\n>>> Target {target}\n___ Final Configuration: {final_configuration}"
    print(message)


def _debug_message_parameters(loaded_module: AbstractModule) -> None:
    print(f"\n>>> Module Loaded: {loaded_module.name}")
    print(f"___ Default Par: {loaded_module.default_parameters}")
    print(f"___ Received Par: {loaded_module.received_config}")
    print(f"___ Run Config Par: {loaded_module.config}")
