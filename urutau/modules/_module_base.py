"""
    Abstract Module for Urutau.
"""

from abc import ABC, abstractmethod
from typing import Any

import astropy.io.fits as fits


class AbstractModule(ABC):
    """
        Abstract Module for Urutau.
        All urutau modules should derive from this one.
    """

    name: str = "MOCK"

    def __init__(self, **kwargs) -> None:
        self._config = dict()
        self._init_config = dict()
        self._default_config = dict()

        self._save_init_config_parameters(**kwargs)
        self._set_init_default_parameters()

        self._config = self._default_config.copy()
        self.setup(**self._init_config)

    @ property
    def default_parameters(self) -> dict:
        """
            Default parameters for the module.
        """
        return self._default_config

    @ property
    def received_config(self) -> dict:
        """
            Initial parameters (used during object creation).
        """
        return self._init_config

    @ property
    def config(self) -> dict:
        """
            Current parameters being used.
        """
        return self._config

    def __getitem__(self, key: str) -> Any:
        return self._config[key]

    def setup(self, **kwargs) -> None:
        """
            Modify the default configuration values for the module.
        """

        for key, value in kwargs.items():
            if key in self._config:
                self._config[key] = value

    def reset_parameters(self) -> None:
        """
            Resets all parameter values to the ones used during object creation.
        """
        self._config.clear()
        self._set_init_default_parameters()
        self.setup(**self._init_config)

    def _save_init_config_parameters(self, **kwargs) -> None:
        """
            Should only be called by the init method during object creation.
        """
        self._init_config.update(kwargs)

    @ abstractmethod
    def _set_init_default_parameters(self) -> None:
        pass

    @ abstractmethod
    def execute(self, input_hdu: fits.HDUList) -> fits.HDUList:
        """
            Execute the module and update the fits_file.
        """
