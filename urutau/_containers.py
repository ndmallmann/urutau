from typing import Type

from .modules._module_base import AbstractModule


class BaseContainer:
    def __init__(self, config: dict) -> None:
        self._config: dict = config

    @property
    def config(self) -> dict:
        return self._config.copy()


class ModuleContainer(BaseContainer):

    def __init__(self, module: Type[AbstractModule], config: dict) -> None:
        self._module: AbstractModule = module
        BaseContainer.__init__(self, config)

    @property
    def module(self) -> Type[AbstractModule]:
        return self._module

    def get_loaded_module(self, target_config: dict = None) -> AbstractModule:
        if target_config is None:
            target_config = dict()

        config_copy = self._config.copy()
        config_copy.update(target_config)

        loaded_module = self._module(**config_copy)

        return loaded_module


class TargetContainer(BaseContainer):

    def __init__(self, target: str, config: dict = None) -> None:
        self._target = target
        BaseContainer.__init__(self, dict() if config is None else config)

    @property
    def target(self) -> str:
        return self._target
