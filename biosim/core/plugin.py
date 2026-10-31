"""Plugin interface and registry.

Every wrapped scientific tool implements Plugin. Workflows only talk to the
registry, so adding or swapping a tool never touches the core.
"""
from __future__ import annotations

import shutil
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any


class Plugin(ABC):
    name: str = ""
    description: str = ""
    required_executables: tuple[str, ...] = ()

    def missing_executables(self) -> list[str]:
        return [exe for exe in self.required_executables if shutil.which(exe) is None]

    @abstractmethod
    def run(self, inputs: dict[str, Any], workdir: Path) -> dict[str, Any]:
        """Run the tool. Take a dict of inputs, return a dict of outputs."""


_REGISTRY: dict[str, Plugin] = {}


def register(cls: type[Plugin]) -> type[Plugin]:
    instance = cls()
    if not instance.name:
        raise ValueError(f"{cls.__name__} must define a non-empty 'name'")
    _REGISTRY[instance.name] = instance
    return cls


def get_plugin(name: str) -> Plugin:
    if name not in _REGISTRY:
        raise KeyError(f"Unknown plugin '{name}'. Available: {sorted(_REGISTRY)}")
    plugin = _REGISTRY[name]
    missing = plugin.missing_executables()
    if missing:
        raise RuntimeError(
            f"Plugin '{name}' needs executables not found on PATH: {', '.join(missing)}"
        )
    return plugin


def list_plugins() -> list[tuple[str, str, list[str]]]:
    return [(p.name, p.description, p.missing_executables()) for p in _REGISTRY.values()]
