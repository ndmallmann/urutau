"""
base_grid_utils.py — Cross-checks the Starlight population bins against the
actual stellar population base used by the default grid file.

The grid file (e.g. reference_grid_muse_newMiles.in) points to a base file
name on the data line right after its header (3rd field, e.g.
"BaseM23BI130SY"). That base file lives next to the Starlight executable
("execute path") and lists, one component per line, the spectrum file name,
its age (yr) and metallicity.

Any base component whose name hints at an AGN component (contains "agn",
case-insensitive) must not fall inside a stellar Population Bin — those
components have their own dedicated fields (Featureless-Continuum
exponents / Black-Body temperatures).
"""

import os
import re

from urutau.starlight_utils import GridReader

AGN_HINT = "agn"

_RE_SPACES = re.compile(r"\s+")
_RE_BRACKETS = re.compile(r"\[.*?\]")


class BaseGridError(Exception):
    """Raised when the grid file or base file can't be resolved/parsed."""


def resolve_base_path(grid_file: str, starlight_path: str) -> str:
    """
    Returns the path to the base file referenced by grid_file, resolved
    relative to the Starlight executable's directory (the "execute path"),
    mirroring how urutau.starlight_utils.StarlightWrapper resolves it.
    """
    if not grid_file or not os.path.isfile(grid_file):
        raise BaseGridError(f"Grid file not found: '{grid_file}'")
    if not starlight_path:
        raise BaseGridError("Starlight executable path is required to locate the base file.")

    grid = GridReader(grid_file)
    entries = grid.input_entries
    if not entries:
        raise BaseGridError(f"Grid file '{grid_file}' has no input entries below its header.")

    base_name = entries[0]["base"]
    base_dir = os.path.dirname(starlight_path)
    base_path = os.path.join(base_dir, base_name)

    if not os.path.isfile(base_path):
        raise BaseGridError(
            f"Base file '{base_name}' (from the grid file) was not found next to the "
            f"Starlight executable, expected at: '{base_path}'"
        )
    return base_path


def read_base_components(grid_file: str, starlight_path: str) -> list:
    """
    Returns a list of {"name", "age", "metal", "is_agn"} dicts, one per
    stellar population component listed in the base file.
    """
    base_path = resolve_base_path(grid_file, starlight_path)

    with open(base_path, "r", encoding="utf-8") as base_handler:
        lines = base_handler.readlines()

    if not lines:
        raise BaseGridError(f"Base file '{base_path}' is empty.")

    number_value = int(_RE_BRACKETS.sub("", lines[0]).strip())

    components = []
    for line in lines[1:1 + number_value]:
        parts = _RE_SPACES.split(line.strip())
        if len(parts) < 2:
            continue
        name = parts[0]
        age = float(parts[1])
        metal = float(parts[2]) if len(parts) > 2 else None
        components.append({
            "name": name,
            "age": age,
            "metal": metal,
            "is_agn": AGN_HINT in name.lower(),
        })
    return components


def find_agn_bin_conflicts(components: list, population_bins: dict) -> list:
    """
    Returns a list of {"bin", "component", "age"} conflicts: AGN-hinted base
    components whose age falls inside one of the given population bins
    ({name: (min, max)}), where they should NOT be, since AGN components use
    their own dedicated fields (fc exps / bb temps) instead.

    Matches Urutau's own bin membership test (see
    urutau/starlight_utils/_starlight_wrapper.py, e.g. _pop_by_light /
    _pop_by_mass): the minimum is EXCLUSIVE and the maximum is INCLUSIVE,
    i.e. a component belongs to a bin when min < age <= max.
    """
    conflicts = []
    for bin_name, (min_age, max_age) in population_bins.items():
        lo, hi = (min_age, max_age) if min_age <= max_age else (max_age, min_age)
        for comp in components:
            if comp["is_agn"] and lo < comp["age"] <= hi:
                conflicts.append({
                    "bin": bin_name, "component": comp["name"], "age": comp["age"],
                })
    return conflicts
