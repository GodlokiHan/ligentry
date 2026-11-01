"""Load measured pharmacodynamic effects from CSV-like text (no Qt, easy to test).

Expected content: one observation per line, time in hours and the measured
effect in whatever unit you use (the unit is only a label; it is not converted).

    time,effect
    0.5,12.0
    1,31.5

Accepted variations match biosim.pkpd.observations:
- a header line (optional). When present, columns are found by name ("time",
  "t", "hours" / "effect", "response", "pd", ...), in any order, with extra
  columns ignored. A header with two unrecognised names is read as time, effect.
- comma, semicolon, tab or whitespace separated. With semicolons or tabs a comma
  is read as the decimal separator (1,5 means 1.5).
- blank lines and lines starting with # are skipped.

Unlike concentrations, an effect may be negative (a drop below baseline).

Without a header the file must have exactly two columns (time, effect). A
headerless file with more columns is refused, because it is impossible to tell
which column is the effect, and a measured concentration column read as the
effect would give a meaningless fit.

Rows are returned sorted by time. Every problem is reported with its line number.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# Shared with the concentration loader so both accept the same file layouts.
from biosim.pkpd.observations import (
    ObservationError,
    _detect_delimiter,
    _is_time_name,
    _norm,
    _split,
    _to_float,
)


@dataclass(frozen=True)
class EffectObservations:
    times: np.ndarray  # hours, ascending
    effect: np.ndarray  # measured effect, in the user's own unit
    source: str = ""

    def __len__(self) -> int:
        return int(self.times.size)


def _is_effect_name(name: str) -> bool:
    n = _norm(name)
    return n.startswith(("effect", "response", "resp")) or n in {"e", "pd", "y"}


def parse_effect_observations(text: str, source: str = "") -> EffectObservations:
    lines = [
        (no, ln)
        for no, ln in enumerate(text.splitlines(), start=1)
        if ln.strip() and not ln.lstrip().startswith("#")
    ]
    if not lines:
        raise ObservationError("No data found: the input is empty.")

    delimiter = _detect_delimiter(lines[0][1])
    decimal_comma = delimiter in (";", "\t")
    rows = [(no, _split(ln, delimiter)) for no, ln in lines]

    t_idx, e_idx = 0, 1
    first_no, first_cells = rows[0]
    has_header = all(_to_float(c, decimal_comma) is None for c in first_cells)
    if has_header:
        t_found = [i for i, n in enumerate(first_cells) if _is_time_name(n)]
        e_found = [i for i, n in enumerate(first_cells) if _is_effect_name(n)]
        if t_found and e_found:
            t_idx, e_idx = t_found[0], e_found[0]
        elif len(first_cells) != 2:
            raise ObservationError(
                f"Line {first_no}: cannot tell which columns are time and effect "
                f"(header: {first_cells}). Name them 'time' and 'effect'."
            )
        rows = rows[1:]
        if not rows:
            raise ObservationError("No data rows found: the file only has a header.")

    need = max(t_idx, e_idx) + 1
    times: list[float] = []
    effects: list[float] = []
    for no, cells in rows:
        if len(cells) < need:
            raise ObservationError(
                f"Line {no}: expected at least {need} columns, found {len(cells)}."
            )
        if not has_header and len(cells) != 2:
            raise ObservationError(
                f"Line {no}: found {len(cells)} columns but there is no header. Add a header "
                "line naming the 'time' and 'effect' columns, or keep only those two columns."
            )
        t = _to_float(cells[t_idx], decimal_comma)
        e = _to_float(cells[e_idx], decimal_comma)
        if t is None:
            raise ObservationError(f"Line {no}: time '{cells[t_idx]}' is not a number.")
        if e is None:
            raise ObservationError(f"Line {no}: effect '{cells[e_idx]}' is not a number.")
        if not (math.isfinite(t) and math.isfinite(e)):
            raise ObservationError(f"Line {no}: values must be finite numbers.")
        if t < 0:
            raise ObservationError(f"Line {no}: time must not be negative (got {t:g}).")
        times.append(t)
        effects.append(e)

    order = np.argsort(np.asarray(times), kind="stable")
    return EffectObservations(
        times=np.asarray(times, dtype=float)[order],
        effect=np.asarray(effects, dtype=float)[order],
        source=source,
    )


def load_effect_observations(path: str | Path) -> EffectObservations:
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8-sig")  # tolerate a BOM from Excel exports
    except (OSError, UnicodeDecodeError) as exc:
        raise ObservationError(f"Cannot read {path}: {exc}") from exc
    return parse_effect_observations(text, source=str(path))
