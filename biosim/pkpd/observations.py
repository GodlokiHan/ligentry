"""Load measured plasma concentrations from CSV-like text (no Qt, easy to test).

Expected content: one observation per line, time in hours and concentration in
mg/L, matching the units used everywhere in biosim.pkpd.

    time,conc
    0.5,1.8
    1,3.1

Accepted variations:
- a header line (optional). When present, columns are found by name ("time",
  "t", "hours" / "conc", "concentration", "cp", ...), in any order, with extra
  columns ignored. A header with two unrecognised names is read as time, conc.
- without a header, the first column is time and the second is concentration.
- comma, semicolon, tab or whitespace separated. With semicolons or tabs a comma
  is read as the decimal separator (1,5 means 1.5).
- blank lines and lines starting with # are skipped.

Rows are returned sorted by time. Every problem is reported with its line number.
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


class ObservationError(ValueError):
    """The data could not be read as time/concentration observations."""


@dataclass(frozen=True)
class Observations:
    times: np.ndarray  # hours, ascending
    conc: np.ndarray  # mg/L
    source: str = ""

    def __len__(self) -> int:
        return int(self.times.size)


def _norm(name: str) -> str:
    return "".join(name.lower().split())


def _is_time_name(name: str) -> bool:
    n = _norm(name)
    return n.startswith("time") or n in {"t", "h", "hr", "hrs", "hour", "hours"}


def _is_conc_name(name: str) -> bool:
    n = _norm(name)
    return n.startswith("conc") or n in {"c", "cp", "cc"}


def _detect_delimiter(line: str) -> str | None:
    """The separator used in `line`, or None for whitespace-separated."""
    counts = {d: line.count(d) for d in (",", ";", "\t")}
    best = max(counts, key=counts.get)
    return best if counts[best] > 0 else None


def _split(line: str, delimiter: str | None) -> list[str]:
    if delimiter is None:
        return line.split()
    return [c.strip() for c in next(csv.reader([line], delimiter=delimiter))]


def _to_float(cell: str, decimal_comma: bool) -> float | None:
    text = cell.strip()
    if decimal_comma:
        text = text.replace(",", ".")
    try:
        return float(text)
    except ValueError:
        return None


def parse_observations(text: str, source: str = "") -> Observations:
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

    t_idx, c_idx = 0, 1
    first_no, first_cells = rows[0]
    has_header = all(_to_float(c, decimal_comma) is None for c in first_cells)
    if has_header:
        t_found = [i for i, n in enumerate(first_cells) if _is_time_name(n)]
        c_found = [i for i, n in enumerate(first_cells) if _is_conc_name(n)]
        if t_found and c_found:
            t_idx, c_idx = t_found[0], c_found[0]
        elif len(first_cells) != 2:
            raise ObservationError(
                f"Line {first_no}: cannot tell which columns are time and concentration "
                f"(header: {first_cells}). Name them 'time' and 'conc'."
            )
        rows = rows[1:]
        if not rows:
            raise ObservationError("No data rows found: the file only has a header.")

    need = max(t_idx, c_idx) + 1
    times: list[float] = []
    concs: list[float] = []
    for no, cells in rows:
        if len(cells) < need:
            raise ObservationError(
                f"Line {no}: expected at least {need} columns, found {len(cells)}."
            )
        t = _to_float(cells[t_idx], decimal_comma)
        c = _to_float(cells[c_idx], decimal_comma)
        if t is None:
            raise ObservationError(f"Line {no}: time '{cells[t_idx]}' is not a number.")
        if c is None:
            raise ObservationError(f"Line {no}: concentration '{cells[c_idx]}' is not a number.")
        if not (math.isfinite(t) and math.isfinite(c)):
            raise ObservationError(f"Line {no}: values must be finite numbers.")
        if t < 0:
            raise ObservationError(f"Line {no}: time must not be negative (got {t:g}).")
        if c < 0:
            raise ObservationError(f"Line {no}: concentration must not be negative (got {c:g}).")
        times.append(t)
        concs.append(c)

    order = np.argsort(np.asarray(times), kind="stable")
    return Observations(
        times=np.asarray(times, dtype=float)[order],
        conc=np.asarray(concs, dtype=float)[order],
        source=source,
    )


def load_observations(path: str | Path) -> Observations:
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8-sig")  # tolerate a BOM from Excel exports
    except (OSError, UnicodeDecodeError) as exc:
        raise ObservationError(f"Cannot read {path}: {exc}") from exc
    return parse_observations(text, source=str(path))
