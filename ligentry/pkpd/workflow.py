"""Glue between plain form values and the PK/PD model (no Qt, easy to test).

The PK/PD window only gathers numbers from its widgets into a ModelInputs, calls
the functions here and draws what comes back. Units are fixed across ligentry.pkpd:
dose mg, volume L, clearance L/h, rates 1/h, time h, concentration mg/L.
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from ligentry.pkpd.fit import FitResult, fit_concentration
from ligentry.pkpd.model import (
    EffectParameters,
    PKParameters,
    Regimen,
    Simulation,
    half_lives,
    hybrid_rate_constants,
    simulate,
)
from ligentry.pkpd.observations import Observations

LINEAR = "linear"
MICHAELIS_MENTEN = "michaelis_menten"
ORAL = "oral"
IV_BOLUS = "iv_bolus"

# parameter name -> (short label, unit)
PARAMETER_INFO = {
    "ka": ("ka", "1/h"),
    "volume": ("V", "L"),
    "clearance": ("CL", "L/h"),
    "vmax": ("Vmax", "mg/h"),
    "km": ("Km", "mg/L"),
    "peripheral_volume": ("V2", "L"),
    "intercompartmental_clearance": ("Q", "L/h"),
}

MAX_TIME_H = 100_000.0

# Below this ratio the two exponential phases of a two-compartment model are too
# close together for data to tell them apart from a single phase.
_MIN_PHASE_RATIO = 3.0


@dataclass
class ModelInputs:
    # regimen
    dose: float = 100.0
    interval: float = 12.0
    n_doses: int = 1
    # pharmacokinetics
    ka: float = 1.0
    volume: float = 40.0
    elimination: str = LINEAR
    clearance: float = 5.0
    vmax: float = 20.0
    km: float = 2.0
    bioavailability: float = 1.0
    # pharmacodynamics
    use_effect: bool = False
    e0: float = 0.0
    emax: float = 100.0
    ec50: float = 1.0
    hill: float = 1.0
    # time window; None = chosen automatically
    t_end: float | None = None
    # dosing route and the optional peripheral compartment (used when compartments == 2)
    route: str = ORAL
    compartments: int = 1
    peripheral_volume: float = 60.0
    intercompartmental_clearance: float = 8.0


@dataclass
class SimulationResult:
    t_end: float
    pk: PKParameters
    regimen: Regimen
    effect: EffectParameters | None
    sim: Simulation
    metrics: dict[str, float]


def free_choices(
    elimination: str, compartments: int = 1, route: str = ORAL
) -> tuple[str, ...]:
    """Parameters that can be fitted for the given model structure.

    ka is left out for IV bolus dosing; the peripheral volume and intercompartmental
    clearance are added for a two-compartment model.
    """
    if elimination == LINEAR:
        names: tuple[str, ...] = ("ka", "volume", "clearance")
    elif elimination == MICHAELIS_MENTEN:
        names = ("ka", "volume", "vmax", "km")
    else:
        raise ValueError(f"unknown elimination model '{elimination}'")
    if compartments == 2:
        names += ("peripheral_volume", "intercompartmental_clearance")
    elif compartments != 1:
        raise ValueError("compartments must be 1 or 2")
    if route == IV_BOLUS:
        names = tuple(n for n in names if n != "ka")
    elif route != ORAL:
        raise ValueError(f"unknown dosing route '{route}'")
    return names


def build_model(
    inputs: ModelInputs,
) -> tuple[PKParameters, Regimen, EffectParameters | None]:
    """Turn form values into validated model objects (ValueError with a clear message)."""
    if inputs.compartments not in (1, 2):
        raise ValueError("compartments must be 1 or 2")
    two = (
        {
            "peripheral_volume": inputs.peripheral_volume,
            "intercompartmental_clearance": inputs.intercompartmental_clearance,
        }
        if inputs.compartments == 2
        else {}
    )
    if inputs.elimination == LINEAR:
        pk = PKParameters(
            ka=inputs.ka, volume=inputs.volume, clearance=inputs.clearance,
            bioavailability=inputs.bioavailability, **two,
        )
    elif inputs.elimination == MICHAELIS_MENTEN:
        pk = PKParameters(
            ka=inputs.ka, volume=inputs.volume, elimination=MICHAELIS_MENTEN,
            vmax=inputs.vmax, km=inputs.km, bioavailability=inputs.bioavailability, **two,
        )
    else:
        raise ValueError(f"unknown elimination model '{inputs.elimination}'")
    regimen = Regimen(
        dose=inputs.dose, interval=inputs.interval, n_doses=int(inputs.n_doses),
        route=inputs.route,
    )
    effect = (
        EffectParameters(
            e0=inputs.e0, emax=inputs.emax, ec50=inputs.ec50, hill=inputs.hill
        )
        if inputs.use_effect
        else None
    )
    return pk, regimen, effect


def _nice_ceiling(x: float) -> float:
    """Round up to a tidy number of hours (1 h steps for short spans, days for long ones)."""
    if x <= 0:
        return 1.0
    if x < 12:
        step = 1.0
    elif x < 24:
        step = 2.0
    elif x < 72:
        step = 6.0
    elif x < 240:
        step = 12.0
    else:
        step = 24.0
    return math.ceil(x / step) * step


def default_time_span(pk: PKParameters, regimen: Regimen) -> float:
    """A window that shows the last dose being absorbed and mostly eliminated (hours).

    The tail is five terminal half-lives, so a two-compartment model gets the longer
    window its slow terminal phase needs. For Michaelis-Menten elimination the
    zero-order phase at Vmax is added first, then a first-order tail at low
    concentration.
    """
    last_dose = float(regimen.dose_times()[-1])
    absorption = 0.0 if regimen.route == IV_BOLUS else 3.0 / pk.ka
    elimination = 5.0 * half_lives(pk)["terminal"]
    if pk.elimination == MICHAELIS_MENTEN:
        amount = regimen.dose if regimen.route == IV_BOLUS else pk.bioavailability * regimen.dose
        elimination += amount / pk.vmax
    return min(_nice_ceiling(last_dose + absorption + elimination), MAX_TIME_H)


def half_life_text(pk: PKParameters) -> str:
    """One line with the half-life(s) of the model, for the readout under the plot."""
    lives = half_lives(pk)
    tail = " (at low concentration)" if pk.elimination == MICHAELIS_MENTEN else ""
    if lives["distribution"] is None:
        return f"Terminal half-life {lives['terminal']:.3g} h{tail}"
    return (
        f"Half-lives: distribution {lives['distribution']:.3g} h, "
        f"terminal {lives['terminal']:.3g} h{tail}"
    )


def simulate_inputs(
    inputs: ModelInputs, n_points: int = 600, min_t_end: float = 0.0
) -> SimulationResult:
    """Simulate the regimen. `min_t_end` stretches an automatic window to cover data."""
    pk, regimen, effect = build_model(inputs)
    if inputs.t_end is None:
        t_end = default_time_span(pk, regimen)
        if min_t_end > t_end:
            t_end = min(_nice_ceiling(min_t_end), MAX_TIME_H)
    else:
        t_end = float(inputs.t_end)
    if not (math.isfinite(t_end) and 0.0 < t_end <= MAX_TIME_H):
        raise ValueError(f"The end time must be between 0 and {MAX_TIME_H:g} h.")

    dose_times = regimen.dose_times()
    dose_times = dose_times[dose_times <= t_end]
    extra = dose_times
    if regimen.route == IV_BOLUS:
        # a point just before each later dose, so the jump in concentration is drawn as a step
        extra = np.concatenate([dose_times, dose_times[dose_times > 0.0] - 1e-9])
    times = np.unique(np.concatenate([np.linspace(0.0, t_end, int(n_points)), extra]))
    sim = simulate(pk, regimen, times, effect)
    return SimulationResult(
        t_end=t_end, pk=pk, regimen=regimen, effect=effect, sim=sim, metrics=sim.metrics()
    )


def run_fit(
    inputs: ModelInputs,
    obs: Observations,
    free: tuple[str, ...] | list[str],
    weighting: str = "proportional",
) -> FitResult:
    """Fit the chosen parameters to measured concentrations, starting from the form values."""
    pk, regimen, _effect = build_model(inputs)
    return fit_concentration(obs.times, obs.conc, regimen, pk, free=tuple(free), weighting=weighting)


def apply_fit(inputs: ModelInputs, result: FitResult) -> ModelInputs:
    """Copy fitted estimates back into the form values."""
    return replace(inputs, **result.estimates)


def fit_table_rows(result: FitResult) -> list[tuple[str, str, str, str]]:
    """(parameter, estimate, unit, CV %) as display strings, in fit order."""
    rows = []
    for name in result.free:
        label, unit = PARAMETER_INFO[name]
        cv = None if result.cv_percent is None else result.cv_percent.get(name)
        rows.append(
            (label, f"{result.estimates[name]:.4g}", unit, "n/a" if cv is None else f"{cv:.1f}")
        )
    return rows


def fit_notes(result: FitResult) -> list[str]:
    """Plain-language warnings to show next to a fit."""
    notes: list[str] = []
    if not result.success:
        notes.append(f"The optimiser did not converge: {result.message}")
    if result.cv_percent is None:
        notes.append(
            "Too few observations to estimate parameter uncertainty "
            "(more points than fitted parameters are needed)."
        )
    else:
        poor = [PARAMETER_INFO[n][0] for n, cv in result.cv_percent.items() if cv > 50.0]
        if poor:
            notes.append("Poorly determined (CV above 50%): " + ", ".join(poor) + ".")
    p = result.params
    oral = result.route == ORAL
    two = p.compartments == 2
    periph_free = any(n in result.free for n in ("peripheral_volume", "intercompartmental_clearance"))

    if p.elimination == LINEAR and "ka" in result.free:
        if two:
            terminal_rate = hybrid_rate_constants(p)[1]
            what = "the terminal rate constant"
        else:
            terminal_rate = p.ke
            what = "the elimination rate constant CL/V"
        if p.ka < terminal_rate:
            notes.append(
                f"Fitted ka is below {what}. Oral data can swap the "
                "two (flip-flop): check this against what you know about the drug, or fix ka."
            )
    if two and periph_free:
        notes.append(
            "Two-compartment fits need several samples from the early, rapidly falling part "
            "of the curve (the distribution phase). Without them V2 and Q are poorly determined."
        )
        if oral:
            notes.append(
                "With oral dosing the absorption phase hides the distribution phase, so V2 and Q "
                "are determined much better by IV bolus data."
            )
        if p.elimination == LINEAR:
            alpha, beta = hybrid_rate_constants(p)
            if alpha / beta < _MIN_PHASE_RATIO:
                notes.append(
                    f"The two exponential phases are not well separated (alpha/beta = "
                    f"{alpha / beta:.2g}), so these data cannot tell a two-compartment model from "
                    "a one-compartment one. Consider the one-compartment model."
                )
    if oral and any(n in result.free for n in ("volume", "clearance", *(
        ("peripheral_volume", "intercompartmental_clearance") if two else ()
    ))):
        if two:
            notes.append(
                "Oral data cannot separate bioavailability from the volumes and clearances, so "
                "V, CL, V2 and Q are apparent values (V/F, CL/F, V2/F and Q/F)."
            )
        else:
            notes.append(
                "Oral data cannot separate bioavailability from volume and clearance, so V and "
                "CL are apparent values (V/F and CL/F)."
            )
    return notes


def simulation_to_csv(path: str | Path, res: SimulationResult) -> Path:
    """Write the simulated curves to a CSV file (time_h, conc_mg_per_L, ...)."""
    path = Path(path)
    sim = res.sim
    header = ["time_h", "conc_mg_per_L", "gut_mg", "central_mg"]
    if sim.peripheral is not None:
        header.append("peripheral_mg")
    if sim.effect is not None:
        header.append("effect")
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for i in range(sim.t.size):
            row = [sim.t[i], sim.conc[i], sim.gut[i], sim.central[i]]
            if sim.peripheral is not None:
                row.append(sim.peripheral[i])
            if sim.effect is not None:
                row.append(sim.effect[i])
            writer.writerow([f"{v:.8g}" for v in row])
    return path
