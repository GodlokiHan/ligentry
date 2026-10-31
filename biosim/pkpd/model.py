"""One-compartment pharmacokinetics with oral dosing and a direct-effect PD model.

Units are fixed across the package: dose in mg, volume in L, clearance in L/h,
rates in 1/h, time in h, concentration in mg/L.

States: the amount in the gut (A_gut) and the amount in the central compartment
(A_c). Plasma concentration is C = A_c / V.

    dA_gut/dt = -ka * A_gut
    dA_c/dt   =  ka * A_gut - elimination
    elimination = (CL / V) * A_c                    (linear)
                = Vmax * C / (Km + C)               (Michaelis-Menten)

Each dose adds bioavailability * dose to A_gut at its dose time.

Effect (direct, no delay), a sigmoid Emax model of the plasma concentration:

    E = E0 + Emax * C**h / (EC50**h + C**h)

Use a negative Emax for an inhibitory effect.

Oral data alone cannot separate bioavailability F from volume V, so with F = 1
a fitted V is really V/F.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.integrate import solve_ivp

_trapezoid = getattr(np, "trapezoid", None) or np.trapz


@dataclass(frozen=True)
class PKParameters:
    ka: float  # absorption rate constant, 1/h
    volume: float  # volume of distribution, L
    elimination: Literal["linear", "michaelis_menten"] = "linear"
    clearance: float | None = None  # L/h, linear elimination
    vmax: float | None = None  # mg/h, Michaelis-Menten elimination
    km: float | None = None  # mg/L, Michaelis-Menten elimination
    bioavailability: float = 1.0  # fraction of each dose that reaches the gut

    def __post_init__(self) -> None:
        if self.ka <= 0 or self.volume <= 0:
            raise ValueError("ka and volume must be positive")
        if not 0 < self.bioavailability <= 1:
            raise ValueError("bioavailability must be in (0, 1]")
        if self.elimination == "linear":
            if self.clearance is None or self.clearance <= 0:
                raise ValueError("linear elimination needs a positive clearance")
        elif self.elimination == "michaelis_menten":
            if self.vmax is None or self.km is None or self.vmax <= 0 or self.km <= 0:
                raise ValueError("Michaelis-Menten elimination needs positive vmax and km")
        else:
            raise ValueError("elimination must be 'linear' or 'michaelis_menten'")

    @property
    def ke(self) -> float:
        """Elimination rate constant (1/h). Only defined for linear elimination."""
        if self.elimination != "linear":
            raise ValueError("ke is only defined for linear elimination")
        return self.clearance / self.volume


@dataclass(frozen=True)
class EffectParameters:
    e0: float = 0.0  # baseline effect
    emax: float = 1.0  # maximum change from baseline (negative for inhibition)
    ec50: float = 1.0  # concentration giving half of emax, mg/L
    hill: float = 1.0  # sigmoidicity

    def __post_init__(self) -> None:
        if self.ec50 <= 0 or self.hill <= 0:
            raise ValueError("ec50 and hill must be positive")


@dataclass(frozen=True)
class Regimen:
    dose: float  # mg per dose
    interval: float = 0.0  # h between doses (needed when n_doses > 1)
    n_doses: int = 1
    start: float = 0.0  # time of the first dose, h

    def __post_init__(self) -> None:
        if self.dose <= 0:
            raise ValueError("dose must be positive")
        if self.n_doses < 1:
            raise ValueError("n_doses must be at least 1")
        if self.n_doses > 1 and self.interval <= 0:
            raise ValueError("interval must be positive when n_doses > 1")
        if self.start < 0:
            raise ValueError("start must not be negative")

    def dose_times(self) -> np.ndarray:
        return self.start + self.interval * np.arange(self.n_doses)


@dataclass
class Simulation:
    t: np.ndarray
    gut: np.ndarray  # amount in gut, mg (post-dose at a dose time)
    central: np.ndarray  # amount in central compartment, mg
    conc: np.ndarray  # plasma concentration, mg/L
    effect: np.ndarray | None  # None when no effect parameters were given

    def metrics(self) -> dict[str, float]:
        """Cmax, Tmax and AUC over the simulated window (trapezoid rule)."""
        i = int(np.argmax(self.conc))
        return {
            "cmax": float(self.conc[i]),
            "tmax": float(self.t[i]),
            "auc": float(_trapezoid(self.conc, self.t)),
        }


def effect_from_conc(conc: np.ndarray | float, eff: EffectParameters) -> np.ndarray:
    c = np.clip(np.asarray(conc, dtype=float), 0.0, None)
    ch = c**eff.hill
    return eff.e0 + eff.emax * ch / (eff.ec50**eff.hill + ch)


def _rhs(params: PKParameters):
    ka, v = params.ka, params.volume
    if params.elimination == "linear":
        ke = params.clearance / v

        def rhs(_t, y):
            return [-ka * y[0], ka * y[0] - ke * y[1]]

    else:
        vmax, km = params.vmax, params.km

        def rhs(_t, y):
            c = max(y[1], 0.0) / v
            return [-ka * y[0], ka * y[0] - vmax * c / (km + c)]

    return rhs


def simulate(
    params: PKParameters,
    regimen: Regimen,
    times: np.ndarray,
    effect: EffectParameters | None = None,
) -> Simulation:
    """Integrate the model and return the state at the requested times.

    `times` (hours, any order) are evaluated as given; the integration always
    starts at time 0 with an empty body. Doses are applied as instantaneous
    additions to the gut; at a dose time the reported state includes that dose.
    Plasma concentration is continuous across doses.
    """
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or times.size == 0:
        raise ValueError("times must be a non-empty 1-D array")
    if np.any(times < 0):
        raise ValueError("times must not be negative")

    order = np.argsort(times, kind="stable")
    ts = times[order]
    t_end = float(ts[-1])

    dose_times = regimen.dose_times()
    dose_times = dose_times[dose_times <= t_end]
    breaks = np.unique(np.concatenate(([0.0], dose_times, [t_end])))
    amount = params.bioavailability * regimen.dose
    rhs = _rhs(params)

    y = np.zeros(2)
    out = np.zeros((2, ts.size))
    dose_set = set(float(d) for d in dose_times)

    if t_end == 0.0:
        if 0.0 in dose_set:
            y[0] += amount
        out[:, :] = y[:, None]
    else:
        for k in range(len(breaks) - 1):
            t0, t1 = float(breaks[k]), float(breaks[k + 1])
            if t0 in dose_set:
                y[0] += amount
            last = k == len(breaks) - 2
            mask = (ts >= t0) & ((ts <= t1) if last else (ts < t1))
            seg_times = ts[mask]
            eval_pts = np.unique(np.concatenate((seg_times, [t1])))
            sol = solve_ivp(
                rhs, (t0, t1), y, t_eval=eval_pts, method="LSODA", rtol=1e-9, atol=1e-12
            )
            if not sol.success:
                raise RuntimeError(f"ODE solver failed: {sol.message}")
            idx = np.searchsorted(eval_pts, seg_times)
            out[:, mask] = sol.y[:, idx]
            y = sol.y[:, -1].copy()

    inv = np.empty_like(order)
    inv[order] = np.arange(order.size)
    gut = out[0][inv]
    central = out[1][inv]
    conc = central / params.volume
    eff = effect_from_conc(conc, effect) if effect is not None else None
    return Simulation(t=times, gut=gut, central=central, conc=conc, effect=eff)
