"""One- and two-compartment pharmacokinetics with oral or IV bolus dosing and a
direct-effect PD model.

Units are fixed across the package: dose in mg, volume in L, clearance in L/h,
rates in 1/h, time in h, concentration in mg/L.

States: the amount in the gut (A_gut), in the central compartment (A_c) and, for
the two-compartment model, in the peripheral compartment (A_p). Plasma
concentration is C = A_c / V, with V the central volume (V1).

    dA_gut/dt = -ka * A_gut
    dA_c/dt   =  ka * A_gut - elimination - distribution
    dA_p/dt   =  distribution                        (two-compartment only)
    distribution = Q * (A_c / V - A_p / V2)
    elimination = (CL / V) * A_c                     (linear)
                = Vmax * C / (Km + C)                (Michaelis-Menten)

Q is the intercompartmental clearance (L/h) and V2 the peripheral volume (L).
Elimination always happens from the central compartment.

Dosing route:
    oral      each dose adds bioavailability * dose to A_gut, then absorbs with ka.
    iv_bolus  each dose is injected straight into the central compartment (the
              whole dose; ka and bioavailability play no part). Plasma
              concentration then jumps at a dose time, and the value reported at
              a dose time includes that dose.

Effect (direct, no delay), a sigmoid Emax model of the plasma concentration:

    E = E0 + Emax * C**h / (EC50**h + C**h)

Use a negative Emax for an inhibitory effect.

Oral data alone cannot separate bioavailability F from the volumes and clearances,
so with F = 1 a fitted V is really V/F (and CL, V2, Q are CL/F, V2/F, Q/F).
IV bolus data do not have that problem.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.integrate import solve_ivp

_trapezoid = getattr(np, "trapezoid", None) or np.trapz

ROUTES = ("oral", "iv_bolus")


@dataclass(frozen=True)
class PKParameters:
    ka: float  # absorption rate constant, 1/h (ignored for IV bolus dosing)
    volume: float  # central volume of distribution V (V1), L
    elimination: Literal["linear", "michaelis_menten"] = "linear"
    clearance: float | None = None  # L/h, linear elimination
    vmax: float | None = None  # mg/h, Michaelis-Menten elimination
    km: float | None = None  # mg/L, Michaelis-Menten elimination
    bioavailability: float = 1.0  # fraction of each oral dose that reaches the gut
    peripheral_volume: float | None = None  # V2, L (two-compartment model)
    intercompartmental_clearance: float | None = None  # Q, L/h (two-compartment model)

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
        v2, q = self.peripheral_volume, self.intercompartmental_clearance
        if (v2 is None) != (q is None):
            raise ValueError(
                "a two-compartment model needs both the peripheral volume and the "
                "intercompartmental clearance"
            )
        if v2 is not None and (v2 <= 0 or q <= 0):
            raise ValueError("peripheral volume and intercompartmental clearance must be positive")

    @property
    def compartments(self) -> int:
        return 2 if self.peripheral_volume is not None else 1

    @property
    def ke(self) -> float:
        """Elimination rate constant CL/V (1/h). Only defined for linear elimination."""
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
    route: Literal["oral", "iv_bolus"] = "oral"

    def __post_init__(self) -> None:
        if self.dose <= 0:
            raise ValueError("dose must be positive")
        if self.n_doses < 1:
            raise ValueError("n_doses must be at least 1")
        if self.n_doses > 1 and self.interval <= 0:
            raise ValueError("interval must be positive when n_doses > 1")
        if self.start < 0:
            raise ValueError("start must not be negative")
        if self.route not in ROUTES:
            raise ValueError(f"route must be one of {ROUTES}")

    def dose_times(self) -> np.ndarray:
        return self.start + self.interval * np.arange(self.n_doses)


@dataclass
class Simulation:
    t: np.ndarray
    gut: np.ndarray  # amount in gut, mg (post-dose at a dose time; zero for IV bolus)
    central: np.ndarray  # amount in central compartment, mg
    conc: np.ndarray  # plasma concentration, mg/L
    effect: np.ndarray | None  # None when no effect parameters were given
    peripheral: np.ndarray | None = None  # amount in peripheral compartment, mg (two-compartment)

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


def _elimination_clearance(params: PKParameters) -> float:
    """Clearance at low concentration: CL for linear elimination, Vmax/Km otherwise."""
    if params.elimination == "linear":
        return params.clearance
    return params.vmax / params.km


def hybrid_rate_constants(params: PKParameters) -> tuple[float, float]:
    """(alpha, beta) of a two-compartment model, alpha > beta, in 1/h.

    They are the roots of s**2 - (k10 + k12 + k21) * s + k10 * k21 = 0 with
    k10 = CL/V, k12 = Q/V, k21 = Q/V2. For Michaelis-Menten elimination they
    describe the low-concentration limit (CL = Vmax/Km).
    """
    if params.compartments != 2:
        raise ValueError("hybrid rate constants need a two-compartment model")
    k10 = _elimination_clearance(params) / params.volume
    k12 = params.intercompartmental_clearance / params.volume
    k21 = params.intercompartmental_clearance / params.peripheral_volume
    total = k10 + k12 + k21
    root = math.sqrt(max(total * total - 4.0 * k10 * k21, 0.0))
    return 0.5 * (total + root), 0.5 * (total - root)


def half_lives(params: PKParameters) -> dict[str, float | None]:
    """Terminal half-life, and the distribution half-life for a two-compartment model.

    For Michaelis-Menten elimination the values describe the low-concentration tail.
    """
    if params.compartments == 2:
        alpha, beta = hybrid_rate_constants(params)
        return {"terminal": math.log(2.0) / beta, "distribution": math.log(2.0) / alpha}
    return {
        "terminal": math.log(2.0) * params.volume / _elimination_clearance(params),
        "distribution": None,
    }


def _rhs(params: PKParameters):
    ka, v = params.ka, params.volume
    if params.compartments == 2:
        v2, q = params.peripheral_volume, params.intercompartmental_clearance
        if params.elimination == "linear":
            cl = params.clearance

            def rhs(_t, y):
                c = y[1] / v
                flow = q * (c - y[2] / v2)
                return [-ka * y[0], ka * y[0] - cl * c - flow, flow]

        else:
            vmax, km = params.vmax, params.km

            def rhs(_t, y):
                c = max(y[1], 0.0) / v
                flow = q * (c - y[2] / v2)
                return [-ka * y[0], ka * y[0] - vmax * c / (km + c) - flow, flow]

    elif params.elimination == "linear":
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
    additions (to the gut for oral dosing, to the central compartment for an IV
    bolus); at a dose time the reported state includes that dose. With oral dosing
    plasma concentration is continuous across doses, with an IV bolus it jumps.
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
    iv = regimen.route == "iv_bolus"
    amount = regimen.dose if iv else params.bioavailability * regimen.dose
    target = 1 if iv else 0  # state that receives the dose
    rhs = _rhs(params)

    n_states = 3 if params.compartments == 2 else 2
    y = np.zeros(n_states)
    out = np.zeros((n_states, ts.size))
    dose_set = set(float(d) for d in dose_times)

    if t_end == 0.0:
        if 0.0 in dose_set:
            y[target] += amount
        out[:, :] = y[:, None]
    else:
        for k in range(len(breaks) - 1):
            t0, t1 = float(breaks[k]), float(breaks[k + 1])
            if t0 in dose_set:
                y[target] += amount
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
        if t_end in dose_set:
            # a dose exactly at the end of the window: no segment starts there, so apply it
            # to the points at that time here (the reported state includes the dose)
            final = y.copy()
            final[target] += amount
            out[:, ts >= t_end] = final[:, None]

    inv = np.empty_like(order)
    inv[order] = np.arange(order.size)
    gut = out[0][inv]
    central = out[1][inv]
    peripheral = out[2][inv] if n_states == 3 else None
    conc = central / params.volume
    eff = effect_from_conc(conc, effect) if effect is not None else None
    return Simulation(
        t=times, gut=gut, central=central, conc=conc, effect=eff, peripheral=peripheral
    )
