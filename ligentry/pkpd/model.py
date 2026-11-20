"""One- and two-compartment pharmacokinetics with oral, IV bolus, or IV infusion dosing
and a sigmoid Emax PD model (direct, or delayed through an effect compartment).

Units are fixed across the package: dose in mg, volume in L, clearance in L/h,
rates in 1/h, time in h, concentration in mg/L.

States: the amount in the gut (A_gut), in the central compartment (A_c) and, for
the two-compartment model, in the peripheral compartment (A_p). Plasma
concentration is C = A_c / V, with V the central volume (V1).

    dA_gut/dt = -ka * A_gut
    dA_c/dt   =  ka * A_gut + R_inf(t) - elimination - distribution
    dA_p/dt   =  distribution                        (two-compartment only)
    distribution = Q * (A_c / V - A_p / V2)
    elimination = (CL / V) * A_c                     (linear)
                = Vmax * C / (Km + C)                (Michaelis-Menten)

Q is the intercompartmental clearance (L/h) and V2 the peripheral volume (L).
Elimination always happens from the central compartment.

Dosing routes:
    oral         each dose adds bioavailability * dose to A_gut at dose_time + tlag,
                 then absorbs with rate ka.
    iv_bolus     each dose is injected straight into the central compartment (the
                 whole dose; ka, tlag and bioavailability play no part). Plasma
                 concentration jumps at a dose time, and the value reported at
                 a dose time includes that dose.
    iv_infusion  each dose is infused into the central compartment at constant zero-order
                 rate dose / infusion_duration over [t_dose, t_dose + infusion_duration).
                 Overlapping infusions sum their rates.

Effect, a sigmoid Emax model:

    E = E0 + Emax * Ce**h / (EC50**h + Ce**h)

By default the effect is direct: Ce is the plasma concentration C and there is no
delay. Give a rate constant ke0 (1/h) to add an effect compartment instead. Ce is
then the concentration at a hypothetical effect site that follows plasma with a
lag:

    dCe/dt = ke0 * (C - Ce)
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.integrate import solve_ivp

_trapezoid = getattr(np, "trapezoid", None) or np.trapz

ROUTES = ("oral", "iv_bolus", "iv_infusion")


@dataclass(frozen=True)
class PKParameters:
    ka: float  # absorption rate constant, 1/h (ignored for IV routes)
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
    ke0: float | None = None  # effect-compartment rate constant, 1/h (None = direct effect)

    def __post_init__(self) -> None:
        if self.ec50 <= 0 or self.hill <= 0:
            raise ValueError("ec50 and hill must be positive")
        if self.ke0 is not None and not self.ke0 > 0:
            raise ValueError("ke0 must be positive (or None for a direct effect)")

    @property
    def delayed(self) -> bool:
        return self.ke0 is not None


@dataclass(frozen=True)
class Regimen:
    dose: float  # mg per dose
    interval: float = 0.0  # h between doses (needed when n_doses > 1)
    n_doses: int = 1
    start: float = 0.0  # time of the first dose, h
    route: Literal["oral", "iv_bolus", "iv_infusion"] = "oral"
    infusion_duration: float = 1.0  # h per dose (used when route == "iv_infusion")
    tlag: float = 0.0  # absorption lag time, h (used when route == "oral")

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
        if self.route == "iv_infusion" and (not math.isfinite(self.infusion_duration) or self.infusion_duration <= 0):
            raise ValueError("infusion_duration must be positive for iv_infusion")
        if not math.isfinite(self.tlag) or self.tlag < 0:
            raise ValueError("tlag must not be negative")

    def dose_times(self) -> np.ndarray:
        return self.start + self.interval * np.arange(self.n_doses)


@dataclass
class Simulation:
    t: np.ndarray
    gut: np.ndarray  # amount in gut, mg (post-dose at an oral absorption time; zero for IV routes)
    central: np.ndarray  # amount in central compartment, mg
    conc: np.ndarray  # plasma concentration, mg/L
    effect: np.ndarray | None  # None when no effect parameters were given
    peripheral: np.ndarray | None = None  # amount in peripheral compartment, mg (two-compartment)
    effect_conc: np.ndarray | None = None  # effect-site concentration Ce, mg/L (effect compartment)

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


def equilibration_half_life(ke0: float) -> float:
    """Half-life (h) with which the effect site follows plasma, ln2/ke0."""
    if not ke0 > 0:
        raise ValueError("ke0 must be positive")
    return math.log(2.0) / ke0


def _elimination_clearance(params: PKParameters) -> float:
    """Clearance at low concentration: CL for linear elimination, Vmax/Km otherwise."""
    if params.elimination == "linear":
        return params.clearance
    return params.vmax / params.km


def hybrid_rate_constants(params: PKParameters) -> tuple[float, float]:
    """(alpha, beta) of a two-compartment model, alpha > beta, in 1/h."""
    if params.compartments != 2:
        raise ValueError("hybrid rate constants need a two-compartment model")
    k10 = _elimination_clearance(params) / params.volume
    k12 = params.intercompartmental_clearance / params.volume
    k21 = params.intercompartmental_clearance / params.peripheral_volume
    total = k10 + k12 + k21
    root = math.sqrt(max(total * total - 4.0 * k10 * k21, 0.0))
    return 0.5 * (total + root), 0.5 * (total - root)


def half_lives(params: PKParameters) -> dict[str, float | None]:
    """Terminal half-life, and the distribution half-life for a two-compartment model."""
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
    """Integrate the model and return the state at the requested times."""
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or times.size == 0:
        raise ValueError("times must be a non-empty 1-D array")
    if np.any(times < 0):
        raise ValueError("times must not be negative")

    order = np.argsort(times, kind="stable")
    ts = times[order]
    t_end = float(ts[-1])

    raw_doses = regimen.dose_times()
    route = regimen.route

    if route == "oral":
        arrival_times = raw_doses + regimen.tlag
        active_arrivals = arrival_times[arrival_times <= t_end]
        breaks = np.unique(np.concatenate(([0.0], active_arrivals, [t_end])))
        amount = params.bioavailability * regimen.dose
        arrival_set = set(float(d) for d in active_arrivals)
    elif route == "iv_bolus":
        active_arrivals = raw_doses[raw_doses <= t_end]
        breaks = np.unique(np.concatenate(([0.0], active_arrivals, [t_end])))
        amount = regimen.dose
        arrival_set = set(float(d) for d in active_arrivals)
    else:  # iv_infusion
        inf_starts = raw_doses[raw_doses <= t_end]
        inf_stops = (raw_doses + regimen.infusion_duration)
        inf_stops = inf_stops[inf_stops <= t_end]
        breaks = np.unique(np.concatenate(([0.0], inf_starts, inf_stops, [t_end])))
        amount = 0.0
        arrival_set = set()
        inf_rate_single = regimen.dose / regimen.infusion_duration
        inf_end_times = raw_doses + regimen.infusion_duration

    pk_rhs = _rhs(params)
    n_pk = 3 if params.compartments == 2 else 2
    ke0 = effect.ke0 if effect is not None else None
    volume = params.volume

    def make_rhs(inf_rate: float):
        if ke0 is None and inf_rate == 0.0:
            return pk_rhs

        def rhs(t, y):
            d = pk_rhs(t, y[:n_pk])
            if inf_rate != 0.0:
                d[1] += inf_rate
            if ke0 is not None:
                d.append(ke0 * (max(y[1], 0.0) / volume - y[n_pk]))
            return d

        return rhs

    n_states = n_pk + (1 if ke0 is not None else 0)
    y = np.zeros(n_states)
    out = np.zeros((n_states, ts.size))
    bolus_target = 1 if route == "iv_bolus" else 0

    if t_end == 0.0:
        if 0.0 in arrival_set:
            y[bolus_target] += amount
        out[:, :] = y[:, None]
    else:
        for k in range(len(breaks) - 1):
            t0, t1 = float(breaks[k]), float(breaks[k + 1])
            if t0 in arrival_set:
                y[bolus_target] += amount
            if route == "iv_infusion":
                t_mid = 0.5 * (t0 + t1)
                n_active = int(np.sum((raw_doses <= t_mid) & (t_mid < inf_end_times)))
                seg_rhs = make_rhs(n_active * inf_rate_single)
            else:
                seg_rhs = make_rhs(0.0)
            last = k == len(breaks) - 2
            mask = (ts >= t0) & ((ts <= t1) if last else (ts < t1))
            seg_times = ts[mask]
            eval_pts = np.unique(np.concatenate((seg_times, [t1])))
            sol = solve_ivp(
                seg_rhs, (t0, t1), y, t_eval=eval_pts, method="LSODA", rtol=1e-9, atol=1e-12
            )
            if not sol.success:
                raise RuntimeError(f"ODE solver failed: {sol.message}")
            idx = np.searchsorted(eval_pts, seg_times)
            out[:, mask] = sol.y[:, idx]
            y = sol.y[:, -1].copy()
        if t_end in arrival_set:
            final = y.copy()
            final[bolus_target] += amount
            out[:, ts >= t_end] = final[:, None]

    inv = np.empty_like(order)
    inv[order] = np.arange(order.size)
    gut = out[0][inv]
    central = out[1][inv]
    peripheral = out[2][inv] if n_pk == 3 else None
    conc = central / params.volume
    site = np.clip(out[n_pk][inv], 0.0, None) if ke0 is not None else None
    eff = None
    if effect is not None:
        eff = effect_from_conc(site if site is not None else conc, effect)
    return Simulation(
        t=times, gut=gut, central=central, conc=conc, effect=eff, peripheral=peripheral,
        effect_conc=site,
    )


def effect_site_concentration(
    params: PKParameters, regimen: Regimen, times: np.ndarray, ke0: float
) -> np.ndarray:
    """Effect-site concentration Ce (mg/L) at `times` for an effect compartment with rate ke0."""
    sim = simulate(params, regimen, times, EffectParameters(ke0=ke0))
    return sim.effect_conc
