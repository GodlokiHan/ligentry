"""Fit the sigmoid Emax effect model to measured effects (no Qt).

The model is the effect relationship from ligentry.pkpd.model:

    E = E0 + Emax * Ce**h / (EC50**h + Ce**h)

Ce is the plasma concentration for a direct effect. With an effect compartment
(ke0 given) Ce is the effect-site concentration, which follows plasma with a lag:

    dCe/dt = ke0 * (C - Ce)

`fit_effect` works on pairs of (plasma concentration, measured effect). The plasma
concentration at each measurement time normally comes from the PK model, see
ligentry.pkpd.effect_workflow. A delayed effect also needs `site`, a function that
returns the effect-site concentration at the measurement times for a given ke0
(the PK model is run again for every ke0 the optimiser tries).

EC50, the Hill slope and ke0 are fitted on a log scale (always positive); E0 and
Emax are fitted as they are, so Emax can be negative for an inhibitory effect.
EC50, the Hill slope and ke0 get default limits so that poorly informative data
cannot send them to absurd values:

    EC50 between 1e-3 x the lowest non-zero concentration and 1e3 x the highest
    Hill slope between 0.1 and 10
    ke0 between 0.001 and 100 per hour (equilibration half-life 0.007 h to 693 h)

When ke0 is fitted without a starting value, a coarse scan over ke0 picks the
starting point, because the error surface of a delayed effect has local minima.

`se` holds standard errors in each parameter's own unit. `cv_percent` expresses them
relative to the estimate, which is only meaningful for EC50, the Hill slope and ke0:
a baseline or Emax near zero gives an enormous, meaningless CV.

A parameter that finishes on a limit is reported in `at_bound`: treat it as
"the data cannot pin this down", not as a measurement. A ke0 on its upper limit
means the data show no detectable delay.

If the measured concentrations never get near EC50, Emax and EC50 trade off
against each other and their uncertainties are large. The fit will still return
numbers, so read the CV% and the notes before trusting them.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Callable

import numpy as np
from scipy.optimize import least_squares

from ligentry.pkpd.model import EffectParameters

EFFECT_FITTABLE = ("e0", "emax", "ec50", "hill", "ke0")
_LOG_PARAMS = ("ec50", "hill", "ke0")
DEFAULT_HILL_BOUNDS = (0.1, 10.0)
DEFAULT_KE0_BOUNDS = (1e-3, 100.0)
_KE0_SCAN = np.geomspace(0.01, 30.0, 12)


@dataclass
class EffectFitResult:
    effect: EffectParameters  # fitted values for all parameters (ke0 is None for a direct effect)
    start: EffectParameters  # values the optimiser started from
    free: tuple[str, ...]
    estimates: dict[str, float]
    se: dict[str, float] | None  # standard errors (None without spare degrees of freedom)
    cv_percent: dict[str, float | None] | None  # None entry when the estimate is exactly 0
    at_bound: tuple[str, ...]  # free parameters that finished on a limit
    rmse: float  # root-mean-square error, in the unit of the measured effect
    effect_range: float  # largest minus smallest measured effect, to judge uncertainties against
    n_obs: int
    success: bool
    message: str


def guess_effect_start(conc: np.ndarray, effect: np.ndarray) -> EffectParameters:
    """Starting values from the data: baseline at the lowest concentrations, Emax from
    the highest, EC50 at the median non-zero concentration, Hill slope 1."""
    conc = np.asarray(conc, dtype=float)
    effect = np.asarray(effect, dtype=float)
    if conc.size == 0 or conc.shape != effect.shape:
        raise ValueError("conc and effect must be non-empty arrays of the same length")
    order = np.argsort(conc, kind="stable")
    k = max(1, conc.size // 4)
    low = float(np.mean(effect[order[:k]]))
    high = float(np.mean(effect[order[-k:]]))
    change = high - low
    if change == 0.0:
        change = float(np.max(effect) - np.min(effect)) or 1.0
    positive = conc[conc > 0]
    ec50 = float(np.median(positive)) if positive.size else 1.0
    return EffectParameters(e0=low, emax=1.5 * change, ec50=ec50, hill=1.0)


def _predict(values: dict[str, float], conc: np.ndarray) -> np.ndarray:
    ch = conc ** values["hill"]
    return values["e0"] + values["emax"] * ch / (values["ec50"] ** values["hill"] + ch)


def _scan_ke0(
    effect: np.ndarray,
    base: EffectParameters,
    others: tuple[str, ...],
    bounds: dict[str, tuple[float, float]] | None,
    site: Callable[[float], np.ndarray],
    ke0_limits: tuple[float, float],
) -> tuple[float, EffectParameters | None]:
    """Try a grid of ke0 values and return the best one with the other parameters fitted for it."""
    best_k, best_fit, best_rmse = None, None, np.inf
    for k in _KE0_SCAN:
        if not ke0_limits[0] <= k <= ke0_limits[1]:
            continue
        ce = site(float(k))
        if not np.all(np.isfinite(ce)) or float(ce.max()) <= 0.0:
            continue
        try:
            if others:
                fit = fit_effect(ce, effect, start=replace(base, ke0=None), free=others, bounds=bounds)
                rmse, params = fit.rmse, fit.effect
            else:
                params = replace(base, ke0=None)
                values = {n: float(getattr(params, n)) for n in ("e0", "emax", "ec50", "hill")}
                rmse = float(np.sqrt(np.mean((_predict(values, ce) - effect) ** 2)))
        except ValueError:
            continue
        if rmse < best_rmse:
            best_k, best_fit, best_rmse = float(k), params, rmse
    if best_k is None:
        return float(np.sqrt(ke0_limits[0] * ke0_limits[1])), None
    return best_k, best_fit


def fit_effect(
    conc: np.ndarray,
    effect: np.ndarray,
    start: EffectParameters | None = None,
    free: tuple[str, ...] | list[str] = ("e0", "emax", "ec50"),
    bounds: dict[str, tuple[float, float]] | None = None,
    max_nfev: int = 400,
    site: Callable[[float], np.ndarray] | None = None,
) -> EffectFitResult:
    """Fit the chosen effect parameters.

    `conc` is the plasma concentration at each measurement. For a delayed effect (ke0
    free, or `start.ke0` given) `site(ke0)` must return the effect-site concentration
    at the same times.
    """
    conc = np.asarray(conc, dtype=float)
    effect = np.asarray(effect, dtype=float)
    free = tuple(free)
    if conc.ndim != 1 or conc.shape != effect.shape:
        raise ValueError("conc and effect must be 1-D arrays of the same length")
    if not (np.all(np.isfinite(conc)) and np.all(np.isfinite(effect))):
        raise ValueError("concentrations and effects must be finite numbers")
    if np.any(conc < 0):
        raise ValueError("concentrations must not be negative")
    if not free:
        raise ValueError("at least one parameter must be free")
    if len(set(free)) != len(free):
        raise ValueError("a parameter is listed twice in 'free'")
    for name in free:
        if name not in EFFECT_FITTABLE:
            raise ValueError(f"cannot fit '{name}'; choose from {EFFECT_FITTABLE}")
    if conc.size < len(free):
        raise ValueError(f"{conc.size} observations cannot determine {len(free)} parameters")
    c_max = float(conc.max())
    if c_max <= 0.0:
        raise ValueError(
            "All concentrations at the measurement times are zero, so there is no "
            "concentration-effect relationship to fit. Check the dosing and the times."
        )
    c_lo = float(conc[conc > 0].min())

    if start is None:
        start = guess_effect_start(conc, effect)

    delayed = "ke0" in free or start.ke0 is not None
    if delayed and site is None:
        raise ValueError(
            "A delayed effect needs the PK model to give the effect-site concentration "
            "(the 'site' function)."
        )

    limits: dict[str, tuple[float, float]] = {
        "e0": (-np.inf, np.inf),
        "emax": (-np.inf, np.inf),
        "ec50": (1e-3 * c_lo, 1e3 * c_max),
        "hill": DEFAULT_HILL_BOUNDS,
        "ke0": DEFAULT_KE0_BOUNDS,
    }
    for name, (lo, hi) in (bounds or {}).items():
        if name in free:
            if not lo < hi:
                raise ValueError(f"bounds for {name} must satisfy low < high")
            if name in _LOG_PARAMS and lo <= 0:
                raise ValueError(f"bounds for {name} must be positive")
            limits[name] = (float(lo), float(hi))

    if "ke0" in free and start.ke0 is None:
        # no starting value: let a coarse scan choose it (the error surface has local minima)
        others = tuple(n for n in free if n != "ke0")
        k0, best = _scan_ke0(effect, start, others, bounds, site, limits["ke0"])
        start = replace(best if best is not None else start, ke0=k0)

    def to_x(name: str, value: float) -> float:
        return float(np.log(value)) if name in _LOG_PARAMS else float(value)

    lower = np.array([to_x(n, limits[n][0]) if np.isfinite(limits[n][0]) else -np.inf for n in free])
    upper = np.array([to_x(n, limits[n][1]) if np.isfinite(limits[n][1]) else np.inf for n in free])
    x0 = np.array([to_x(n, getattr(start, n)) for n in free])
    pad = 1e-9 * np.maximum(1.0, np.abs(x0))
    x0 = np.clip(x0, lower + pad, upper - pad)

    fixed = {
        n: float(getattr(start, n))
        for n in EFFECT_FITTABLE
        if n not in free and getattr(start, n) is not None
    }
    fixed_drive = site(fixed["ke0"]) if "ke0" in fixed else None  # ke0 held at a known value

    def native(x: np.ndarray) -> dict[str, float]:
        values = dict(fixed)
        for n, v in zip(free, x):
            values[n] = float(np.exp(v)) if n in _LOG_PARAMS else float(v)
        return values

    def drive(values: dict[str, float]) -> np.ndarray:
        if "ke0" not in values:
            return conc
        return site(values["ke0"]) if "ke0" in free else fixed_drive

    def predict(values: dict[str, float]) -> np.ndarray:
        return _predict(values, drive(values))

    def residuals(x: np.ndarray) -> np.ndarray:
        return predict(native(x)) - effect

    kwargs = {"diff_step": 1e-6} if delayed else {}  # ODE noise swamps the default tiny step
    sol = least_squares(
        residuals, x0, bounds=(lower, upper), max_nfev=max_nfev, x_scale="jac", **kwargs
    )
    values = native(sol.x)
    fitted = EffectParameters(
        e0=values["e0"], emax=values["emax"], ec50=values["ec50"], hill=values["hill"],
        ke0=values.get("ke0"),
    )
    rmse = float(np.sqrt(np.mean((predict(values) - effect) ** 2)))

    at_bound = tuple(
        n
        for i, n in enumerate(free)
        if (np.isfinite(lower[i]) and abs(sol.x[i] - lower[i]) <= 1e-6 * max(1.0, abs(lower[i])))
        or (np.isfinite(upper[i]) and abs(sol.x[i] - upper[i]) <= 1e-6 * max(1.0, abs(upper[i])))
    )

    se = None
    cv: dict[str, float | None] | None = None
    dof = conc.size - len(free)
    if dof > 0:
        s2 = float(np.sum(sol.fun**2) / dof)
        cov_x = s2 * np.linalg.pinv(sol.jac.T @ sol.jac)
        # delta method: a log-scale parameter has d(native)/dx = native value
        scale = np.array([values[n] if n in _LOG_PARAMS else 1.0 for n in free])
        cov = cov_x * np.outer(scale, scale)
        se = {n: float(np.sqrt(max(cov[i, i], 0.0))) for i, n in enumerate(free)}
        cv = {
            n: (None if values[n] == 0.0 else float(100.0 * se[n] / abs(values[n])))
            for n in free
        }

    return EffectFitResult(
        effect=fitted,
        start=replace(start),
        free=free,
        estimates={n: values[n] for n in free},
        se=se,
        cv_percent=cv,
        at_bound=at_bound,
        rmse=rmse,
        effect_range=float(effect.max() - effect.min()),
        n_obs=int(conc.size),
        success=bool(sol.success),
        message=str(sol.message),
    )
