"""Fit pharmacokinetic parameters to measured plasma concentrations.

Parameters (except non-negative oral lag time tlag) are fitted on a log scale,
so they always stay positive. With oral data, ka and the elimination rate can
trade places and still fit the curve equally well (the "flip-flop" problem), so
give a starting guess that is physically sensible, and use `bounds` when needed.

ka and tlag do not apply to IV bolus or IV infusion dosing and cannot be fitted there.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import least_squares

from ligentry.pkpd.model import PKParameters, Regimen, simulate

FITTABLE = (
    "ka",
    "volume",
    "clearance",
    "vmax",
    "km",
    "peripheral_volume",
    "intercompartmental_clearance",
    "tlag",
)
_LINEAR_PK_PARAMS = ("tlag",)
WEIGHTINGS = ("none", "proportional")

_WEIGHT_FLOOR = 0.05


@dataclass
class FitResult:
    params: PKParameters
    free: tuple[str, ...]
    estimates: dict[str, float]
    cv_percent: dict[str, float | None] | None  # None when there are no spare degrees of freedom
    rmse: float  # root-mean-square error in concentration units (mg/L)
    n_obs: int
    success: bool
    message: str
    route: str = "oral"  # dosing route the fit was made with
    regimen: Regimen | None = None  # updated regimen (carries fitted tlag if fitted)


def fit_concentration(
    times: np.ndarray,
    conc: np.ndarray,
    regimen: Regimen,
    start_params: PKParameters,
    free: tuple[str, ...] | list[str] = ("ka", "volume", "clearance"),
    weighting: str = "proportional",
    bounds: dict[str, tuple[float, float]] | None = None,
    max_nfev: int = 400,
) -> FitResult:
    times = np.asarray(times, dtype=float)
    conc = np.asarray(conc, dtype=float)
    free = tuple(free)
    if times.shape != conc.shape or times.ndim != 1:
        raise ValueError("times and conc must be 1-D arrays of the same length")
    if np.any(conc < 0):
        raise ValueError("concentrations must not be negative")
    if weighting not in WEIGHTINGS:
        raise ValueError(f"weighting must be one of {WEIGHTINGS}")
    if not free:
        raise ValueError("at least one parameter must be free")
    for name in free:
        if name not in FITTABLE:
            raise ValueError(f"cannot fit '{name}'; choose from {FITTABLE}")
        if name != "tlag" and getattr(start_params, name) is None:
            raise ValueError(
                f"'{name}' has no starting value (it does not apply to "
                f"{start_params.elimination} elimination)"
            )
    if regimen.route in ("iv_bolus", "iv_infusion"):
        for invalid in ("ka", "tlag"):
            if invalid in free:
                raise ValueError(
                    f"{invalid} does not apply to {regimen.route} dosing, so it cannot be fitted"
                )
    if times.size < len(free):
        raise ValueError(f"{times.size} observations cannot determine {len(free)} parameters")

    max_tlag = float(regimen.interval) if regimen.n_doses > 1 else max(float(times.max()), 1.0)
    lower = np.full(len(free), -np.inf)
    upper = np.full(len(free), np.inf)
    for i, name in enumerate(free):
        if name == "tlag":
            lower[i], upper[i] = 0.0, max_tlag

    for name, (lo, hi) in (bounds or {}).items():
        if name in free:
            i = free.index(name)
            if name in _LINEAR_PK_PARAMS:
                lower[i], upper[i] = float(lo), float(hi)
            else:
                lower[i], upper[i] = np.log(lo), np.log(hi)

    x0_list: list[float] = []
    for name in free:
        if name == "tlag":
            val = float(regimen.tlag)
            if val <= 0.0:
                positive_t = times[times > 0]
                val = min(0.5 * float(positive_t[0]) if positive_t.size else 0.1, 0.5 * max_tlag)
            x0_list.append(val)
        else:
            x0_list.append(float(np.log(float(getattr(start_params, name)))))
    x0 = np.array(x0_list, dtype=float)
    x0 = np.clip(x0, lower + 1e-9, upper - 1e-9)

    if weighting == "proportional":
        scale = np.maximum(conc, _WEIGHT_FLOOR * max(conc.max(), 1e-12))
    else:
        scale = np.ones_like(conc)

    def make(x: np.ndarray) -> tuple[PKParameters, Regimen]:
        pk_updates = {}
        reg_updates = {}
        for n, v in zip(free, x):
            if n == "tlag":
                reg_updates["tlag"] = max(0.0, float(v))
            else:
                pk_updates[n] = float(np.exp(v))
        pk_obj = replace(start_params, **pk_updates) if pk_updates else start_params
        reg_obj = replace(regimen, **reg_updates) if reg_updates else regimen
        return pk_obj, reg_obj

    def residuals(x: np.ndarray) -> np.ndarray:
        pk_obj, reg_obj = make(x)
        return (simulate(pk_obj, reg_obj, times).conc - conc) / scale

    kwargs = {"diff_step": 1e-6} if "tlag" in free else {}
    sol = least_squares(
        residuals, x0, bounds=(lower, upper), max_nfev=max_nfev, x_scale=1.0, **kwargs
    )
    fitted, fitted_reg = make(sol.x)
    pred = simulate(fitted, fitted_reg, times).conc
    rmse = float(np.sqrt(np.mean((pred - conc) ** 2)))

    estimates: dict[str, float] = {}
    for n in free:
        estimates[n] = float(fitted_reg.tlag) if n == "tlag" else float(getattr(fitted, n))

    cv: dict[str, float | None] | None = None
    dof = times.size - len(free)
    if dof > 0:
        s2 = float(np.sum(sol.fun**2) / dof)
        cov = s2 * np.linalg.pinv(sol.jac.T @ sol.jac)
        cv = {}
        for i, n in enumerate(free):
            sd = float(np.sqrt(max(cov[i, i], 0.0)))
            if n in _LINEAR_PK_PARAMS:
                est = estimates[n]
                cv[n] = None if est == 0.0 else float(100.0 * sd / abs(est))
            else:
                cv[n] = float(100.0 * sd)

    return FitResult(
        params=fitted,
        free=free,
        estimates=estimates,
        cv_percent=cv,
        rmse=rmse,
        n_obs=int(times.size),
        success=bool(sol.success),
        message=str(sol.message),
        route=regimen.route,
        regimen=fitted_reg,
    )
