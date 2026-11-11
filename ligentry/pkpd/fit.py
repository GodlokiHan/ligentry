"""Fit pharmacokinetic parameters to measured plasma concentrations.

Parameters are fitted on a log scale, so they always stay positive. With oral
data, ka and the elimination rate can trade places and still fit the curve
equally well (the "flip-flop" problem), so give a starting guess that is
physically sensible, and use `bounds` (for example a lower bound on ka) when
the absorption is known to be faster than the elimination.

For a two-compartment model the peripheral volume and the intercompartmental
clearance can be fitted as well. They are only well determined when the data
include several samples from the early, rapidly falling distribution phase, and
IV bolus data determine them far better than oral data, where the absorption
phase hides the distribution phase. ka does not apply to IV bolus dosing and
cannot be fitted there.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import least_squares

from ligentry.pkpd.model import PKParameters, Regimen, simulate

FITTABLE = (
    "ka", "volume", "clearance", "vmax", "km", "peripheral_volume", "intercompartmental_clearance",
)
WEIGHTINGS = ("none", "proportional")

# Proportional weights divide by max(observed, this fraction of the largest observed
# concentration), so near-zero samples cannot dominate the fit.
_WEIGHT_FLOOR = 0.05


@dataclass
class FitResult:
    params: PKParameters
    free: tuple[str, ...]
    estimates: dict[str, float]
    cv_percent: dict[str, float] | None  # None when there are no spare degrees of freedom
    rmse: float  # root-mean-square error in concentration units (mg/L)
    n_obs: int
    success: bool
    message: str
    route: str = "oral"  # dosing route the fit was made with


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
        if getattr(start_params, name) is None:
            raise ValueError(
                f"'{name}' has no starting value (it does not apply to "
                f"{start_params.elimination} elimination)"
            )
    if regimen.route == "iv_bolus" and "ka" in free:
        raise ValueError("ka does not apply to IV bolus dosing, so it cannot be fitted")
    if times.size < len(free):
        raise ValueError(f"{times.size} observations cannot determine {len(free)} parameters")

    lower = np.full(len(free), -np.inf)
    upper = np.full(len(free), np.inf)
    for name, (lo, hi) in (bounds or {}).items():
        if name in free:
            i = free.index(name)
            lower[i], upper[i] = np.log(lo), np.log(hi)

    x0 = np.log([float(getattr(start_params, n)) for n in free])
    x0 = np.clip(x0, lower + 1e-9, upper - 1e-9)

    if weighting == "proportional":
        scale = np.maximum(conc, _WEIGHT_FLOOR * max(conc.max(), 1e-12))
    else:
        scale = np.ones_like(conc)

    def make(x: np.ndarray) -> PKParameters:
        return replace(start_params, **{n: float(np.exp(v)) for n, v in zip(free, x)})

    def residuals(x: np.ndarray) -> np.ndarray:
        return (simulate(make(x), regimen, times).conc - conc) / scale

    sol = least_squares(residuals, x0, bounds=(lower, upper), max_nfev=max_nfev, x_scale=1.0)
    fitted = make(sol.x)
    pred = simulate(fitted, regimen, times).conc
    rmse = float(np.sqrt(np.mean((pred - conc) ** 2)))

    cv = None
    dof = times.size - len(free)
    if dof > 0:
        s2 = float(np.sum(sol.fun**2) / dof)
        cov = s2 * np.linalg.pinv(sol.jac.T @ sol.jac)
        cv = {n: float(100.0 * np.sqrt(max(cov[i, i], 0.0))) for i, n in enumerate(free)}

    return FitResult(
        params=fitted,
        free=free,
        estimates={n: float(getattr(fitted, n)) for n in free},
        cv_percent=cv,
        rmse=rmse,
        n_obs=int(times.size),
        success=bool(sol.success),
        message=str(sol.message),
        route=regimen.route,
    )
