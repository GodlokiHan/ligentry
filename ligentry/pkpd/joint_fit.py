"""Joint (simultaneous) pharmacokinetic and pharmacodynamic parameter fitting."""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import least_squares

from ligentry.pkpd.effect_fit import (
    DEFAULT_HILL_BOUNDS,
    DEFAULT_KE0_BOUNDS,
    EFFECT_FITTABLE,
    fit_effect,
    guess_effect_start,
)
from ligentry.pkpd.fit import (
    FITTABLE as PK_FITTABLE,
    WEIGHTINGS,
    _WEIGHT_FLOOR,
    fit_concentration,
)
from ligentry.pkpd.model import (
    EffectParameters,
    PKParameters,
    Regimen,
    effect_from_conc,
    effect_site_concentration,
    simulate,
)

_PK_LOG_PARAMS = tuple(n for n in PK_FITTABLE if n != "tlag")
_PD_LOG_PARAMS = ("ec50", "hill", "ke0")


@dataclass
class JointFitResult:
    pk: PKParameters
    effect: EffectParameters
    free_pk: tuple[str, ...]
    free_pd: tuple[str, ...]
    pk_estimates: dict[str, float]
    pd_estimates: dict[str, float]
    se: dict[str, float] | None
    cv_percent: dict[str, float | None] | None
    at_bound_pd: tuple[str, ...]
    pk_rmse: float  # mg/L
    pd_rmse: float  # effect units
    effect_range: float
    n_pk_obs: int
    n_pd_obs: int
    success: bool
    message: str
    route: str = "oral"
    regimen: Regimen | None = None


def fit_joint_pkpd(
    pk_times: np.ndarray,
    pk_conc: np.ndarray,
    pd_times: np.ndarray,
    pd_effect: np.ndarray,
    regimen: Regimen,
    start_pk: PKParameters,
    start_effect: EffectParameters | None = None,
    free_pk: tuple[str, ...] | list[str] = ("ka", "volume", "clearance"),
    free_pd: tuple[str, ...] | list[str] = ("e0", "emax", "ec50"),
    weighting: str = "proportional",
    use_form_start: bool = False,
    pk_bounds: dict[str, tuple[float, float]] | None = None,
    pd_bounds: dict[str, tuple[float, float]] | None = None,
    max_nfev: int = 600,
) -> JointFitResult:
    """Simultaneously fit PK and PD parameters to concentration and effect data."""
    pk_times = np.asarray(pk_times, dtype=float)
    pk_conc = np.asarray(pk_conc, dtype=float)
    pd_times = np.asarray(pd_times, dtype=float)
    pd_effect = np.asarray(pd_effect, dtype=float)
    free_pk = tuple(free_pk)
    free_pd = tuple(free_pd)

    if pk_times.shape != pk_conc.shape or pk_times.ndim != 1 or pk_times.size == 0:
        raise ValueError("pk_times and pk_conc must be non-empty 1-D arrays of the same length")
    if pd_times.shape != pd_effect.shape or pd_times.ndim != 1 or pd_times.size == 0:
        raise ValueError("pd_times and pd_effect must be non-empty 1-D arrays of the same length")
    if np.any(pk_conc < 0):
        raise ValueError("concentrations must not be negative")
    if weighting not in WEIGHTINGS:
        raise ValueError(f"weighting must be one of {WEIGHTINGS}")
    if not free_pk or not free_pd:
        raise ValueError("joint fitting requires at least one free PK and one free PD parameter")
    if len(set(free_pk)) != len(free_pk) or len(set(free_pd)) != len(free_pd):
        raise ValueError("duplicate parameter names in free_pk or free_pd")

    for name in free_pk:
        if name not in PK_FITTABLE:
            raise ValueError(f"cannot fit PK parameter '{name}'")
        if name != "tlag" and getattr(start_pk, name) is None:
            raise ValueError(
                f"'{name}' has no starting value (not applicable to {start_pk.elimination} elimination)"
            )
    if regimen.route in ("iv_bolus", "iv_infusion"):
        for invalid in ("ka", "tlag"):
            if invalid in free_pk:
                raise ValueError(
                    f"{invalid} does not apply to {regimen.route} dosing, so it cannot be fitted"
                )

    for name in free_pd:
        if name not in EFFECT_FITTABLE:
            raise ValueError(f"cannot fit PD parameter '{name}'")

    total_obs = pk_times.size + pd_times.size
    total_free = len(free_pk) + len(free_pd)
    if total_obs < total_free:
        raise ValueError(f"{total_obs} total observations cannot determine {total_free} parameters")

    # Stage 1: warm-start & residual scale estimation via fast sequential pass
    init_pk = start_pk
    init_reg = regimen
    if not use_form_start:
        try:
            seq_pk = fit_concentration(
                pk_times, pk_conc, regimen, start_pk,
                free=free_pk, weighting=weighting, bounds=pk_bounds, max_nfev=250,
            )
            init_pk = seq_pk.params
            if seq_pk.regimen is not None:
                init_reg = seq_pk.regimen
        except Exception:
            init_pk = start_pk
            init_reg = regimen

    plasma_pd_init = simulate(init_pk, init_reg, pd_times).conc
    c_max = float(plasma_pd_init.max())
    if c_max <= 0.0:
        raise ValueError(
            "All predicted concentrations at the effect times are zero. Check the dosing and times."
        )
    c_lo = float(plasma_pd_init[plasma_pd_init > 0].min())

    delayed = "ke0" in free_pd or (start_effect is not None and start_effect.ke0 is not None)
    site_fn = (
        (lambda k: effect_site_concentration(init_pk, init_reg, pd_times, k)) if delayed else None
    )

    if start_effect is None:
        base_eff = guess_effect_start(plasma_pd_init, pd_effect)
    else:
        base_eff = start_effect
    if delayed and base_eff.ke0 is None:
        base_eff = replace(base_eff, ke0=0.5 if use_form_start else None)

    init_eff = base_eff
    if not use_form_start:
        try:
            if delayed and "ke0" not in free_pd and init_eff.ke0 is not None:
                guide = site_fn(init_eff.ke0)
            else:
                guide = plasma_pd_init
            guess = guess_effect_start(guide, pd_effect)
            keep = {n: getattr(base_eff, n) for n in EFFECT_FITTABLE if n not in free_pd}
            warm_start = replace(guess, **keep)
            seq_pd = fit_effect(
                plasma_pd_init, pd_effect, start=warm_start,
                free=free_pd, bounds=pd_bounds, max_nfev=250, site=site_fn,
            )
            init_eff = seq_pd.effect
        except Exception:
            if init_eff.ke0 is None and delayed:
                init_eff = replace(init_eff, ke0=0.5)

    if delayed and init_eff.ke0 is None:
        init_eff = replace(init_eff, ke0=0.5)

    # Estimate relative noise scales (sigma_pk, sigma_pd) for balanced joint weighting
    if weighting == "proportional":
        pk_point_scale = np.maximum(pk_conc, _WEIGHT_FLOOR * max(pk_conc.max(), 1e-12))
    else:
        pk_point_scale = np.ones_like(pk_conc)

    pk_pred_init = simulate(init_pk, init_reg, pk_times).conc
    pk_w_res = (pk_pred_init - pk_conc) / pk_point_scale
    dof_pk = max(1, pk_times.size - len(free_pk))
    sigma_pk = float(np.sqrt(np.sum(pk_w_res**2) / dof_pk))
    if not np.isfinite(sigma_pk) or sigma_pk <= 1e-8:
        sigma_pk = 0.05 if weighting == "proportional" else max(1e-3, 0.02 * float(pk_conc.max()))

    if delayed:
        drive_init = effect_site_concentration(init_pk, init_reg, pd_times, init_eff.ke0)
    else:
        drive_init = simulate(init_pk, init_reg, pd_times).conc
    pd_pred_init = effect_from_conc(drive_init, init_eff)
    dof_pd = max(1, pd_times.size - len(free_pd))
    sigma_pd = float(np.sqrt(np.sum((pd_pred_init - pd_effect) ** 2) / dof_pd))
    eff_span = float(pd_effect.max() - pd_effect.min())
    if not np.isfinite(sigma_pd) or sigma_pd <= 1e-8:
        sigma_pd = max(1e-3, 0.02 * (eff_span if eff_span > 0 else 1.0))

    max_tlag = float(regimen.interval) if regimen.n_doses > 1 else max(float(pk_times.max()), 1.0)
    pk_lower = np.full(len(free_pk), -np.inf)
    pk_upper = np.full(len(free_pk), np.inf)
    for i, name in enumerate(free_pk):
        if name == "tlag":
            pk_lower[i], pk_upper[i] = 0.0, max_tlag

    for name, (lo, hi) in (pk_bounds or {}).items():
        if name in free_pk:
            i = free_pk.index(name)
            if name == "tlag":
                pk_lower[i], pk_upper[i] = float(lo), float(hi)
            else:
                pk_lower[i], pk_upper[i] = np.log(lo), np.log(hi)

    pd_limits: dict[str, tuple[float, float]] = {
        "e0": (-np.inf, np.inf),
        "emax": (-np.inf, np.inf),
        "ec50": (1e-3 * c_lo, 1e3 * c_max),
        "hill": DEFAULT_HILL_BOUNDS,
        "ke0": DEFAULT_KE0_BOUNDS,
    }
    for name, (lo, hi) in (pd_bounds or {}).items():
        if name in free_pd:
            if not lo < hi:
                raise ValueError(f"bounds for {name} must satisfy low < high")
            if name in _PD_LOG_PARAMS and lo <= 0:
                raise ValueError(f"bounds for {name} must be positive")
            pd_limits[name] = (float(lo), float(hi))

    def to_pd_x(name: str, val: float) -> float:
        return float(np.log(val)) if name in _PD_LOG_PARAMS else float(val)

    pd_lower = np.array(
        [to_pd_x(n, pd_limits[n][0]) if np.isfinite(pd_limits[n][0]) else -np.inf for n in free_pd]
    )
    pd_upper = np.array(
        [to_pd_x(n, pd_limits[n][1]) if np.isfinite(pd_limits[n][1]) else np.inf for n in free_pd]
    )

    x0_pk_list: list[float] = []
    for n in free_pk:
        if n == "tlag":
            val = float(init_reg.tlag)
            if val <= 0.0:
                pos_t = pk_times[pk_times > 0]
                val = min(0.5 * float(pos_t[0]) if pos_t.size else 0.1, 0.5 * max_tlag)
            x0_pk_list.append(val)
        else:
            x0_pk_list.append(float(np.log(float(getattr(init_pk, n)))))
    x0_pk = np.array(x0_pk_list, dtype=float)
    x0_pd = np.array([to_pd_x(n, float(getattr(init_eff, n))) for n in free_pd])
    lower = np.concatenate([pk_lower, pd_lower])
    upper = np.concatenate([pk_upper, pd_upper])
    x0 = np.concatenate([x0_pk, x0_pd])
    pad = 1e-9 * np.maximum(1.0, np.abs(x0))
    x0 = np.clip(x0, lower + pad, upper - pad)

    n_pk = len(free_pk)
    fixed_pd = {
        n: float(getattr(init_eff, n))
        for n in EFFECT_FITTABLE
        if n not in free_pd and getattr(init_eff, n) is not None
    }

    def unpack(x: np.ndarray) -> tuple[PKParameters, Regimen, EffectParameters]:
        pk_vals = {}
        reg_vals = {}
        for n, v in zip(free_pk, x[:n_pk]):
            if n == "tlag":
                reg_vals["tlag"] = max(0.0, float(v))
            else:
                pk_vals[n] = float(np.exp(v))
        pk_obj = replace(init_pk, **pk_vals) if pk_vals else init_pk
        reg_obj = replace(init_reg, **reg_vals) if reg_vals else init_reg

        pd_vals = dict(fixed_pd)
        for n, v in zip(free_pd, x[n_pk:]):
            pd_vals[n] = float(np.exp(v)) if n in _PD_LOG_PARAMS else float(v)
        eff_obj = EffectParameters(
            e0=pd_vals["e0"],
            emax=pd_vals["emax"],
            ec50=pd_vals["ec50"],
            hill=pd_vals["hill"],
            ke0=pd_vals.get("ke0"),
        )
        return pk_obj, reg_obj, eff_obj

    all_times, inv_idx = np.unique(np.concatenate([pk_times, pd_times]), return_inverse=True)
    pk_idx = inv_idx[: pk_times.size]
    pd_idx = inv_idx[pk_times.size :]

    def joint_residuals(x: np.ndarray) -> np.ndarray:
        pk_obj, reg_obj, eff_obj = unpack(x)
        sim = simulate(pk_obj, reg_obj, all_times, eff_obj)
        r_pk = (sim.conc[pk_idx] - pk_conc) / (pk_point_scale * sigma_pk)
        r_pd = (sim.effect[pd_idx] - pd_effect) / sigma_pd
        return np.concatenate([r_pk, r_pd])

    sol = least_squares(
        joint_residuals,
        x0,
        bounds=(lower, upper),
        max_nfev=max_nfev,
        ftol=1e-10,
        xtol=1e-10,
        x_scale="jac",
        diff_step=1e-6,
    )
    fitted_pk, fitted_reg, fitted_eff = unpack(sol.x)

    sim_final = simulate(fitted_pk, fitted_reg, all_times, fitted_eff)
    pk_rmse = float(np.sqrt(np.mean((sim_final.conc[pk_idx] - pk_conc) ** 2)))
    pd_rmse = float(np.sqrt(np.mean((sim_final.effect[pd_idx] - pd_effect) ** 2)))

    at_bound_pd = tuple(
        n
        for i, n in enumerate(free_pd)
        if (
            np.isfinite(pd_lower[i])
            and abs(sol.x[n_pk + i] - pd_lower[i]) <= 1e-6 * max(1.0, abs(pd_lower[i]))
        )
        or (
            np.isfinite(pd_upper[i])
            and abs(sol.x[n_pk + i] - pd_upper[i]) <= 1e-6 * max(1.0, abs(pd_upper[i]))
        )
    )

    pk_estimates = {
        n: (float(fitted_reg.tlag) if n == "tlag" else float(getattr(fitted_pk, n)))
        for n in free_pk
    }
    pd_estimates = {n: float(getattr(fitted_eff, n)) for n in free_pd}

    all_free = free_pk + free_pd
    se: dict[str, float] | None = None
    cv: dict[str, float | None] | None = None
    dof = total_obs - total_free
    if dof > 0:
        s2 = float(np.sum(sol.fun**2) / dof)
        cov_x = s2 * np.linalg.pinv(sol.jac.T @ sol.jac)
        native_vals = dict(pk_estimates)
        native_vals.update(pd_estimates)
        log_set = set(_PK_LOG_PARAMS) | set(_PD_LOG_PARAMS)
        scale = np.array([native_vals[n] if n in log_set else 1.0 for n in all_free])
        cov = cov_x * np.outer(scale, scale)
        se = {n: float(np.sqrt(max(cov[i, i], 0.0))) for i, n in enumerate(all_free)}
        cv = {
            n: (None if native_vals[n] == 0.0 else float(100.0 * se[n] / abs(native_vals[n])))
            for n in all_free
        }

    return JointFitResult(
        pk=fitted_pk,
        effect=fitted_eff,
        free_pk=free_pk,
        free_pd=free_pd,
        pk_estimates=pk_estimates,
        pd_estimates=pd_estimates,
        se=se,
        cv_percent=cv,
        at_bound_pd=at_bound_pd,
        pk_rmse=pk_rmse,
        pd_rmse=pd_rmse,
        effect_range=eff_span,
        n_pk_obs=int(pk_times.size),
        n_pd_obs=int(pd_times.size),
        success=bool(sol.success),
        message=str(sol.message),
        route=regimen.route,
        regimen=fitted_reg,
    )
