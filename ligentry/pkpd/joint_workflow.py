"""Glue between ModelInputs and simultaneous (joint) PK/PD fitting (no Qt)."""
from __future__ import annotations

from dataclasses import replace

import numpy as np

from ligentry.pkpd.effect_data import EffectObservations
from ligentry.pkpd.effect_workflow import EFFECT_PARAMETER_INFO, form_effect
from ligentry.pkpd.fit import FitResult
from ligentry.pkpd.joint_fit import JointFitResult, fit_joint_pkpd
from ligentry.pkpd.model import effect_site_concentration, equilibration_half_life, simulate
from ligentry.pkpd.observations import Observations
from ligentry.pkpd.workflow import ModelInputs, PARAMETER_INFO, build_model, fit_notes

_LOG_SCALE_PD = ("ec50", "hill", "ke0")
_POOR_CV_PERCENT = 50.0
_POOR_SE_FRACTION = 0.25


def run_joint_fit(
    inputs: ModelInputs,
    pk_obs: Observations,
    pd_obs: EffectObservations,
    free_pk: tuple[str, ...] | list[str],
    free_pd: tuple[str, ...] | list[str] = ("e0", "emax", "ec50"),
    weighting: str = "proportional",
    use_form_start: bool = False,
) -> tuple[JointFitResult, np.ndarray]:
    """Simultaneously fit free PK and PD parameters to concentration and effect datasets.

    Returns (result, driving_concentrations_at_pd_times).
    """
    free_pk = tuple(free_pk)
    free_pd = tuple(free_pd)
    pk, regimen, _eff = build_model(inputs)
    delayed = "ke0" in free_pd or inputs.delayed_effect

    current_eff = form_effect(inputs)
    if delayed and current_eff.ke0 is None:
        current_eff = replace(current_eff, ke0=inputs.ke0)

    result = fit_joint_pkpd(
        pk_obs.times,
        pk_obs.conc,
        pd_obs.times,
        pd_obs.effect,
        regimen=regimen,
        start_pk=pk,
        start_effect=current_eff,
        free_pk=free_pk,
        free_pd=free_pd,
        weighting=weighting,
        use_form_start=use_form_start,
    )
    if result.effect.ke0 is not None:
        drive = effect_site_concentration(result.pk, regimen, pd_obs.times, result.effect.ke0)
    else:
        drive = simulate(result.pk, regimen, pd_obs.times).conc
    return result, drive


def apply_joint_fit(inputs: ModelInputs, result: JointFitResult) -> ModelInputs:
    """Copy both fitted PK and PD estimates back into ModelInputs and enable the effect."""
    e = result.effect
    changes: dict = dict(
        **result.pk_estimates,
        use_effect=True,
        e0=e.e0,
        emax=e.emax,
        ec50=e.ec50,
        hill=e.hill,
    )
    if e.ke0 is not None:
        changes.update(delayed_effect=True, ke0=e.ke0)
    return replace(inputs, **changes)


def joint_fit_table_rows(result: JointFitResult) -> list[tuple[str, str, str, str]]:
    """(parameter, estimate, unit, uncertainty) rows for the unified PK + PD results table."""
    rows: list[tuple[str, str, str, str]] = []
    for name in result.free_pk:
        label, unit = PARAMETER_INFO[name]
        cv = None if result.cv_percent is None else result.cv_percent.get(name)
        unc = "n/a" if cv is None else f"CV {cv:.1f}%"
        rows.append((label, f"{result.pk_estimates[name]:.4g}", unit, unc))

    for name in result.free_pd:
        label, unit = EFFECT_PARAMETER_INFO[name]
        se = None if result.se is None else result.se.get(name)
        unc = "n/a" if se is None else f"\u00b1{se:.3g}"
        rows.append((label, f"{result.pd_estimates[name]:.4g}", unit, unc))
    return rows


def joint_fit_notes(result: JointFitResult, drive_conc: np.ndarray) -> list[str]:
    """Combined PK + PD diagnostic notes for a joint fit."""
    notes: list[str] = []
    if not result.success:
        notes.append(f"The joint optimiser did not converge: {result.message}")

    # Reuse PK structural & identifiability checks via a lightweight FitResult adapter
    pk_cv = (
        {n: result.cv_percent[n] for n in result.free_pk if result.cv_percent.get(n) is not None}
        if result.cv_percent is not None
        else None
    )
    pk_proxy = FitResult(
        params=result.pk,
        free=result.free_pk,
        estimates=result.pk_estimates,
        cv_percent=pk_cv,
        rmse=result.pk_rmse,
        n_obs=result.n_pk_obs + result.n_pd_obs,
        success=result.success,
        message=result.message,
        route=result.route,
    )
    for note in fit_notes(pk_proxy):
        if "Poorly determined" not in note and "Too few observations" not in note:
            notes.append(note)

    if result.cv_percent is None or result.se is None:
        notes.append(
            "Too few combined observations to estimate joint parameter uncertainty "
            "(more total points than fitted PK + PD parameters are needed)."
        )
    else:
        poor: list[str] = []
        for name in result.free_pk:
            cv = result.cv_percent.get(name)
            if cv is not None and cv > _POOR_CV_PERCENT:
                poor.append(f"{PARAMETER_INFO[name][0]} (CV {cv:.0f}%)")
        for name in result.free_pd:
            label = EFFECT_PARAMETER_INFO[name][0]
            if name in _LOG_SCALE_PD:
                cv = result.cv_percent.get(name)
                if cv is not None and cv > _POOR_CV_PERCENT:
                    poor.append(f"{label} (CV {cv:.0f}%)")
            elif (
                result.effect_range > 0
                and result.se[name] > _POOR_SE_FRACTION * result.effect_range
            ):
                poor.append(f"{label} (uncertainty > 25% of effect range)")
        if poor:
            notes.append("Poorly determined: " + "; ".join(poor) + ".")

    if "ke0" in result.at_bound_pd:
        ke0 = result.effect.ke0
        if ke0 is not None and ke0 >= 1.0:
            notes.append(
                "ke0 finished on its upper limit: the data show no detectable delay between "
                "concentration and effect."
            )
        else:
            notes.append(
                "ke0 finished on its lower limit: the delay is longer than these data can resolve."
            )

    other_bound = tuple(n for n in result.at_bound_pd if n != "ke0")
    if other_bound:
        names = ", ".join(EFFECT_PARAMETER_INFO[n][0] for n in other_bound)
        notes.append(f"{names} finished on its allowed limit, so the data cannot determine it.")

    drive_conc = np.asarray(drive_conc, dtype=float)
    if "ec50" in result.free_pd and drive_conc.size and float(drive_conc.max()) > 0:
        c_max = float(drive_conc.max())
        positive = drive_conc[drive_conc > 0]
        ec50 = result.effect.ec50
        what = "effect-site concentration" if result.effect.ke0 is not None else "concentration"
        if ec50 > c_max:
            notes.append(
                f"EC50 ({ec50:.4g} mg/L) is above the highest {what} ({c_max:.4g} mg/L): "
                "Emax and EC50 are only weakly determined."
            )
        elif positive.size and ec50 < float(positive.min()):
            notes.append(
                f"EC50 ({ec50:.4g} mg/L) is below the lowest non-zero {what} "
                f"({float(positive.min()):.4g} mg/L): the effect is saturated."
            )

    if result.effect.ke0 is not None:
        notes.append(
            f"Joint effect-compartment fit: equilibration half-life "
            f"{equilibration_half_life(result.effect.ke0):.3g} h (ke0 {result.effect.ke0:.3g} 1/h), "
            "with PK uncertainty propagated into PD estimates."
        )
    return notes
