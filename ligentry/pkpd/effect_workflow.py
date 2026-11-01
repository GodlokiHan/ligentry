"""Glue between the PK/PD form values and the effect fit (no Qt, easy to test).

The effect model is driven by plasma concentration, so the concentration at each
effect measurement time comes from the PK model currently entered in the form.
Fit the PK parameters first (from concentration data) or enter known ones, then
fit the effect parameters here.

The model is a DIRECT effect: the effect follows the plasma concentration with no
delay. If the real effect lags behind the concentration (hysteresis), the fit
will be poor and the residuals will show a systematic pattern.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np

from ligentry.pkpd.effect_data import EffectObservations
from ligentry.pkpd.effect_fit import (
    EFFECT_FITTABLE,
    EffectFitResult,
    fit_effect,
    guess_effect_start,
)
from ligentry.pkpd.model import EffectParameters, effect_from_conc, simulate
from ligentry.pkpd.workflow import ModelInputs, build_model

# parameter name -> (short label, unit); the effect unit is whatever the user measured
EFFECT_PARAMETER_INFO = {
    "e0": ("E0", "effect"),
    "emax": ("Emax", "effect"),
    "ec50": ("EC50", "mg/L"),
    "hill": ("Hill", ""),
}

_POOR_CV_PERCENT = 50.0  # for EC50 and the Hill slope
_POOR_SE_FRACTION = 0.25  # for E0 and Emax, as a fraction of the measured effect range


def effect_free_choices() -> tuple[str, ...]:
    return EFFECT_FITTABLE


def concentrations_at(inputs: ModelInputs, obs: EffectObservations) -> np.ndarray:
    """Plasma concentration (mg/L) at each effect measurement time, from the form's PK model."""
    pk, regimen, _effect = build_model(inputs)
    return simulate(pk, regimen, obs.times).conc


def form_effect(inputs: ModelInputs) -> EffectParameters:
    return EffectParameters(e0=inputs.e0, emax=inputs.emax, ec50=inputs.ec50, hill=inputs.hill)


def run_effect_fit(
    inputs: ModelInputs,
    obs: EffectObservations,
    free: tuple[str, ...] | list[str] = ("e0", "emax", "ec50"),
    use_form_start: bool = False,
) -> tuple[EffectFitResult, np.ndarray]:
    """Fit the chosen effect parameters. Returns (result, concentrations at the data times).

    By default the free parameters start from values guessed from the data and the
    fixed ones keep the values in the form. With use_form_start=True every parameter
    starts from the form values.
    """
    free = tuple(free)
    conc = concentrations_at(inputs, obs)
    current = form_effect(inputs)
    if use_form_start:
        start = current
    else:
        guess = guess_effect_start(conc, obs.effect)
        start = replace(guess, **{n: getattr(current, n) for n in EFFECT_FITTABLE if n not in free})
    return fit_effect(conc, obs.effect, start=start, free=free), conc


def apply_effect_fit(inputs: ModelInputs, result: EffectFitResult) -> ModelInputs:
    """Copy the fitted effect parameters back into the form values and switch the effect on."""
    e = result.effect
    return replace(inputs, use_effect=True, e0=e.e0, emax=e.emax, ec50=e.ec50, hill=e.hill)


def predict_effects(inputs: ModelInputs, obs: EffectObservations) -> tuple[np.ndarray, np.ndarray]:
    """(concentration, predicted effect) at the data times for the form's current model."""
    conc = concentrations_at(inputs, obs)
    return conc, effect_from_conc(conc, form_effect(inputs))


def effect_fit_table_rows(result: EffectFitResult) -> list[tuple[str, str, str, str]]:
    """(parameter, estimate, unit, uncertainty) as display strings, in fit order.

    The uncertainty is the standard error in the parameter's own unit, not a CV %:
    a baseline or Emax close to zero would give a meaningless huge percentage.
    """
    rows = []
    for name in result.free:
        label, unit = EFFECT_PARAMETER_INFO[name]
        se = None if result.se is None else result.se.get(name)
        rows.append(
            (label, f"{result.estimates[name]:.4g}", unit, "n/a" if se is None else f"\u00b1{se:.3g}")
        )
    return rows


def effect_fit_notes(result: EffectFitResult, conc: np.ndarray) -> list[str]:
    """Plain-language warnings to show next to an effect fit."""
    notes: list[str] = []
    if not result.success:
        notes.append(f"The optimiser did not converge: {result.message}")
    if result.cv_percent is None:
        notes.append(
            "Too few observations to estimate parameter uncertainty "
            "(more points than fitted parameters are needed)."
        )
    else:
        poor = []
        for name in result.free:
            label = EFFECT_PARAMETER_INFO[name][0]
            if name in ("ec50", "hill"):
                cv = result.cv_percent.get(name)
                if cv is not None and cv > _POOR_CV_PERCENT:
                    poor.append(f"{label} (CV {cv:.0f}%)")
            elif result.effect_range > 0 and result.se[name] > _POOR_SE_FRACTION * result.effect_range:
                poor.append(f"{label} (uncertainty above a quarter of the measured effect range)")
        if poor:
            notes.append("Poorly determined: " + "; ".join(poor) + ".")
    if result.at_bound:
        names = ", ".join(EFFECT_PARAMETER_INFO[n][0] for n in result.at_bound)
        notes.append(
            f"{names} finished on its allowed limit, so the data cannot determine it. "
            "Fix it at a known value instead of trusting the fitted one."
        )
    conc = np.asarray(conc, dtype=float)
    if "ec50" in result.free and conc.size and float(conc.max()) > 0:
        c_max = float(conc.max())
        positive = conc[conc > 0]
        ec50 = result.effect.ec50
        if ec50 > c_max:
            notes.append(
                f"EC50 ({ec50:.4g} mg/L) is above the highest concentration at the data times "
                f"({c_max:.4g} mg/L): the effect never approaches saturation, so Emax and EC50 "
                "are only weakly determined."
            )
        elif positive.size and ec50 < float(positive.min()):
            notes.append(
                f"EC50 ({ec50:.4g} mg/L) is below the lowest non-zero concentration "
                f"({float(positive.min()):.4g} mg/L): the effect is saturated everywhere "
                "that was measured."
            )
    notes.append(
        "The effect is modelled as following the plasma concentration with no delay. A delayed "
        "effect shows as a systematic pattern in the residuals."
    )
    return notes
