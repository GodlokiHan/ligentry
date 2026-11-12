"""Glue between the PK/PD form values and the effect fit (no Qt, easy to test).

The effect model is driven by concentration, so the concentration at each effect
measurement time comes from the PK model currently entered in the form. Fit the PK
parameters first (from concentration data) or enter known ones, then fit the
effect parameters here.

Two effect models are supported:

* DIRECT (default): the effect follows the plasma concentration with no delay. If
  the real effect lags behind the concentration (hysteresis), the fit will be poor
  and the residuals will show a systematic pattern.
* EFFECT COMPARTMENT (the form's "delayed effect" option, parameter ke0): the effect
  follows a hypothetical effect-site concentration Ce that chases plasma,
  dCe/dt = ke0 * (C - Ce). Then ke0 is fitted along with E0, Emax, EC50 and Hill.
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
from ligentry.pkpd.model import (
    EffectParameters,
    effect_from_conc,
    effect_site_concentration,
    equilibration_half_life,
    simulate,
)
from ligentry.pkpd.workflow import ModelInputs, build_model

# parameter name -> (short label, unit); the effect unit is whatever the user measured
EFFECT_PARAMETER_INFO = {
    "e0": ("E0", "effect"),
    "emax": ("Emax", "effect"),
    "ec50": ("EC50", "mg/L"),
    "hill": ("Hill", ""),
    "ke0": ("ke0", "1/h"),
}

_DIRECT_CHOICES = ("e0", "emax", "ec50", "hill")
_POOR_CV_PERCENT = 50.0  # for EC50, the Hill slope and ke0
_POOR_SE_FRACTION = 0.25  # for E0 and Emax, as a fraction of the measured effect range
_LOG_SCALE_PARAMS = ("ec50", "hill", "ke0")


def effect_free_choices(delayed: bool = False) -> tuple[str, ...]:
    """Effect parameters that can be fitted; ke0 is added for a delayed effect."""
    return EFFECT_FITTABLE if delayed else _DIRECT_CHOICES


def concentrations_at(inputs: ModelInputs, obs: EffectObservations) -> np.ndarray:
    """Plasma concentration (mg/L) at each effect measurement time, from the form's PK model."""
    pk, regimen, _effect = build_model(inputs)
    return simulate(pk, regimen, obs.times).conc


def driving_concentrations_at(inputs: ModelInputs, obs: EffectObservations) -> np.ndarray:
    """Concentration that drives the effect at each measurement time (mg/L).

    The plasma concentration for a direct effect, the effect-site concentration Ce
    when the form has a delayed effect.
    """
    if not inputs.delayed_effect:
        return concentrations_at(inputs, obs)
    pk, regimen, _effect = build_model(inputs)
    return effect_site_concentration(pk, regimen, obs.times, inputs.ke0)


def form_effect(inputs: ModelInputs) -> EffectParameters:
    return EffectParameters(
        e0=inputs.e0, emax=inputs.emax, ec50=inputs.ec50, hill=inputs.hill,
        ke0=inputs.ke0 if inputs.delayed_effect else None,
    )


def run_effect_fit(
    inputs: ModelInputs,
    obs: EffectObservations,
    free: tuple[str, ...] | list[str] = ("e0", "emax", "ec50"),
    use_form_start: bool = False,
) -> tuple[EffectFitResult, np.ndarray]:
    """Fit the chosen effect parameters.

    Returns (result, concentrations that drove the fit at the data times): plasma
    concentrations for a direct effect, effect-site concentrations (at the fitted
    ke0) for a delayed one.

    By default the free parameters start from values guessed from the data and the
    fixed ones keep the values in the form. With use_form_start=True every parameter
    starts from the form values. Putting "ke0" in `free` fits an effect compartment
    even when the form's delayed-effect box is off; without a form start ke0 is then
    chosen by a scan, with a form start it begins at the form's ke0.
    """
    free = tuple(free)
    pk, regimen, _effect = build_model(inputs)
    delayed = "ke0" in free or inputs.delayed_effect
    plasma = simulate(pk, regimen, obs.times).conc
    site = (lambda ke0: effect_site_concentration(pk, regimen, obs.times, ke0)) if delayed else None

    current = form_effect(inputs)
    if delayed and current.ke0 is None:
        current = replace(current, ke0=inputs.ke0)
    if use_form_start:
        start = current
    else:
        # the data-derived guess uses the effect-site concentration when ke0 is fixed
        guide = site(current.ke0) if (delayed and "ke0" not in free) else plasma
        guess = guess_effect_start(guide, obs.effect)
        keep = {n: getattr(current, n) for n in EFFECT_FITTABLE if n not in free}
        start = replace(guess, **keep)  # ke0: form value if fixed, None (scan) if free
    result = fit_effect(plasma, obs.effect, start=start, free=free, site=site)
    drive = site(result.effect.ke0) if result.effect.ke0 is not None else plasma
    return result, drive


def apply_effect_fit(inputs: ModelInputs, result: EffectFitResult) -> ModelInputs:
    """Copy the fitted effect parameters back into the form values and switch the effect on.

    A fit that included ke0 also switches the delayed-effect option on.
    """
    e = result.effect
    changes = dict(use_effect=True, e0=e.e0, emax=e.emax, ec50=e.ec50, hill=e.hill)
    if e.ke0 is not None:
        changes.update(delayed_effect=True, ke0=e.ke0)
    return replace(inputs, **changes)


def predict_effects(inputs: ModelInputs, obs: EffectObservations) -> tuple[np.ndarray, np.ndarray]:
    """(driving concentration, predicted effect) at the data times for the form's current model.

    The driving concentration is plasma for a direct effect and the effect-site
    concentration for a delayed one.
    """
    conc = driving_concentrations_at(inputs, obs)
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
    """Plain-language warnings to show next to an effect fit.

    `conc` is the concentration that drove the fit (what run_effect_fit returned).
    """
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
            if name in _LOG_SCALE_PARAMS:
                cv = result.cv_percent.get(name)
                if cv is not None and cv > _POOR_CV_PERCENT:
                    poor.append(f"{label} (CV {cv:.0f}%)")
            elif result.effect_range > 0 and result.se[name] > _POOR_SE_FRACTION * result.effect_range:
                poor.append(f"{label} (uncertainty above a quarter of the measured effect range)")
        if poor:
            notes.append("Poorly determined: " + "; ".join(poor) + ".")
    if "ke0" in result.at_bound:
        ke0 = result.effect.ke0
        if ke0 is not None and ke0 >= 1.0:
            notes.append(
                "ke0 finished on its upper limit: the data show no detectable delay between "
                "concentration and effect. A direct effect describes them as well."
            )
        else:
            notes.append(
                "ke0 finished on its lower limit: the delay is longer than these data can "
                "resolve, so ke0 is not determined."
            )
    other_bound = tuple(n for n in result.at_bound if n != "ke0")
    if other_bound:
        names = ", ".join(EFFECT_PARAMETER_INFO[n][0] for n in other_bound)
        notes.append(
            f"{names} finished on its allowed limit, so the data cannot determine it. "
            "Fix it at a known value instead of trusting the fitted one."
        )
    conc = np.asarray(conc, dtype=float)
    if "ec50" in result.free and conc.size and float(conc.max()) > 0:
        c_max = float(conc.max())
        positive = conc[conc > 0]
        ec50 = result.effect.ec50
        what = "effect-site concentration" if result.effect.ke0 is not None else "concentration"
        if ec50 > c_max:
            notes.append(
                f"EC50 ({ec50:.4g} mg/L) is above the highest {what} at the data times "
                f"({c_max:.4g} mg/L): the effect never approaches saturation, so Emax and EC50 "
                "are only weakly determined."
            )
        elif positive.size and ec50 < float(positive.min()):
            notes.append(
                f"EC50 ({ec50:.4g} mg/L) is below the lowest non-zero {what} "
                f"({float(positive.min()):.4g} mg/L): the effect is saturated everywhere "
                "that was measured."
            )
    ke0 = result.effect.ke0
    if ke0 is None:
        notes.append(
            "The effect is modelled as following the plasma concentration with no delay. A delayed "
            "effect shows as a systematic pattern in the residuals."
        )
    else:
        notes.append(
            f"Effect-compartment model: the effect follows an effect-site concentration that lags "
            f"plasma with an equilibration half-life of {equilibration_half_life(ke0):.3g} h "
            f"(ke0 {ke0:.3g} 1/h). ke0 and EC50 are only as reliable as the PK model that "
            "supplies the plasma concentrations: a wrong PK model can look like a delay."
        )
    return notes
