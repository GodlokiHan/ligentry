"""Residuals of the current effect model against measured effects (no Qt).

A residual is observed minus predicted effect, in whatever unit the effect data
use. The predicted effect comes from the form's PK model (which gives the
concentration at each effect time) and the form's effect parameters, whether or
not the effect box is ticked.

The model is a direct effect with no delay. A good fit scatters around zero; a
systematic curve (for example residuals that are positive on the way up and
negative on the way down) points to a delayed effect that this model cannot
describe.

The result is the same `Residuals` container used for concentrations: its fields
are just times, observed and predicted values, so the unit is the caller's.
"""
from __future__ import annotations

from ligentry.pkpd.effect_data import EffectObservations
from ligentry.pkpd.effect_workflow import predict_effects
from ligentry.pkpd.residuals import Residuals
from ligentry.pkpd.workflow import ModelInputs


def effect_residuals_for(inputs: ModelInputs, obs: EffectObservations) -> Residuals:
    """Predict the effect at the data times with the form's model and compare."""
    _conc, predicted = predict_effects(inputs, obs)
    return Residuals(times=obs.times.copy(), observed=obs.effect.copy(), predicted=predicted)
