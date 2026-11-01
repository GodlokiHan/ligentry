"""Residuals of the current model against measured concentrations (no Qt).

A residual is observed minus predicted concentration, in mg/L. The PK/PD window
plots them under the concentration curve so a poor fit is visible at a glance: a
good fit scatters around zero, a systematic curve in the residuals means the
model (or the dosing entered) does not match the data.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ligentry.pkpd.model import simulate
from ligentry.pkpd.observations import Observations
from ligentry.pkpd.workflow import ModelInputs, build_model


@dataclass(frozen=True)
class Residuals:
    times: np.ndarray  # h
    observed: np.ndarray  # mg/L
    predicted: np.ndarray  # mg/L

    @property
    def residual(self) -> np.ndarray:
        return self.observed - self.predicted

    @property
    def rmse(self) -> float:
        return float(np.sqrt(np.mean(self.residual**2)))

    def __len__(self) -> int:
        return int(self.times.size)


def residuals_for(inputs: ModelInputs, obs: Observations) -> Residuals:
    """Simulate the model at the observation times and compare with the data."""
    pk, regimen, _effect = build_model(inputs)
    predicted = simulate(pk, regimen, obs.times).conc
    return Residuals(times=obs.times.copy(), observed=obs.conc.copy(), predicted=predicted)
