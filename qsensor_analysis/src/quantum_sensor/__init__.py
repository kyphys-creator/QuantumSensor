"""quantum_sensor: dark-matter quantum-sensor response analysis."""

from .analysis import DarkMatterQuantumAnalysis
from .config import RunConfig, BackgroundModel, BACKGROUND_SCENARIOS
from .statistics import (fit_toys, flux_band, asimov_significance,
                         mc_significance, ToyEnsemble,
                         find_confidence_band, pointwise_band,
                         save_pointwise_band, load_pointwise_band)

__all__ = [
    "DarkMatterQuantumAnalysis",
    "RunConfig",
    "BackgroundModel",
    "BACKGROUND_SCENARIOS",
    "fit_toys",
    "flux_band",
    "asimov_significance",
    "mc_significance",
    "ToyEnsemble",
    "find_confidence_band",
    "pointwise_band",
    "save_pointwise_band",
    "load_pointwise_band",
]
