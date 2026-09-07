"""
===
Setup types and the two public entry points for response-function evaluation.
===

rsp_evaluator


"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from wilson_suite.wilson_experiment.indep_vars_and_axes import SpectralAxisSet
from wilson_suite.wilson_derive.response_terms import VibPerturbedTerm
from wilson_suite.wilson_main.abstractions import MolPropsCollection
from wilson_suite.wilson_utils.unit_convertor import convNu2Ene

from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import SpectralFeature
from wilson_suite.wilson_intensities.amplitudes.vibene_differences import VibDiffCache
from wilson_suite.wilson_intensities.amplitudes.term_parts import (
    EvaluationDataAndConfigs,
    VibStatesData,
)
from wilson_suite.wilson_intensities.amplitudes.evaluators import (
    process_resonance_motifs,
    evaluate_terms_coeffs,
    get_features_to_draw,
)
from wilson_suite.wilson_intensities.amplitudes.full_amplitude_coeff import (
    precalculate_unique_coeff_parts,
    identify_precalc_unique_coeff_parts,
)

# FIXME: these three were used but never imported in the original module.
# Point these at their real locations.
from wilson_suite.wilson_intensities.amplitudes.grid_manager_evaluator import GridManager
# from wilson_suite.wilson_intensities.amplitudes. EvaluatedResult
from wilson_suite.wilson_intensities.amplitudes.evaluation_wf import evaluate_region as _eval_one_region


# -w1 + w2 is always > 0 ==> magn_conds = ((-1, 2),)
MagnConditions = tuple[tuple[int, ...], ...]


# --- configuration -------------------------------------------------------

"""
TermsInAxes.need_what() (rps_evaluation.py:65, currently NotImplementedError) is the work manifest. It enables "verify all data present before computing anything", and batching of QC jobs from the manifest — which MolPropsCollection.group_by_calc_setup / build_request_dict are approaching from the wrong end.
"""

@dataclass(frozen=True)
class TermsInAxes:
    """Terms expressed in axis variables. Caller is responsible for flattening."""

    axis_choice: SpectralAxisSet
    terms: list[VibPerturbedTerm]
    magn_conds: MagnConditions | None = None

    def __post_init__(self) -> None:
        bad = [t for t in self.terms if not isinstance(t, VibPerturbedTerm)]
        if bad:
            raise TypeError(f"TermsInAxes.terms holds {len(bad)} non-VibPerturbedTerm entries")

    def need_what(self) -> set[str]:
        raise NotImplementedError


@dataclass(frozen=True)
class SpectralGrid:
    window: Any  # spectral window bounds
    resolution: float


@dataclass(frozen=True)
class BoxParams:
    dyn_range: float
    box_range_safety_margin: float
    minimum_box_padding: float
    scale_wrt_max_intensity: bool = True


@dataclass(frozen=True)
class RspFunEvalSetup:
    terms: TermsInAxes
    gamma_cm1: float
    polarization: np.ndarray
    include_states: tuple
    grid: SpectralGrid
    boxes: BoxParams
    mask_forbidden_region: bool = False
    magn_conds_margin: float = 0.0

    def __post_init__(self) -> None:
        if self.mask_forbidden_region and self.terms.magn_conds is None:
            raise ValueError("mask_forbidden_region set but terms carry no magn_conds")

    @property
    def gamma_au(self) -> float:
        return convNu2Ene(self.gamma_cm1)

    @property
    def axis_choice(self) -> SpectralAxisSet:
        return self.terms.axis_choice

    @property
    def magn_conds(self) -> MagnConditions | None:
        return self.terms.magn_conds


# --- external data -------------------------------------------------------


@dataclass(frozen=True)
class MolSystemData:
    """Everything obtained externally. No configuration."""

    name: str
    states: tuple
    eigenvals: np.ndarray
    eigenvecs: np.ndarray
    mol_props: MolPropsCollection
    natoms: int | None = None
    geo: Any = None
    geo_extra: Any = None
    linear: bool = False
    conformer: str = "conf1"


# --- public API ----------------------------------------------------------


def compute_features(setup: RspFunEvalSetup, data: MolSystemData) -> list[SpectralFeature]:
    """
    Locate resonance motifs, evaluate term coefficients, build features.
    
    1. terms for eval
    2. data for eval
    3. res motifs - finding locations
    4. precalculate parts for term coeffs
    5. coeffs
    6. features from locs and coeffs


    """
    terms = setup.terms.terms
    configs = _prep_data(data, setup)

    motif_locs, terms_for_motifs = process_resonance_motifs(
        terms, configs.vibstates_data, VibDiffCache()
    )

    need_precalc = identify_precalc_unique_coeff_parts(terms=terms)
    precalc = precalculate_unique_coeff_parts(
        need_to_precalc=need_precalc, data_and_configs=configs
    )
    coeffs = evaluate_terms_coeffs(terms, motif_locs, configs, precalc)

    features, _zero = get_features_to_draw(
        motif_locs, terms_for_motifs, coeffs, setup.gamma_cm1
    )
    return features


def render_grid(
    features: list[SpectralFeature],
    setup: RspFunEvalSetup,
    executor=None,
) -> EvaluatedResult:
    """
    Box, filter, and evaluate features onto the spectral grid.
    
    1. 
    2. 

    """
    features = _dress_with_boxes(features, setup.boxes)

    if setup.mask_forbidden_region:
        features = SpectralFeature.apply_magn_cond_filter(
            features,
            magn_conditions=setup.magn_conds,
            magn_conditions_margin=setup.magn_conds_margin,
        )

    window = SpectralFeature.filter_to_spec_window(features, setup.grid.window)
    if not window.full_features:
        raise ValueError(
            "This SpectralWindow contains no features. Change the window bounds or the terms."
        )

    gm = GridManager(window)
    gm.make_fullgrid(setup.grid.resolution)

    regions = gm.create_regions()
    if not regions:
        raise ValueError("No regions were created")

    gm.place_results_into_grid(_evaluate_regions(regions, setup.gamma_au, executor))
    return EvaluatedResult(spec=gm.full_grid, axes=gm.axes, setup=setup)


# --- internals -----------------------------------------------------------


def _prep_data(data: MolSystemData, setup: RspFunEvalSetup) -> EvaluationDataAndConfigs:
    n_modes = len(data.eigenvals)
    vib_data = VibStatesData(
        allstates=tuple(data.states),
        harmonic_osc_states_labels=setup.include_states,
        number_of_nmodes=n_modes,
    )
    return EvaluationDataAndConfigs(
        props_data=data.mol_props,
        vibstates_data=vib_data,
        number_of_nmodes=n_modes,
        nm_inds_choices=setup.include_states,
        pulse_polarization_vector=setup.polarization,
        nc_sqrt_eigval=data.eigenvals,
    )


def _dress_with_boxes(
    features: list[SpectralFeature], boxes: BoxParams
) -> list[SpectralFeature]:
    max_intensity = SpectralFeature.get_max_intensity_feat(features).get_intensity()
    return SpectralFeature.dress_these_with_boxes(
        features,
        max_intensity,
        max_intensity / boxes.dyn_range,
        box_range_safety_margin=boxes.box_range_safety_margin,
        scale_wrt_max_intensity=boxes.scale_wrt_max_intensity,
        minimum_box_padding=boxes.minimum_box_padding,
    )


def _evaluate_regions(regions, 
                    vib_data: "VibStatesData", 
                    vibdiff_cache: "VibDiffCache", 
                      gamma_au: float, executor=None) -> list:
    if executor is None:
        return [_eval_one_region(r, vib_data, vibdiff_cache, gamma_au) for r in regions]
    return list(executor.map(_eval_one_region, regions, [gamma_au] * len(regions)))