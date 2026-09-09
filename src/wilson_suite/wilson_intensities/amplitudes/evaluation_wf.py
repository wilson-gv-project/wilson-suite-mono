import time
import numpy as np

from dataclasses import dataclass, field

from typing import TYPE_CHECKING, Any, Dict
if TYPE_CHECKING:
    from wilson_suite.wilson_main.workflow_abstractions import WilsonSimulation
    from wilson_suite.wilson_main.abstractions import VibAnaSetup, MolecularProperty
    from wilson_suite.wilson_main.spectrum_abstractions import SpecEvalSetup
    from wilson_suite.wilson_intensities.amplitudes.grid_manager_evaluator import GridRegion
    from wilson_suite.wilson_intensities.amplitudes.numerical_abstractions import CompiledTermGroup, NumericalResonanceMotif
    from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import SpectralWindow, ResLocGeoObject
    from wilson_suite.wilson_derive.response_terms import VibPerturbedTerm

from wilson_suite.wilson_intensities.amplitudes.vibene_differences import VibDiffCache
from wilson_suite.wilson_intensities.amplitudes.term_parts import PrecalculatedData, VibStatesData
from wilson_suite.wilson_intensities.amplitudes.term_parts import EvaluationDataAndConfigs, ParameterSet, ResonanceMotif

from wilson_suite.wilson_utils.unit_convertor import convNu2Ene

from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import SpectralFeature
from wilson_suite.wilson_intensities.amplitudes.grid_manager_evaluator import GridManager

from wilson_suite.wilson_intensities.amplitudes.evaluators import prepTermsForEval
from wilson_suite.wilson_intensities.amplitudes.evaluators import prepDataForEval
from wilson_suite.wilson_intensities.amplitudes.evaluators import process_resonance_motifs
from wilson_suite.wilson_intensities.amplitudes.evaluators import evaluate_terms_coeffs
from wilson_suite.wilson_intensities.amplitudes.full_amplitude_coeff import (precalculate_unique_coeff_parts, 
                                                                             identify_precalc_unique_coeff_parts)
from wilson_suite.wilson_intensities.amplitudes.evaluators import get_features_to_draw
from wilson_suite.wilson_experiment.indep_vars_and_axes import SpectralAxisSet
from wilson_suite.wilson_derive.term_var_translate import translate_magn_conditions_to_axisvars, translate_terms_to_axis_variables
from wilson_suite.wilson_main.abstractions import DataOriginInfo
from wilson_suite.wilson_main.abstractions import VibAnaSetup, MolPropsCollection

from contextlib import contextmanager

import logging
logger = logging.getLogger("wilson")

@dataclass
class EvaluationContext:
    """
    Execution metadata for a single EvaluationWorkflow run.
    This class is intentionally non-scientific.
    """
    verbose: bool = False

    # name -> elapsed seconds
    timing: dict[str, float] = field(default_factory=dict)

    # name of the step currently executing / last failed
    failed_at: str | None = None

    # optional: step name -> arbitrary object
    intermediates: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class EvaluationInputs:
    terms: Any
    number_of_modes: int
    props: list["MolecularProperty"]
    spec_eval_setup: 'SpecEvalSetup'
    vib_ana_setup: 'VibAnaSetup'
    pulse_polarization_vector: tuple[float, float, float]

"""
terms [props, max_state_lvl] -- from get_fully_enhanced_terms(VibExperiment, ...)
magn_conds - from VibExperiment
pulse_polarization_vector -- terms should have that as a parameter???
gamma
axis_choice

vib_ana_configs (harm/ahnarm - anharm_inhouse or not)

dyn_range --- render config
spec_window

data:
number of normal modes (total and choices here)
states and energies (all)

"""

# -w1 + w2 is always significantly > 0 ==> magn_conditions=((-1, 2),)
MagnConditions = tuple[tuple[int]]

@dataclass(frozen=True)
class DerivedTerms:
    """Expressions and the conditions they were derived under. Experiment coordinates."""
    terms: dict[int, dict[tuple, VibPerturbedTerm]]
    magn_conds: MagnConditions | None
    # because this may be derived with more info than just the terms
    available_axes: tuple[SpectralAxisSet, ...]

    def in_axes(self, axis_choice: SpectralAxisSet) -> 'TermsInAxes':
        if axis_choice not in self.available_axes:
            raise ValueError(f"{axis_choice} not available for these terms")
        return TermsInAxes(
            axis_choice=axis_choice,
            terms=translate_terms_to_axis_variables(self.terms, axis_choice),
            magn_conds=(None if self.magn_conds is None else
                        translate_magn_conditions_to_axisvars(self.magn_conds, axis_choice)),
        )


@dataclass(frozen=True)
class TermsInAxes:
    """Everything in axis variables. Ready to evaluate."""
    axis_choice: SpectralAxisSet
    terms: dict[int, dict[tuple, VibPerturbedTerm]]
    magn_conds: MagnConditions | None

    def need_what(self) -> set[str]: ...


@dataclass(frozen=True)
class RspFunEvalSetup:
    terms: TermsInAxes # attributes: magn_conds, axis_choices(all possible), axes
    # methods: input_vars -- what data is needed for evaluation (states vib ene, max_state_lvs, axis_choice, gamma, pulse_polarization_vector)
    # this defines experiment represented
    # here should be with axis_choice
    # terms post axis selection

    # how to eval
    gamma_cm1: float
    polarization: np.ndarray
    grid: SpectralGrid
    vib_ana: VibAnaConfig
    calc_config: DataOriginInfo
    mask_forbidden_region: bool = False

    @property
    def gamma_au(self) -> float:
        return convNu2Ene(self.gamma_cm1)
    
    @property
    def axis_choice(self):
        return self.terms.axis_choice

    @property
    def magn_conds(self):
        return self.terms.magn_conds

    def __post_init__(self):
        if self.mask_forbidden_region and self.terms.magn_conds is None:
            raise ValueError("mask_forbidden_region set but terms carry no magn_conds")

    def make_datarequest(self) -> dict[str, DataOriginInfo]:
        return {name: self.calc_config for name in self.terms.need_what()}


@dataclass(frozen=True)
class EvalData:
    """Everything obtained externally. No configuration."""
    mol_props: MolPropsCollection
    states: tuple
    eigenvals: np.ndarray
    eigenvecs: np.ndarray

def stage_prep_terms(terms):
    terms = prepTermsForEval(terms)
    return terms

def stage_prep_data(eval_data: EvalData, rsp_eval_setup: RspFunEvalSetup, include_states_list):
    vibdiff_cache = VibDiffCache()

    number_of_nmodes = len(eval_data.eigenvals)
    vib_data = VibStatesData(allstates=tuple(eval_data.states), 
                             harmonic_osc_states_labels=include_states_list,
                             number_of_nmodes=number_of_nmodes)
    
    data_configs = EvaluationDataAndConfigs(props_data=eval_data.mol_props,
                                            vibstates_data=vib_data,
                                            number_of_nmodes=number_of_nmodes,
                                            nm_inds_choices=include_states_list,
                                            pulse_polarization_vector=rsp_eval_setup.polarization,
                                            nc_sqrt_eigval=eval_data.eigenvals)

    return vibdiff_cache, vib_data, data_configs

def stage_process_resonances(terms_flat, vib_data, vibdiffs) -> tuple[dict, dict]:
    motif_locs, terms_for_motifs = process_resonance_motifs(terms_flat, vib_data, vibdiffs)

    return motif_locs, terms_for_motifs


def stage_precalculations(terms_flat, data_configs):
    need_precalc = identify_precalc_unique_coeff_parts(terms=terms_flat)
    precalculated = precalculate_unique_coeff_parts(
        need_to_precalc=need_precalc, data_and_configs=data_configs)
    return precalculated


def stage_term_coefficients(terms_flat, motif_locs, data_configs, precalc) -> dict:

    return evaluate_terms_coeffs(terms_flat, motif_locs, data_configs, precalc)


def stage_get_allfeats(motif_locs, terms_for_motifs, coefficients, gamma):
    """
    gamma = rsp_eval_setup.gamma_cm1
    """
    # lineshape_parameter here is goint to be a single float now and be the same(uniform) for all features
    features, zero_feats = get_features_to_draw(motif_res_loc=motif_locs, 
                                                terms_for_motifs=terms_for_motifs,
                                                term_coeffs_per_index=coefficients,
                                                lineshape_parameter=gamma)

    return features, zero_feats


def stage_dress_with_featboxes(features, dyn_range):
    """
    box_range_safety_margin, scale_wrt_max_intensity, minimum_box_padding

    """

    max_intensity_in_window = SpectralFeature.get_max_intensity_feat(features).get_intensity()
    min_intensity_in_window = max_intensity_in_window / dyn_range

    features = SpectralFeature.dress_these_with_boxes(features,
                                                      max_intensity_in_window, 
                                                      min_intensity_in_window,
                                                      box_range_safety_margin=box_range_safety_margin,
                                                      scale_wrt_max_intensity=scale_wrt_max_intensity,
                                                      minimum_box_padding=minimum_box_padding,
                                                                        )

    return features

def stage_filter_magn_conds(features, rsp_eval_setup):
    features = SpectralFeature.apply_magn_cond_filter(features,
                                                      magn_conditions=rsp_eval_setup.terms.magn_conds,
                                                      magn_conditions_margin=magn_conditions_margin)

    return features

def stage_place_in_specwindow(features, spec_window):
    spec_window = SpectralFeature.filter_to_spec_window(features, spec_window)
    if not spec_window.full_features:
        raise ValueError("This SpectralWindow does not contain any features. " \
        "Change the bounds of the window or use different terms.")

    return spec_window


def stage_make_grid_manager(spec_window, grid_resolution):
    grid_manager = GridManager(spec_window)
    grid_manager.make_fullgrid(grid_resolution)

    return grid_manager


def stage_make_regions(grid_manager):
    regions = grid_manager.create_regions()
    if not regions:
        raise ValueError("No regions were created")

    return regions

def stage_regions_results(regions, vib_data, vibdiffs):

    regions_results = evaluate_regions(regions, vib_data, vibdiffs, gamma)

    return regions_results

def stage_place_results(grid_manager, regions_results):
    grid_manager.place_results_into_grid(regions_results)

    return grid_manager.full_grid


class Pipeline:
    STAGES = [stage_prep_terms, stage_prep_data, 
              stage_process_resonances,
              stage_precalculations, stage_term_coefficients, 
              stage_get_allfeats,
              stage_dress_with_featboxes, 
              # stage_filter_magn_conds,
              stage_place_in_specwindow,
              stage_make_grid_manager,
              stage_make_regions,
              stage_regions_results,
              stage_place_results]

    def run(self, setup, data, *, upto=None, resume_from=None): ...


@dataclass(frozen=True)
class EvaluationInputsExtended:
    """
        # .ev_info.Gamma, .ev_info.Gamma_unit, 
		# .ev_info.dynamic_range, .ev_info.box_range_safety_margin, 
		# .ev_info.scale_wrt_max_intensity, .ev_info.minimum_box_padding, 
		# .ev_info.exp_magn_conditions, .ev_info.magn_conditions_margin, 
		# .ev_info.spectral_window, .ev_info.grid_resolution, 
		# vib_ana_setup.states, vib_ana_setup.include_list,
		# vib_ana_setup.number_of_modes, vib_ana_setup.nc_sqrt_eigval
		# experiment.polarization_avg_vector
		# system.Nnmodes
    """
    terms: Any
    number_of_modes: int
    props: list["MolecularProperty"]
    Gamma: float
    Gamma_unit: str
    dynamic_range: int
    box_range_safety_margin: float
    scale_wrt_max_intensity: bool
    minimum_box_padding: float
    exp_magn_conditions: bool
    magn_conditions_margin: float
    spectral_window: None
    grid_resolution: None
    states: list
    include_list: list
    number_of_modes: int
    nc_sqrt_eigval: dict
    pulse_polarization_vector: tuple[float, float, float]


@dataclass
class EvaluationArtifacts:
    terms: list = None
    vib_data: 'VibStatesData' = None
    vibdiff_cache: 'VibDiffCache' = None
    data_configs: 'EvaluationDataAndConfigs' = None
    motif_locs: dict['ResonanceMotif', 'ResLocGeoObject'] = None
    terms_for_motifs: dict['ResonanceMotif', list['VibPerturbedTerm']] = None
    need_precalc: dict[str, Any] = None
    precalculated: 'PrecalculatedData' = None
    coefficients: dict['VibPerturbedTerm', dict['ParameterSet', float]] = None
    features: list[SpectralFeature] = None
    zero_feats: list[SpectralFeature] = None
    spec_window: 'SpectralWindow' = None
    grid_manager: GridManager = None
    regions: list['GridRegion'] = None
    regions_results: Dict[str, np.ndarray] = None


def make_evaluation_inputs(
                            *,
                            simulation: "WilsonSimulation" = None,
                            terms=None,
                            number_of_modes: int = None,
                            props: list["MolecularProperty"] = None,
                            spec_eval_setup: 'SpecEvalSetup' = None,
                            vib_ana_setup: 'VibAnaSetup' = None,
                            pulse_polarization_vector=None,
                        ) -> EvaluationInputs:
    """
    Create EvaluationInputs either from a WilsonSimulation
    or from explicitly provided components.

    Exactly one source of truth must be used.
    """

    if simulation is not None:
        if any(x is not None for x in (
            terms, number_of_modes, props,
            spec_eval_setup, vib_ana_setup,
            pulse_polarization_vector,
        )):
            raise ValueError(
                "Provide either simulation OR explicit inputs, not both"
            )

        # Extract and normalize from simulation
        spec_eval_setup = simulation.spec_eval_setup 
        # MR: Here assuming that spectral axes were set, so changed to use translated terms
        terms = simulation.terms_in_axis_choice
        number_of_modes = simulation.system.Nnmodes
        props = simulation.props
        vib_ana_setup = simulation.vib_ana_setup
        pulse_polarization_vector = tuple(
            simulation.exp.polarization_avg_vector
        )

    # ---- validation of completeness ----

    missing = [
        name for name, value in {
            "terms": terms,
            "number_of_modes": number_of_modes,
            "props": props,
            "spec_eval_setup": spec_eval_setup,
            "vib_ana_setup": vib_ana_setup,
            "pulse_polarization_vector": pulse_polarization_vector,
        }.items()
        if value is None
    ]

    if missing:
        raise ValueError(
            f"Missing required inputs: {', '.join(missing)}"
        )

    # ---- light normalization / sanity checks ----

    if len(pulse_polarization_vector) != 3:
        raise ValueError(
            "pulse_polarization_vector must be length 3"
        )

    return EvaluationInputs(
        terms=terms,
        number_of_modes=number_of_modes,
        props=props,
        spec_eval_setup=spec_eval_setup,
        vib_ana_setup=vib_ana_setup,
        pulse_polarization_vector=tuple(pulse_polarization_vector),
    )


class EvaluationWorkflow:
    """
    A workflow that tracks steps and captures intermediates on error
    
    Can work with a WilsonSimulation in a prepared state (READY) or stanalone with provided necessary inputs.
    
    """
    def __init__(self, inputs: EvaluationInputs, parallel=None, verbose: bool = False):
        """
        ctx = EvaluationContext which would hold timing, failures, intermediates(?) saved during the run
        artifacts = EvaluationArtifacts holds the intermediate artifacts of the run which are used at other points of the run

        Parameters:
            inputs: EvaluationInputs - has ( terms, number_of_modes, props(with vals), 
                                             spec_eval_setup, vib_ana_setup, 
                                             pulse_polarization_vector(for orientational avrg) )

        run() method - a sequence of executed steps; intermediate results are needed for further steps and are saved in self.artifacts
        """
        self.ctx = EvaluationContext(verbose=verbose)
        self.artifacts = EvaluationArtifacts()

        self.inputs = inputs
        self.parallel = parallel


    def _validate_inputs(self):
        """
        WilsonSimulation at this point should have necessary data.

        """
        if not self.inputs.terms:
            raise ValueError("Non-empty 'terms' should be provided")
        if not self.inputs.number_of_modes:
            raise ValueError("'number_of_modes' should be provided")
        if not self.inputs.props:
            raise ValueError("Non-empty 'props'  should be provided")
        else:
            for p in self.inputs.props:
                if p.vals is None:
                    raise ValueError(f"Property {p.trivial_name} has vals=None")

        if not self.inputs.vib_ana_setup.isAllSet:
            raise ValueError("'vib_ana_setup' should be all set (vib_ana_setup.isAllSet)")
        # validate spec_eval_setup, pulse_polarization_vector, number_of_modes(?)

        if not self.inputs.pulse_polarization_vector or len(self.inputs.pulse_polarization_vector) != 3:
            raise ValueError("'pulse_polarization_vector' should be a length 3 vector") # -- is it true?

    @contextmanager
    def step(self, name):
        self.ctx.failed_at = name
        start = time.time()
        try:
            yield
        finally:
            self.ctx.timing[name] = time.time() - start


    def run_get_feats(self):
        """
        Workflow with these steps:
            prep_terms
            prep_data
            process_resonances
            term_coefficients
            all_features -- this? should be accessible on top of list of feats
        """
        self._validate_inputs()
        
        try:
            # Part 1: Preparation
            with self.step("prep_terms"):
                self.artifacts.terms = prepTermsForEval(self.inputs.terms)

            with self.step("prep_data"): # could be in data inputs
                _data, _cache, _configs = prepDataForEval(self.inputs.pulse_polarization_vector, 
                                                          self.inputs.vib_ana_setup, 
                                                          self.inputs.props)
                self.artifacts.vib_data = _data
                self.artifacts.vibdiff_cache = _cache
                self.artifacts.data_configs = _configs

            # self._save_checkpoint('Step1')  # Save checkpoint

            # Part 2: Process resonances and calculate coefficients
            # get resonances locations for all terms
            with self.step("process_resonances"):
                self.artifacts.motif_locs, self.artifacts.terms_for_motifs = process_resonance_motifs(self.artifacts.terms,
                                                                                                            self.artifacts.vib_data,
                                                                                                            self.artifacts.vibdiff_cache)
            with self.step("term_coefficients"):
                self.artifacts.need_precalc = identify_precalc_unique_coeff_parts(terms=self.artifacts.terms)
                self.artifacts.precalculated = precalculate_unique_coeff_parts(
                    need_to_precalc=self.artifacts.need_precalc, data_and_configs=self.artifacts.data_configs)
                self.artifacts.coefficients = evaluate_terms_coeffs(self.artifacts.terms,
                                                                       self.artifacts.motif_locs,
                                                                       self.artifacts.data_configs,
                                                                       self.artifacts.precalculated)

            # self._save_checkpoint('Step2')  # Save checkpoint

            # Part 3: Extract features and place them in the spectral window
            with self.step("all_features"):
                
                # maybe should check unit by the value here as well? or somewhere before taking unit flag
                if self.inputs.spec_eval_setup.ev_info.Gamma_unit == 'au':
                    gamma = convNu2Ene(self.inputs.spec_eval_setup.ev_info.Gamma, reverse=True)
                elif self.inputs.spec_eval_setup.ev_info.Gamma_unit == 'cm-1':
                    gamma = self.inputs.spec_eval_setup.ev_info.Gamma
                else:
                    raise ValueError('Gamma cannot be converted from the given unit to au')
                
                # lineshape_parameter here is goint to be a single float now and be the same(uniform) for all features
                self.artifacts.features, self.artifacts.zero_feats = get_features_to_draw(motif_res_loc=self.artifacts.motif_locs, 
                                                                  terms_for_motifs=self.artifacts.terms_for_motifs,
                                                                  term_coeffs_per_index=self.artifacts.coefficients,
                                                                  lineshape_parameter=gamma)
                # print('\nall_features step')
                # print(f' There are {len(self.artifacts.features)} features')

        except Exception as e:
            from wilson_suite.wilson_utils.serialization import pickle_this_to
            filename_pkl = 'eval_wf.pkl'
            pickle_this_to(self, filename_pkl)
            
            raise type(e)(
                f"Failed at '{self.ctx.failed_at}': {e} EvaluationWorkflow instanse was saved to `{filename_pkl}`."
            ) from e


    def run_specwindow_feats(self, features=None):
        """
        
        """
        # features_to_use will be used only in first step - further steps are chained
        if self.artifacts.features is None:
            if features is None:
                raise ValueError("This workflow does not have self.artifacts.features values nor input features!")
            features_to_use = features
        else:
            features_to_use = self.artifacts.features
        
        try:
            with self.step("dress_with_featboxes"):
                max_intensity_in_window = SpectralFeature.get_max_intensity_feat(features_to_use).get_intensity()
                min_intensity_in_window = max_intensity_in_window / self.inputs.spec_eval_setup.ev_info.dynamic_range

                self.artifacts.features = SpectralFeature.dress_these_with_boxes(features_to_use,
                                                                                 max_intensity_in_window, 
                                                                                 min_intensity_in_window,
                                                                                 box_range_safety_margin=
                                                                                 self.inputs.spec_eval_setup.ev_info.box_range_safety_margin,
                                                                                 scale_wrt_max_intensity=
                                                                                 self.inputs.spec_eval_setup.ev_info.scale_wrt_max_intensity,
                                                                                 minimum_box_padding=
                                                                                 self.inputs.spec_eval_setup.ev_info.minimum_box_padding,
                                                                                 )
                # print('\ndress_with_featboxes step')
                # print(f' There are {len(self.artifacts.features)} features')
                # SpectralFeature.print_list_features(self.artifacts.features)

            if self.inputs.spec_eval_setup.ev_info.apply_exp_magn_conditions_eval:
                with self.step("filter_magn_conds"):
                    self.artifacts.features = SpectralFeature.apply_magn_cond_filter(self.artifacts.features,
                                                                                    magn_conditions=self.inputs.spec_eval_setup.ev_info.exp_magn_conditions,
                                                                                    magn_conditions_margin=self.inputs.spec_eval_setup.ev_info.magn_conditions_margin)
                    # print('\nfilter_magn_conds step')
                    # print(f' There are {len(self.artifacts.features)} features')
                    # SpectralFeature.print_list_features(self.artifacts.features)

            with self.step("place_in_specwindow"):
                self.artifacts.spec_window = SpectralFeature.filter_to_spec_window(self.artifacts.features, self.inputs.spec_eval_setup.ev_info.spectral_window)
                if not self.artifacts.spec_window.full_features:
                    raise ValueError("This SpectralWindow does not contain any features. Change the bounds of the window or use different terms.")

        except Exception as e:
            from wilson_suite.wilson_utils.serialization import pickle_this_to
            filename_pkl = 'eval_wf.pkl'
            pickle_this_to(self, filename_pkl)
            
            raise type(e)(
                f"Failed at '{self.ctx.failed_at}': {e} EvaluationWorkflow instanse was saved to `{filename_pkl}`."
            ) from e
    

    def run(self, custom_grid=None, verbose=False):
        """
        Run evaluation, return dict with axes and results grid
        """
        self._validate_inputs()
        
        try:
            # Part 1: Preparation
            with self.step("prep_terms"):
                self.artifacts.terms = prepTermsForEval(self.inputs.terms)

            with self.step("prep_data"): # could be in data inputs
                _data, _cache, _configs = prepDataForEval(self.inputs.pulse_polarization_vector, 
                                                          self.inputs.vib_ana_setup, 
                                                          self.inputs.props)
                self.artifacts.vib_data = _data
                self.artifacts.vibdiff_cache = _cache
                self.artifacts.data_configs = _configs

            # self._save_checkpoint('Step1')  # Save checkpoint

            # Part 2: Process resonances and calculate coefficients
            # get resonances locations for all terms
            with self.step("process_resonances"):
                self.artifacts.motif_locs, self.artifacts.terms_for_motifs = process_resonance_motifs(self.artifacts.terms,
                                                                                                            self.artifacts.vib_data,
                                                                                                            self.artifacts.vibdiff_cache)
            with self.step("term_coefficients"):
                self.artifacts.need_precalc = identify_precalc_unique_coeff_parts(terms=self.artifacts.terms)
                self.artif8u7acts.precalculated = precalculate_unique_coeff_parts(
                    need_to_precalc=self.artifacts.need_precalc, data_and_configs=self.artifacts.data_configs)
                self.artifacts.coefficients = evaluate_terms_coeffs(self.artifacts.terms,
                                                                       self.artifacts.motif_locs,
                                                                       self.artifacts.data_configs,
                                                                       self.artifacts.precalculated)

            # self._save_checkpoint('Step2')  # Save checkpoint

            # Part 3: Extract features and place them in the spectral window
            with self.step("all_features"):
                
                # maybe should check unit by the value here as well? or somewhere before taking unit flag
                if self.inputs.spec_eval_setup.ev_info.Gamma_unit == 'au':
                    gamma = convNu2Ene(self.inputs.spec_eval_setup.ev_info.Gamma, reverse=True)
                elif self.inputs.spec_eval_setup.ev_info.Gamma_unit == 'cm-1':
                    gamma = self.inputs.spec_eval_setup.ev_info.Gamma
                else:
                    raise ValueError('Gamma cannot be converted from the given unit to au')
                
                # lineshape_parameter here is goint to be a single float now and be the same(uniform) for all features
                self.artifacts.features, self.artifacts.zero_feats = get_features_to_draw(motif_res_loc=self.artifacts.motif_locs, 
                                                                  terms_for_motifs=self.artifacts.terms_for_motifs,
                                                                  term_coeffs_per_index=self.artifacts.coefficients,
                                                                  lineshape_parameter=gamma)
                print('\nall_features step')
                print(f' There are {len(self.artifacts.features)} features')
                if verbose:
                    SpectralFeature.print_list_features(self.artifacts.features)


            with self.step("dress_with_featboxes"):
                max_intensity_in_window = SpectralFeature.get_max_intensity_feat(self.artifacts.features).get_intensity()
                min_intensity_in_window = max_intensity_in_window / self.inputs.spec_eval_setup.ev_info.dynamic_range

                self.artifacts.features = SpectralFeature.dress_these_with_boxes(self.artifacts.features,
                                                                                 max_intensity_in_window, 
                                                                                 min_intensity_in_window,
                                                                                 box_range_safety_margin=
                                                                                 self.inputs.spec_eval_setup.ev_info.box_range_safety_margin,
                                                                                 scale_wrt_max_intensity=
                                                                                 self.inputs.spec_eval_setup.ev_info.scale_wrt_max_intensity,
                                                                                 minimum_box_padding=
                                                                                 self.inputs.spec_eval_setup.ev_info.minimum_box_padding,
                                                                                 )
                print('\ndress_with_featboxes step')
                print(f' There are {len(self.artifacts.features)} features')
                if verbose:
                    SpectralFeature.print_list_features(self.artifacts.features)

            if self.inputs.spec_eval_setup.ev_info.apply_exp_magn_conditions_eval:
                with self.step("filter_magn_conds"):
                    self.artifacts.features = SpectralFeature.apply_magn_cond_filter(self.artifacts.features,
                                                                                    magn_conditions=self.inputs.spec_eval_setup.ev_info.exp_magn_conditions,
                                                                                    magn_conditions_margin=self.inputs.spec_eval_setup.ev_info.magn_conditions_margin)
                    print('\nfilter_magn_conds step')
                    print(f' There are {len(self.artifacts.features)} features')
                    if verbose:
                        SpectralFeature.print_list_features(self.artifacts.features)

            with self.step("place_in_specwindow"):
                self.artifacts.spec_window = SpectralFeature.filter_to_spec_window(self.artifacts.features, self.inputs.spec_eval_setup.ev_info.spectral_window)
                if not self.artifacts.spec_window.full_features:
                    raise ValueError("This SpectralWindow does not contain any features. Change the bounds of the window or use different terms.")
                
            # self._save_checkpoint('Step3')  # Save checkpoint

            # Part 4: Grid management and region evaluation
            with self.step("make_grid_manager"):
                self.artifacts.grid_manager = GridManager(self.artifacts.spec_window)
                self.artifacts.grid_manager.make_fullgrid(self.inputs.spec_eval_setup.ev_info.grid_resolution)

            with self.step("make_regions"):
                '''
                in create_regions:
                    - formal_domains = self.spec_window.find_clusters_by_featboxes():  --- wilson_intensities/amplitudes/grid_manager_evaluator.py
                        - clusters = domains.features_to_clusters(features=all_features) --- wilson_intensities/amplitudes/spectrum_composition.py
                        - features: list = clusters[c]
                        - tuple(RectangularDomain.from_features(features) for c in clusters) -- this returns (clusters are turned into RectangularDomains)
                
                [clusters are grouped features based on overalps of feature boxes here; so input features should be dressed with boxes]

                featured with default boxes from postinit are created in GridManager.spec_window.find_clusters_by_featboxes, 
                    but actaully features are initialized before - in get_features_to_draw(), in "all_features" step, 
                    and spectral window holds them.
                '''
                self.artifacts.regions = self.artifacts.grid_manager.create_regions()
                if not self.artifacts.regions:
                    raise ValueError("No regions were created")

            with self.step("regions_results"):
                if self.inputs.spec_eval_setup.ev_info.Gamma_unit == 'cm-1':
                    gamma = convNu2Ene(self.inputs.spec_eval_setup.ev_info.Gamma)
                elif self.inputs.spec_eval_setup.ev_info.Gamma_unit == 'au':
                    gamma = self.inputs.spec_eval_setup.ev_info.Gamma
                else:
                    raise ValueError('Gamma cannot be converted from the given unit to au')
                self.artifacts.regions_results = evaluate_regions(self.artifacts.regions, 
                                                                     self.artifacts.vib_data, 
                                                                     self.artifacts.vibdiff_cache,
                                                                     gamma,
                                                                     self.ctx.verbose)

            # self._save_checkpoint('Step4')  # Save checkpoint

            # Part 5: Assemble the full grid
            with self.step("place_results"):
                self.artifacts.grid_manager.place_results_into_grid(self.artifacts.regions_results)

            # Return results
            return self.artifacts.grid_manager.full_grid
            
        except Exception as e:
            from wilson_suite.wilson_utils.serialization import pickle_this_to
            filename_pkl = 'eval_wf.pkl'
            pickle_this_to(self, filename_pkl)
            
            raise type(e)(
                f"Failed at '{self.ctx.failed_at}': {e} EvaluationWorkflow instanse was saved to `{filename_pkl}`."
            ) from e
 
    def _save_checkpoint(self, name: str):
        """
        FIXME: self.results are not serializable yet
        """
        raise NotImplementedError('_save_checkpoint')
        import json

        # Save intermediate results to a file or log them
        with open(f'checkpoint_{name}.json', 'w') as f:
            json.dump(self.results, f)


class EvaluationWorkflow_NEW:
    """
    A workflow that tracks steps and captures intermediates on error
    
    Can work with a WilsonSimulation in a prepared state (READY) or stanalone with provided necessary inputs.
    
    """
    from wilson_suite.wilson_main.workflow_abstractions_updated import SealedSetup

    def __init__(self, 
                 setup_inputs: SealedSetup, 
                 parallel=None, verbose: bool = False):
        """
        ctx = EvaluationContext which would hold timing, failures, intermediates(?) saved during the run
        artifacts = EvaluationArtifacts holds the intermediate artifacts of the run which are used at other points of the run

        Parameters:
            inputs: EvaluationInputs - has ( terms, number_of_modes, props(with vals), 
                                             spec_eval_setup, vib_ana_setup, 
                                             pulse_polarization_vector(for orientational avrg) )

        run() method - a sequence of executed steps; intermediate results are needed for further steps and are saved in self.artifacts
        """
        self.ctx = EvaluationContext(verbose=verbose)
        self.artifacts = EvaluationArtifacts()

        self.setup_inputs = setup_inputs
        self.terms_flat = self.setup_inputs.terms_in_axes.make_flat(as_list=True)

        self.parallel = parallel


    @contextmanager
    def step(self, name):
        self.ctx.failed_at = name
        start = time.time()
        try:
            yield
        finally:
            self.ctx.timing[name] = time.time() - start



    def run(self, custom_grid=None, verbose=False):
        """
        Run evaluation, return dict with axes and results grid

        prep_data - vibana, props, polarization -- makes more user-friendly structures -- IS IT NEEDED?
        process_resonances - identify motif_locs and terms_for_motifs
        term_coefficients - needs terms and motif_locs
        all_features - 
        dress_with_featboxes - 
        filter_magn_conds - 
        place_in_specwindow - 
        make_grid_manager - 
        make_regions - 
        regions_results - 
        place_results - 
        """
        
        # Part 1: Preparation

        with self.step("prep_data"): # could be in data inputs
            self.artifacts.vibdiff_cache = VibDiffCache()
            self.artifacts.vib_data = VibStatesData(allstates=tuple(self.setup_inputs.vib_ana.states), 
                                        harmonic_osc_states_labels=self.setup_inputs.vib_ana.include_list,
                                        number_of_nmodes=self.setup_inputs.vib_ana.number_of_modes)
            self.artifacts.data_configs = EvaluationDataAndConfigs(props_data=self.setup_inputs.prop_order.props_coll,
                                                        vibstates_data=self.artifacts.vib_data,
                                                        number_of_nmodes=self.setup_inputs.vib_ana.number_of_modes,
                                                        nm_inds_choices=self.setup_inputs.vib_ana.include_list,
                                                        pulse_polarization_vector=self.setup_inputs.experiment.polarization_avg_vector,
                                                        nc_sqrt_eigval=self.setup_inputs.vib_ana.nc_sqrt_eigval)

        # self._save_checkpoint('Step1')  # Save checkpoint

        # Part 2: Process resonances and calculate coefficients
        # get resonances locations for all terms
        with self.step("process_resonances"):
            self.artifacts.motif_locs, self.artifacts.terms_for_motifs = process_resonance_motifs(self.terms_flat,
                                                                                                        self.artifacts.vib_data,
                                                                                                        self.artifacts.vibdiff_cache)
        with self.step("term_coefficients"):
            self.artifacts.need_precalc = identify_precalc_unique_coeff_parts(terms=self.terms_flat)
            self.artifacts.precalculated = precalculate_unique_coeff_parts(
                need_to_precalc=self.artifacts.need_precalc, data_and_configs=self.artifacts.data_configs)
            self.artifacts.coefficients = evaluate_terms_coeffs(self.terms_flat,
                                                                    self.artifacts.motif_locs,
                                                                    self.artifacts.data_configs,
                                                                    self.artifacts.precalculated)

        # self._save_checkpoint('Step2')  # Save checkpoint

        # Part 3: Extract features and place them in the spectral window
        with self.step("all_features"):
            
            # maybe should check unit by the value here as well? or somewhere before taking unit flag
            if self.setup_inputs.spec_eval.ev_info.Gamma_unit == 'au':
                gamma = convNu2Ene(self.setup_inputs.spec_eval.ev_info.Gamma, reverse=True)
            elif self.setup_inputs.spec_eval.ev_info.Gamma_unit == 'cm-1':
                gamma = self.setup_inputs.spec_eval.ev_info.Gamma
            else:
                raise ValueError('Gamma cannot be converted from the given unit to au')
            
            # lineshape_parameter here is goint to be a single float now and be the same(uniform) for all features
            self.artifacts.features, self.artifacts.zero_feats = get_features_to_draw(motif_res_loc=self.artifacts.motif_locs, 
                                                                terms_for_motifs=self.artifacts.terms_for_motifs,
                                                                term_coeffs_per_index=self.artifacts.coefficients,
                                                                lineshape_parameter=gamma)
            print('\nall_features step')
            print(f' There are {len(self.artifacts.features)} features')
            if verbose:
                SpectralFeature.print_list_features(self.artifacts.features)

        with self.step("dress_with_featboxes"):
            max_intensity_in_window = SpectralFeature.get_max_intensity_feat(self.artifacts.features).get_intensity()
            min_intensity_in_window = max_intensity_in_window / self.setup_inputs.spec_eval.ev_info.dynamic_range

            self.artifacts.features = SpectralFeature.dress_these_with_boxes(self.artifacts.features,
                                                                                max_intensity_in_window, 
                                                                                min_intensity_in_window,
                                                                                box_range_safety_margin=
                                                                                self.setup_inputs.spec_eval.ev_info.box_range_safety_margin,
                                                                                scale_wrt_max_intensity=
                                                                                self.setup_inputs.spec_eval.ev_info.scale_wrt_max_intensity,
                                                                                minimum_box_padding=
                                                                                self.setup_inputs.spec_eval.ev_info.minimum_box_padding,
                                                                                )
            print('\ndress_with_featboxes step')
            print(f' There are {len(self.artifacts.features)} features')
            if verbose:
                SpectralFeature.print_list_features(self.artifacts.features)

        if self.setup_inputs.spec_eval.ev_info.apply_exp_magn_conditions_eval:
            with self.step("filter_magn_conds"):
                self.artifacts.features = SpectralFeature.apply_magn_cond_filter(self.artifacts.features,
                                                                                magn_conditions=self.setup_inputs.terms_in_axes.magn_conditions,
                                                                                magn_conditions_margin=self.setup_inputs.spec_eval.ev_info.magn_conditions_margin)
                print('\nfilter_magn_conds step')
                print(f' There are {len(self.artifacts.features)} features')
                if verbose:
                    SpectralFeature.print_list_features(self.artifacts.features)

        with self.step("place_in_specwindow"):
            self.artifacts.spec_window = SpectralFeature.filter_to_spec_window(self.artifacts.features, self.setup_inputs.spec_eval.ev_info.spectral_window)
            if not self.artifacts.spec_window.full_features:
                raise ValueError("This SpectralWindow does not contain any features. Change the bounds of the window or use different terms.")
            
        # self._save_checkpoint('Step3')  # Save checkpoint

        # Part 4: Grid management and region evaluation
        with self.step("make_grid_manager"):
            self.artifacts.grid_manager = GridManager(self.artifacts.spec_window)
            self.artifacts.grid_manager.make_fullgrid(self.setup_inputs.spec_eval.ev_info.grid_resolution)

        with self.step("make_regions"):
            '''
            in create_regions:
                - formal_domains = self.spec_window.find_clusters_by_featboxes():  --- wilson_intensities/amplitudes/grid_manager_evaluator.py
                    - clusters = domains.features_to_clusters(features=all_features) --- wilson_intensities/amplitudes/spectrum_composition.py
                    - features: list = clusters[c]
                    - tuple(RectangularDomain.from_features(features) for c in clusters) -- this returns (clusters are turned into RectangularDomains)
            
            [clusters are grouped features based on overalps of feature boxes here; so input features should be dressed with boxes]

            featured with default boxes from postinit are created in GridManager.spec_window.find_clusters_by_featboxes, 
                but actaully features are initialized before - in get_features_to_draw(), in "all_features" step, 
                and spectral window holds them.
            '''
            self.artifacts.regions = self.artifacts.grid_manager.create_regions()
            if not self.artifacts.regions:
                raise ValueError("No regions were created")

        with self.step("regions_results"):
            if self.setup_inputs.spec_eval.ev_info.Gamma_unit == 'cm-1':
                gamma = convNu2Ene(self.setup_inputs.spec_eval.ev_info.Gamma)
            elif self.setup_inputs.spec_eval.ev_info.Gamma_unit == 'au':
                gamma = self.setup_inputs.spec_eval.ev_info.Gamma
            else:
                raise ValueError('Gamma cannot be converted from the given unit to au')
            self.artifacts.regions_results = evaluate_regions(self.artifacts.regions, 
                                                                    self.artifacts.vib_data, 
                                                                    self.artifacts.vibdiff_cache,
                                                                    gamma,
                                                                    self.ctx.verbose)

        # self._save_checkpoint('Step4')  # Save checkpoint

        # Part 5: Assemble the full grid
        with self.step("place_results"):
            self.artifacts.grid_manager.place_results_into_grid(self.artifacts.regions_results)

        # Return results
        return self.artifacts.grid_manager.full_grid
        
 
    def _save_checkpoint(self, name: str):
        """
        FIXME: self.results are not serializable yet
        """
        raise NotImplementedError('_save_checkpoint')
        import json

        # Save intermediate results to a file or log them
        with open(f'checkpoint_{name}.json', 'w') as f:
            json.dump(self.results, f)




def evaluate_regions(regions: list["GridRegion"], 
                     vib_data: "VibStatesData", 
                     vibdiff_cache: "VibDiffCache", 
                     gamma: float):

    # Step 2: Evaluate each region
    region_results = {}
    for region in regions:
        
        region_results[region] = evaluate_region(region, vib_data, vibdiff_cache, gamma)
        

    return region_results

def evaluate_region(region: "GridRegion",
                    vib_data: "VibStatesData", 
                    vibdiff_cache: "VibDiffCache", 
                    gamma: float) -> np.ndarray:
    """Evaluate all features in a single grid region."""
    # Initialize result array
    target_shape = np.broadcast(*(arr for arr in region.coords.values())).shape
    result = np.zeros(target_shape, dtype=complex)
    
    # Sum contributions from all features
    for feature in region.features:
        
        result += evaluate_feature(feature, vib_data, vibdiff_cache, gamma, region.coords_au)
        
    return result

def evaluate_feature(feature: 'SpectralFeature', 
                     vib_data: "VibStatesData", 
                     vibdiff_cache: "VibDiffCache", 
                     gamma: float,
                     coords: dict[str, np.ndarray]) -> np.ndarray:
    """Evaluate a single feature on grid coordinates."""
    # Compile feature to numerical form
    from .numerical_abstractions import compile_feature
    compiled_groups = compile_feature(feature, vib_data, vibdiff_cache)
    
    # Sum all compiled groups
    target_shape = np.broadcast(*(arr for arr in coords.values())).shape
    feature_sum = np.zeros(target_shape, dtype=complex)

    
    for group in compiled_groups:
        feature_sum += evaluate_compiled_group(group, coords, gamma)
    
    # Apply amplitude coefficient
    return feature.amplitude_coeff * feature_sum


def evaluate_resonance_motif(motif: 'NumericalResonanceMotif',
                             coords: Dict[str, np.ndarray],
                             gamma: float) -> np.ndarray:
    """
    Calculate resonance motif contribution at grid points.
    
    Args:
        motif: Compiled resonance motif with conditions
        coords: Dict of axis_label -> meshgrid array
        
    Returns:
        Complex array with resonance contributions
    """
    target_shape = np.broadcast(*(arr for arr in coords.values())).shape
    total = np.ones(target_shape, dtype=complex)
    
    for res_cond in motif.res_conds:
        # Calculate photon frequency: sum over axes
        pfreq = sum(coords[ax] * res_cond.pf_dict[ax] 
                    for ax in res_cond.pf_dict)
        # Resonance denominator
        z = res_cond.vib_energy_diff - pfreq - 1j * gamma
        total *= 1.0 / z
        
    return total

def evaluate_compiled_group(group: 'CompiledTermGroup',
                            coords: Dict[str, np.ndarray],
                            gamma: float) -> np.ndarray:
    """Sum all resonance motifs in a compiled group."""
    target_shape = np.broadcast(*(arr for arr in coords.values())).shape
    result = np.zeros(target_shape, dtype=complex)
    # result = np.zeros_like(next(iter(coords.values())), dtype=complex)
    
    for motif in group.resonance_motifs:
        result += evaluate_resonance_motif(motif, coords, gamma)
    # print('evaluate_compiled_group', result)
    return result

