"""
Evaluator functions for WilsonSimulation
"""
from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import ResLocGeoObject, SpectralFeature
from ..amplitudes.full_amplitude_coeff import evaluate_term_coeffs
from ..amplitudes.resonances import find_resonance_locations_wrt_index_choices, identify_unique_resmotifs

from wilson_suite.wilson_main.abstractions import VibAnaSetup, MolecularSystem, MolPropsCollection
from wilson_suite.wilson_intensities.amplitudes.term_parts import ResonanceMotif, VibStatesData
from wilson_suite.wilson_intensities.amplitudes.vibene_differences import VibDiffCache
from wilson_suite.wilson_intensities.amplitudes.term_parts import (EvaluationDataAndConfigs, ParameterSet,
                                                                   TermParametersChoice, PrecalculatedData)
from ...wilson_derive.response_terms import VibPerturbedTerm

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from wilson_suite.wilson_main.abstractions import VibAnaSetup, MolecularProperty
    from wilson_suite.wilson_main.spectrum_abstractions import SpecEvalSetup
    from wilson_suite.wilson_experiment.experiment_abstractions import VibExperiment

import numpy as np

import logging
logger = logging.getLogger("wilson."+__name__)

def prepTermsForEval(terms: dict | list) -> list:
    """
    put data in a form for use on the evaluation step
    """
    if isinstance(terms, type([])):
        for t in terms:
            if not isinstance(t, VibPerturbedTerm):
                raise ValueError("Smth that is not a VibPerturbedTerm was given in a list to prepTermsForEval()")
        return terms

    if isinstance(terms, type({})):
        for t_key in terms:
            if isinstance(terms[t_key], type({})):

                terms_as_list = []
                for i in terms:
                    for j in terms[i]:
                        for t in terms[i][j]:
                            terms_as_list.append(t)
                return terms_as_list
            else:
                if not isinstance(terms[t_key], VibPerturbedTerm):
                    raise ValueError("A flat dictionary but has smth other than VibPerturbedTerm as a value")
                return list(terms.values())

def prepDataForEval(pulse_polarization_vector: np.ndarray,
                    vib_ana_setup: 'VibAnaSetup',
                    props: list['MolecularProperty']) -> tuple[VibStatesData, VibDiffCache, EvaluationDataAndConfigs]:
    """
    put data in a form for use on the evaluation step
    """
    include_list = tuple([v for v in list(vib_ana_setup.nc_sqrt_eigval.keys()) if v not in vib_ana_setup.exclude_modes])
    if include_list == tuple():
        raise ValueError("include_list of included normal modes labels is empty")
    

    vibstates_data = VibStatesData(allstates=tuple(vib_ana_setup.states), 
                                   harmonic_osc_states_labels=vib_ana_setup.include_list,
                                   number_of_nmodes=vib_ana_setup.number_of_modes)
    vibdiff_cache = VibDiffCache()
    props = MolPropsCollection(properties=props)
    
    data_and_configs = EvaluationDataAndConfigs(props_data=props,
                                                vibstates_data=vibstates_data,
                                                number_of_nmodes=vib_ana_setup.number_of_modes,
                                                nm_inds_choices=include_list,
                                                pulse_polarization_vector=pulse_polarization_vector,
                                                nc_sqrt_eigval=vib_ana_setup.nc_sqrt_eigval)

    return vibstates_data, vibdiff_cache, data_and_configs


def process_resonance_motifs(derived_terms: list['VibPerturbedTerm'],
                            vibstates_data: VibStatesData,
                            vibdiff_cache: VibDiffCache) -> tuple[dict[ResonanceMotif,dict[ResLocGeoObject,list]], 
                                                                  dict[ResonanceMotif, list['VibPerturbedTerm']]]:
    """
    Take a list of terms and process resonance motifs and find their locations.
    """
    unique_res_motifs = identify_unique_resmotifs(derived_terms)
    motif_res_loc: dict[ResonanceMotif,dict[ResLocGeoObject,list]] = {}
    
    for res_motif in unique_res_motifs:
        this_motif_res_locs = find_resonance_locations_wrt_index_choices(
            motif=res_motif,
            vibstates_data=vibstates_data,
            vibdiff_cache=vibdiff_cache,
            spec_window=None
        )
        motif_res_loc.update(this_motif_res_locs)

    terms_for_motifs: dict[ResonanceMotif, list[VibPerturbedTerm]] = {}

    for vibterm in derived_terms:
        terms_for_motifs.setdefault(ResonanceMotif(vibterm.res), []).append(vibterm)

    return motif_res_loc, terms_for_motifs

def evaluate_terms_coeffs(derived_terms: list['VibPerturbedTerm'],
                  motif_res_loc: dict[ResonanceMotif, dict[ResLocGeoObject]],
                  data_and_configs: EvaluationDataAndConfigs,
                  precalculated: PrecalculatedData) -> dict['VibPerturbedTerm', dict[ParameterSet, float]]:
    """
    Evaluate coefficients for all terms.

    terms to dict with (params to coeff)
    """
    term_coeffs_per_index = {}
    for vibterm in derived_terms:
        res_motif = ResonanceMotif(vibterm.res)
        term_coeffs_per_index[vibterm] = evaluate_term_coeffs(
            term=vibterm,
            relevant_indices=[k for i in motif_res_loc[res_motif].values() for k in i],
            necessary_data=(data_and_configs, precalculated)
        )
    
    return term_coeffs_per_index

def evaluate_coeff_for_feat(feature: SpectralFeature, 
                            terms_hash_map: dict[int, 'VibPerturbedTerm'],
                            data_and_configs: EvaluationDataAndConfigs,
                            precalculated: PrecalculatedData) -> dict['VibPerturbedTerm', dict[ParameterSet, float]]:
    """
    Using feature information to evaluate its coefficient.

    A wrapper around evaluate_terms_coeffs()
    """
    derived_terms = [terms_hash_map[id] for id in feature.term_contributions[0].term_ids]
    assert feature.term_contributions[0].res_motif == ResonanceMotif(derived_terms[0].res)

    motif_res_loc = {feature.term_contributions[0].res_motif: 
                     {feature.location: [dict(p._parameters) for p in feature.term_contributions[0].states_parameters]}}

    term_coeffs_per_index = evaluate_terms_coeffs(derived_terms=derived_terms,
                                                  motif_res_loc=motif_res_loc,
                                                  data_and_configs=data_and_configs,
                                                  precalculated=precalculated)

    return term_coeffs_per_index


def get_features_from_terms_for_eval(derived_terms: list['VibPerturbedTerm'],
                                     vibstates_data: VibStatesData,
                                     vibdiff_cache: VibDiffCache,
                                     lineshape_parameter: float | None=None) -> list[SpectralFeature]:
    """
    SpectralFeature:
        location
        term_contributions = None
        lineshape_parameter = None
        amplitude_coeff = None
        feat_type: str = None - updated when sorted
        feat_box: Box = None - post-init if lineshape_parameter
    """
    motif_res_loc, terms_for_motifs = process_resonance_motifs(derived_terms, vibstates_data, vibdiff_cache)

    return get_features_to_draw(motif_res_loc, terms_for_motifs, lineshape_parameter=lineshape_parameter)[0]

def get_features_to_draw(motif_res_loc: dict[ResonanceMotif, dict[ResLocGeoObject, list]],
                         terms_for_motifs: dict[ResonanceMotif, list['VibPerturbedTerm']], 
                         term_coeffs_per_index: dict['VibPerturbedTerm', 
                                                     dict[ParameterSet, tuple[float, dict]]] | None=None,
                         lineshape_parameter: float | None=None) -> tuple[list[SpectralFeature], list[SpectralFeature]]:
    """

    lineshape_parameter - uniform lineshape parameters (for all axes) for each feature for now
    lineshape_parameter unit -- will be au - follows from the workflow in step("all_features")
    """
    # a SpectralFeature instanse holds a res_location and list of states parameters that give this res_location; 
    #       the amplitude coefficient is a value in the dict
    features_to_draw: list[SpectralFeature] = []
    zero_coeff_feats: list[SpectralFeature] = []

    for res_motif in motif_res_loc:

        for res_geo_obj, list_state_dicts in motif_res_loc[res_motif].items():

            lst_params = tuple(ParameterSet(states_dict) for states_dict in list_state_dicts)
            
            term_contributions=(
                                TermParametersChoice(
                                    res_motif=res_motif,
                                    states_parameters=lst_params,
                                    term_ids=tuple(t.to_str() for t in terms_for_motifs[res_motif]),
                                ),
                            )

            if term_coeffs_per_index is not None:
                list_to_sum = [term_coeffs_per_index[term][ParameterSet(states_dict)][0] for term in terms_for_motifs[res_motif] for states_dict in list_state_dicts]
                dict_of_contribs = {term.to_str(): term_coeffs_per_index[term][ParameterSet(states_dict)] for term in terms_for_motifs[res_motif] for states_dict in list_state_dicts}
                amplitude_coeff = sum(list_to_sum)
            else:
                amplitude_coeff = None
            
            # disregard locations where coefficient is zero
            spec_feature = SpectralFeature(location=res_geo_obj, 
                                        term_contributions=term_contributions,
                                        term_contrib_by_id = dict_of_contribs,
                                        lineshape_parameter=lineshape_parameter, # uniform lineshape parameters (for all axes) for each feature
                                        amplitude_coeff=amplitude_coeff)
            if amplitude_coeff != 0.:

                if spec_feature not in features_to_draw:
                    features_to_draw.append(spec_feature)
                else:
                    new_specfeat = spec_feature.union(features_to_draw[res_geo_obj][1])
                    features_to_draw.append(new_specfeat)
            else:
                zero_coeff_feats.append(spec_feature)

    return features_to_draw, zero_coeff_feats

