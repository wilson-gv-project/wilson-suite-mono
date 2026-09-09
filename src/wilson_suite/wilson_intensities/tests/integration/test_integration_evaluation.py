import wilson_suite as ws
import numpy as np


def test_evaluation_general_customdata_1elterm():
    print()
    from ..unit.test_domains import get_data_evaluators_tests
    datadict = get_data_evaluators_tests()
    # 'system', 'vib_ana_setup', 'derived_terms', 'props', 
    # 'experiment', 'spec_eval_setup', 'domain_distance_thresholds'
    
    np.set_printoptions(linewidth=180, precision=3)

    from wilson_suite.wilson_main.workflow_abstractions import WilsonSimulation

    from wilson_suite.wilson_derive.derive import get_fully_enhanced_terms
    from ....fixtures import evv_experiment

    evv_exp = evv_experiment()

    axes_choice = evv_exp.valid_axis_combs[0].valid_axis_combs[3] # {'A': [(2,)], 'B': [(-1,), (2,)]}

    bounds_dict = {'B': (900., 900.), 'A': (1864., 1864.)}
    from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import SpectralWindow, Box
    spectral_window = SpectralWindow(box=Box(bounds_dict))

    evi = ws.main.spectrum_abstractions.EvaluationInfo(**{'spectral_window': spectral_window,
                                                          'Gamma': 1., 'Gamma_unit': 'cm-1',
                                                          'grid_resolution': {'A': 1, 'B': 1}})
    mock_sim = WilsonSimulation()
    mock_sim.terms = get_fully_enhanced_terms(experiment=evv_exp)

    mock_sim.exp = evv_exp

    mock_sim.setAxisChoiceAndTranslateTerms(axes_choice)

    mock_sim.spec_eval_setup = ws.main.spectrum_abstractions.SpecEvalSetup(ev_info=evi)

    # use simple model data
    mock_sim.system = datadict['system']
    mock_sim.props = datadict['props']
    mock_sim.vib_ana_setup = datadict['vib_ana_setup']
    mock_sim.vib_ana_setup.max_state_lvl = 3 # there is an issue for the underlying reason for this

    from wilson_suite.wilson_utils.termdict_from_symb_term import derived_terms_flat
    flat_dict = derived_terms_flat(mock_sim.terms, tolistonly=False)

    # TODO: how to update the terms for evaluation? how to make a selection of them after derivation?
    mock_sim.terms = [flat_dict['1_(1, 0)']]

    # FIXME: This printing appears to need update wrt. changes in wilson-derive, made issue
    #print('\n', flat_dict['1_(1, 0)'].to_latex())

    print(mock_sim.vib_ana_setup.max_state_lvl)
    
    mock_sim.evaluate()
    
    for f in mock_sim._workflow.artifacts.features:
        print(f.location, f.term_contributions[0].term_ids)
    
    region = mock_sim._workflow.artifacts.regions[0]
    feat1 = region.domain.full_features[0]
    feat_coeff = feat1.amplitude_coeff
    term_contributions = feat1.term_contributions

    print('\nterm_contributions[0].term_ids', term_contributions[0].term_ids, '\n')
    print('feat_coeff', feat_coeff)

    np.set_printoptions(linewidth=280, precision=1)
    for k,v in mock_sim._workflow.artifacts.grid_manager.full_grid.items():
        print(k,v)
    print('\n==========')
    
    from wilson_suite.wilson_utils.unit_convertor import convNu2Ene
    r_res = ws.intensities.amplitudes.evaluation_wf.evaluate_region(region, 
                                                            mock_sim._workflow.artifacts.vib_data, 
                                                            mock_sim._workflow.artifacts.vibdiff_cache, 
                                                            convNu2Ene(mock_sim.spec_eval_setup.ev_info.Gamma))
    ref_res = np.array([1/(-1j*convNu2Ene(1.))/(-1j*convNu2Ene(1.)) * feat_coeff])

    assert np.allclose(r_res, ref_res)
    assert np.allclose(ref_res, mock_sim.spec)

def test_evaluation_general_customdata_1mechterm():
    print()
    from ..unit.test_domains import get_data_evaluators_tests
    datadict = get_data_evaluators_tests()
    # 'system', 'vib_ana_setup', 'derived_terms', 'props', 
    # 'experiment', 'spec_eval_setup', 'domain_distance_thresholds'
    
    np.set_printoptions(linewidth=180, precision=3)

    from wilson_suite.wilson_main.workflow_abstractions import WilsonSimulation

    from wilson_suite.wilson_derive.derive import get_fully_enhanced_terms
    from ....fixtures import evv_experiment
    
    evv_exp = evv_experiment()
    terms = get_fully_enhanced_terms(experiment=evv_exp)
    axes_choice = evv_exp.valid_axis_combs[0].valid_axis_combs[3] # {'A': [(2,)], 'B': [(-1,), (2,)]}

    bounds_dict = {'B': (900., 900.), 'A': (1864., 1864.)}
    from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import SpectralWindow, Box
    spectral_window = SpectralWindow(box=Box(bounds_dict))

    evi = ws.main.spectrum_abstractions.EvaluationInfo(**{'spectral_window': spectral_window,
                                                          'Gamma': 1., 'Gamma_unit': 'cm-1',
                                                          'grid_resolution': {'A': 1, 'B': 1}})
    mock_sim = WilsonSimulation()
    mock_sim.terms = terms

    mock_sim.exp = evv_exp
    mock_sim.setAxisChoiceAndTranslateTerms(axes_choice)


    mock_sim.spec_eval_setup = ws.main.spectrum_abstractions.SpecEvalSetup(ev_info=evi)

    # use simple model data
    mock_sim.system = datadict['system']
    mock_sim.props = datadict['props']
    mock_sim.vib_ana_setup = datadict['vib_ana_setup']
    mock_sim.vib_ana_setup.max_state_lvl = 3 # there is an issue for the underlying reason for this

    print('\nmock_sim.is_ready', mock_sim.is_ready)

    from wilson_suite.wilson_utils.termdict_from_symb_term import derived_terms_flat
    flat_dict = derived_terms_flat(mock_sim.terms, tolistonly=False)

    # TODO: how to update the terms for evaluation? how to make a selection of them after derivation?

    # FIXME: Cannot find this term
    #mock_sim.terms = [flat_dict['8_(0, 1)']]
    #print('\n', flat_dict['8_(0, 1)'].to_latex())

    mock_sim.evaluate()

    for f in mock_sim._workflow.artifacts.features:
        print(f.location, f.term_contributions[0].term_ids)

    region = mock_sim._workflow.artifacts.regions[0]
    feat1 = region.domain.full_features[0]
    feat_coeff = feat1.amplitude_coeff
    term_contributions = feat1.term_contributions
    
    print('\nterm_contributions[0].term_ids', term_contributions[0].term_ids, '\n')
    print('feat_coeff', feat_coeff)

    np.set_printoptions(linewidth=280, precision=1)
    for k,v in mock_sim._workflow.artifacts.grid_manager.full_grid.items():
        print(k,v)
    print('\n==========')
    
    from wilson_suite.wilson_utils.unit_convertor import convNu2Ene
    r_res = ws.intensities.amplitudes.evaluation_wf.evaluate_region(region, 
                                                            mock_sim._workflow.artifacts.vib_data, 
                                                            mock_sim._workflow.artifacts.vibdiff_cache, 
                                                            convNu2Ene(mock_sim.spec_eval_setup.ev_info.Gamma))
    ref_res = np.array([1/(-1j*convNu2Ene(1.))/(-1j*convNu2Ene(1.)) * feat_coeff])
    assert np.allclose(r_res, ref_res)
    assert np.allclose(ref_res, mock_sim.spec)


def test_full_integration():
    print()
    from ....fixtures import evv_experiment
    from wilson_suite.wilson_utils.paths import SUITE_ROOT

    evv_exp = evv_experiment()
    terms = ws.derive.derive.get_fully_enhanced_terms(experiment=evv_exp)
    axes_choice = evv_exp.valid_axis_combs[0].valid_axis_combs[3] # {'A': [(2,)], 'B': [(-1,), (2,)]}

    calc_setup = ws.main.abstractions.DataOriginInfo(source_type='gaussian', 
                                                     lvl_theory='B3LYP', 
                                                     basis_set='cc-pVQZ', 
                                                     base_file_loc=SUITE_ROOT+'/../data_for_tests/g16_formaldehyde_B3LYPcc_pVQZ.out')

    sim = ws.main.workflow_abstractions.WilsonSimulation()
    sim.addExperiment(evv_exp)
    sim.addTerms(terms=terms) # terms

    mol_system = ws.main.abstractions.MolecularSystem(name='FORM', natoms=4)

    vib_ana = ws.main.abstractions.VibAnaSetup(system=mol_system, regime='GVPT2', vibana_own_analysis='none')
    
    sim.addSystem(mol_system)
    sim.addVibAnaSetup(vib_ana)
    sim.addPropEvalSetup(eval_uniform=calc_setup)
    
    sim.setPropsAndMaxStateLvl() # setting up self.props/sim.props
    sim.dressPropsWithSetup()

    sim.setAxisChoiceAndTranslateTerms(axes_choice)

    from wilson_suite.wilson_utils.wilson_data_obtainer import wilson_data_obtainer
    sim.getResults(obtainer=wilson_data_obtainer)

    from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import SpectralWindow, Box
    
    bounds_dict = {'A': (0., 5000.), 'B': (0., 5000.)}

    spectral_window = SpectralWindow(box=Box(bounds_dict))

    evi = ws.main.spectrum_abstractions.EvaluationInfo(**{'spectral_window': spectral_window,
                                                          'Gamma': 4.7, 'Gamma_unit': 'cm-1',
                                                          'grid_resolution': {'A': 7, 'B': 10}})
    
    eval_setup = ws.main.spectrum_abstractions.SpecEvalSetup(ev_info=evi)

    sim.addSpecEvalSetup(eval_setup)

    sim.vib_ana_setup.set_include_modes_list()

    print('simulation.exp.polarization_avg_vector', sim.exp.polarization_avg_vector)
    sim.evaluate()

    print(len(sim._workflow.artifacts.features))

    np.set_printoptions(linewidth=280, precision=1)

    import matplotlib.pyplot as plt

    Z = np.log(np.abs(sim.spec)**2)
    x = np.unique(sim.spec_eval_setup.grid['A'])
    y = np.unique(sim.spec_eval_setup.grid['B'])

    # if Z.shape == (len(y), len(x)) -> no transpose; if Z.shape == (len(x), len(y)) -> transpose
    # matplotlib expects [y, x] ordering for images
    toplot = Z.T

    plt.pcolormesh(x, y, toplot, shading="auto")
    plt.xlabel('A')
    plt.ylabel('B')
    plt.colorbar(label='log intensity')
    # plt.show()



def test_full_integration_EVV_axes():
    print()
    from ....fixtures import evv_experiment
    from wilson_suite.wilson_utils.paths import SUITE_ROOT

    evv_exp = evv_experiment()
    terms = ws.derive.derive.get_fully_enhanced_terms(experiment=evv_exp)

    from wilson_suite.wilson_utils.builder_functions import make_SpectralAxisSet
    axes_choice: ws.main.spectrum_abstractions.SpectralAxisSet = make_SpectralAxisSet({'A': [1], 'B': [-1,2]}) # {'A': [(1,)], 'B': [(-1,), (2,)]}

    calc_setup = ws.main.abstractions.DataOriginInfo(source_type='gaussian', 
                                                     lvl_theory='B3LYP', 
                                                     basis_set='cc-pVQZ', 
                                                     base_file_loc=SUITE_ROOT+'/../data_for_tests/g16_formaldehyde_B3LYPcc_pVQZ.out')

    sim = ws.main.workflow_abstractions.WilsonSimulation()
    sim.addExperiment(evv_exp)
    sim.addTerms(terms=terms) # terms

    mol_system = ws.main.abstractions.MolecularSystem(name='FORM', natoms=4)

    vib_ana = ws.main.abstractions.VibAnaSetup(system=mol_system, regime='GVPT2', vibana_own_analysis='none')
    
    sim.addSystem(mol_system)
    sim.addVibAnaSetup(vib_ana)
    sim.addPropEvalSetup(eval_uniform=calc_setup)
    
    sim.setPropsAndMaxStateLvl() # setting up self.props/sim.props
    sim.dressPropsWithSetup()

    sim.setAxisChoiceAndTranslateTerms(axes_choice)

    from wilson_suite.wilson_utils.wilson_data_obtainer import wilson_data_obtainer
    sim.getResults(obtainer=wilson_data_obtainer)

    from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import SpectralWindow, Box
    
    bounds_dict = {'A': (0., 5000.), 'B': (0., 5000.)}

    spectral_window = SpectralWindow(box=Box(bounds_dict))

    evi = ws.main.spectrum_abstractions.EvaluationInfo(**{'spectral_window': spectral_window,
                                                          'Gamma': 4.7, 'Gamma_unit': 'cm-1',
                                                          'grid_resolution': {'A': 7, 'B': 10}})
    
    eval_setup = ws.main.spectrum_abstractions.SpecEvalSetup(ev_info=evi)

    sim.addSpecEvalSetup(eval_setup)

    sim.vib_ana_setup.set_include_modes_list()

    print('simulation.exp.polarization_avg_vector', sim.exp.polarization_avg_vector)
    sim.evaluate()

    print(len(sim._workflow.artifacts.features))
    assert len(sim._workflow.artifacts.features) == 9

    np.set_printoptions(linewidth=280, precision=1)

    import matplotlib.pyplot as plt

    Z = np.log(np.abs(sim.spec)**2)
    x = np.unique(sim.spec_eval_setup.grid['A'])
    y = np.unique(sim.spec_eval_setup.grid['B'])

    # if Z.shape == (len(y), len(x)) -> no transpose; if Z.shape == (len(x), len(y)) -> transpose
    # matplotlib expects [y, x] ordering for images
    toplot = Z.T

    plt.pcolormesh(x, y, toplot, shading="auto")
    plt.xlabel('A')
    plt.ylabel('B')
    plt.colorbar(label='log intensity')
    # plt.show()


def test_full_integration__EVV_axes_with_apply_exp_magn_conditions():
    print()
    from ....fixtures import evv_experiment
    from wilson_suite.wilson_utils.paths import SUITE_ROOT

    evv_exp = evv_experiment()
    terms = ws.derive.derive.get_fully_enhanced_terms(experiment=evv_exp)

    from wilson_suite.wilson_utils.builder_functions import make_SpectralAxisSet
    axes_choice: ws.main.spectrum_abstractions.SpectralAxisSet = make_SpectralAxisSet({'A': [1], 'B': [-1,2]}) # {'A': [(1,)], 'B': [(-1,), (2,)]}

    calc_setup = ws.main.abstractions.DataOriginInfo(source_type='gaussian', 
                                                     lvl_theory='B3LYP', 
                                                     basis_set='cc-pVQZ', 
                                                     base_file_loc=SUITE_ROOT+'/../data_for_tests/g16_formaldehyde_B3LYPcc_pVQZ.out')

    sim = ws.main.workflow_abstractions.WilsonSimulation()
    sim.addExperiment(evv_exp)
    sim.addTerms(terms=terms) # terms

    mol_system = ws.main.abstractions.MolecularSystem(name='FORM', natoms=4)

    vib_ana = ws.main.abstractions.VibAnaSetup(system=mol_system, regime='GVPT2', vibana_own_analysis='none')
    
    sim.addSystem(mol_system)
    sim.addVibAnaSetup(vib_ana)
    sim.addPropEvalSetup(eval_uniform=calc_setup)
    
    sim.setPropsAndMaxStateLvl() # setting up self.props/sim.props
    sim.dressPropsWithSetup()

    sim.setAxisChoiceAndTranslateTerms(axes_choice)

    from wilson_suite.wilson_utils.wilson_data_obtainer import wilson_data_obtainer
    sim.getResults(obtainer=wilson_data_obtainer)

    from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import SpectralWindow, Box
    
    bounds_dict = {'A': (0., 5000.), 'B': (0., 5000.)}

    spectral_window = SpectralWindow(box=Box(bounds_dict))

    evi = ws.main.spectrum_abstractions.EvaluationInfo(**{'spectral_window': spectral_window,
                                                          'Gamma': 4.7, 'Gamma_unit': 'cm-1',
                                                          'grid_resolution': {'A': 7, 'B': 10}})
    
    eval_setup = ws.main.spectrum_abstractions.SpecEvalSetup(ev_info=evi)

    sim.addSpecEvalSetup(eval_setup)

    sim.vib_ana_setup.set_include_modes_list()

    print('simulation.exp.polarization_avg_vector', sim.exp.polarization_avg_vector)
    
    # this should filter features to draw
    sim.apply_exp_magn_conditions(where='eval')
    
    sim.evaluate()
    
    print(len(sim._workflow.artifacts.features))
    assert len(sim._workflow.artifacts.features) == 4

    np.set_printoptions(linewidth=280, precision=1)

    import matplotlib.pyplot as plt

    Z = np.log(np.abs(sim.spec)**2)
    x = np.unique(sim.spec_eval_setup.grid['A'])
    y = np.unique(sim.spec_eval_setup.grid['B'])

    # if Z.shape == (len(y), len(x)) -> no transpose; if Z.shape == (len(x), len(y)) -> transpose
    # matplotlib expects [y, x] ordering for images
    toplot = Z.T

    plt.pcolormesh(x, y, toplot, shading="auto")
    plt.xlabel('A')
    plt.ylabel('B')
    plt.colorbar(label='log intensity')
    # plt.show()


def test_full_integration_EVV_axes_dress_these_with_boxes_minimum_box_padding():
    print()
    from ....fixtures import evv_experiment
    from wilson_suite.wilson_utils.paths import SUITE_ROOT

    evv_exp = evv_experiment()
    terms = ws.derive.derive.get_fully_enhanced_terms(experiment=evv_exp)

    from wilson_suite.wilson_utils.builder_functions import make_SpectralAxisSet
    axes_choice: ws.main.spectrum_abstractions.SpectralAxisSet = make_SpectralAxisSet({'A': [1], 'B': [-1,2]}) # {'A': [(1,)], 'B': [(-1,), (2,)]}

    calc_setup = ws.main.abstractions.DataOriginInfo(source_type='gaussian', 
                                                     lvl_theory='B3LYP', 
                                                     basis_set='cc-pVQZ', 
                                                     base_file_loc=SUITE_ROOT+'/../data_for_tests/g16_formaldehyde_B3LYPcc_pVQZ.out')

    sim = ws.main.workflow_abstractions.WilsonSimulation()
    sim.addExperiment(evv_exp)
    sim.addTerms(terms=terms) # terms

    mol_system = ws.main.abstractions.MolecularSystem(name='FORM', natoms=4)

    vib_ana = ws.main.abstractions.VibAnaSetup(system=mol_system, regime='GVPT2', vibana_own_analysis='none')
    
    sim.addSystem(mol_system)
    sim.addVibAnaSetup(vib_ana)
    sim.addPropEvalSetup(eval_uniform=calc_setup)
    
    sim.setPropsAndMaxStateLvl() # setting up self.props/sim.props
    sim.dressPropsWithSetup()

    sim.setAxisChoiceAndTranslateTerms(axes_choice)

    from wilson_suite.wilson_utils.wilson_data_obtainer import wilson_data_obtainer
    sim.getResults(obtainer=wilson_data_obtainer)

    from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import SpectralWindow, Box
    
    bounds_dict = {'A': (0., 5000.), 'B': (0., 5000.)}

    spectral_window = SpectralWindow(box=Box(bounds_dict))

    evi = ws.main.spectrum_abstractions.EvaluationInfo(**{'spectral_window': spectral_window,
                                                          'Gamma': 4.7, 'Gamma_unit': 'cm-1',
                                                          'grid_resolution': {'A': 7, 'B': 10},
                                                          'minimum_box_padding': 10.}) # all features will have boxes now
    
    eval_setup = ws.main.spectrum_abstractions.SpecEvalSetup(ev_info=evi)

    sim.addSpecEvalSetup(eval_setup)

    sim.vib_ana_setup.set_include_modes_list()

    print('simulation.exp.polarization_avg_vector', sim.exp.polarization_avg_vector)
    sim.evaluate()

    print(len(sim._workflow.artifacts.features))
    assert len(sim._workflow.artifacts.features) == 60 # all of them now

    np.set_printoptions(linewidth=280, precision=1)

    import matplotlib.pyplot as plt

    Z = np.log(np.abs(sim.spec)**2)
    x = np.unique(sim.spec_eval_setup.grid['A'])
    y = np.unique(sim.spec_eval_setup.grid['B'])

    # if Z.shape == (len(y), len(x)) -> no transpose; if Z.shape == (len(x), len(y)) -> transpose
    # matplotlib expects [y, x] ordering for images
    toplot = Z.T

    plt.pcolormesh(x, y, toplot, shading="auto")
    plt.xlabel('A')
    plt.ylabel('B')
    plt.colorbar(label='log intensity')
    # plt.show()



def test_full_integration_other_axes_choice():
    print()
    from ....fixtures import evv_experiment
    from wilson_suite.wilson_utils.paths import SUITE_ROOT

    evv_exp = evv_experiment()
    terms = ws.derive.derive.get_fully_enhanced_terms(experiment=evv_exp)
    axes_choice = evv_exp.valid_axis_combs[0].valid_axis_combs[0] # {'A': [(2,)], 'B': [(-1,), (2,)]}

    calc_setup = ws.main.abstractions.DataOriginInfo(source_type='gaussian', 
                                                     lvl_theory='B3LYP', 
                                                     basis_set='cc-pVQZ', 
                                                     base_file_loc=SUITE_ROOT+'/../data_for_tests/g16_formaldehyde_B3LYPcc_pVQZ.out')

    sim = ws.main.workflow_abstractions.WilsonSimulation()
    sim.addExperiment(evv_exp)
    sim.addTerms(terms=terms) # terms

    mol_system = ws.main.abstractions.MolecularSystem(name='FORM', natoms=4)

    vib_ana = ws.main.abstractions.VibAnaSetup(system=mol_system, regime='GVPT2', vibana_own_analysis='none')
    
    sim.addSystem(mol_system)
    sim.addVibAnaSetup(vib_ana)
    sim.addPropEvalSetup(eval_uniform=calc_setup)
    
    sim.setPropsAndMaxStateLvl() # setting up self.props/sim.props
    sim.dressPropsWithSetup()

    sim.setAxisChoiceAndTranslateTerms(axes_choice)

    from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import SpectralWindow, Box
    
    bounds_dict = {'A': (0., 5000.), 'B': (0., 5000.)}

    spectral_window = SpectralWindow(box=Box(bounds_dict))

    evi = ws.main.spectrum_abstractions.EvaluationInfo(**{'spectral_window': spectral_window,
                                                          'Gamma': 4.7, 'Gamma_unit': 'cm-1',
                                                          'grid_resolution': {'A': 7, 'B': 10}})
    
    eval_setup = ws.main.spectrum_abstractions.SpecEvalSetup(ev_info=evi)

    sim.addSpecEvalSetup(eval_setup)

    from wilson_suite.wilson_utils.wilson_data_obtainer import wilson_data_obtainer
    sim.getResults(obtainer=wilson_data_obtainer)
    sim.vib_ana_setup.set_include_modes_list()

    import pytest
    with pytest.raises(ValueError) as error:
        sim.evaluate()
    assert str(error.value) == "Failed at 'place_in_specwindow': This SpectralWindow does not contain any features. Change the bounds of the window or use different terms. EvaluationWorkflow instanse was saved to `eval_wf.pkl`."

    from wilson_suite.wilson_utils.serialization import unpickle_smth_from
    eval_wf: ws.intensities.amplitudes.evaluation_wf.EvaluationWorkflow = unpickle_smth_from('eval_wf.pkl')

    assert eval_wf.inputs.spec_eval_setup == sim.spec_eval_setup
    # assert eval_wf.inputs.props == sim.props # doesn't work because cff has extra data now. TODO fix later
    assert eval_wf.inputs.vib_ana_setup == sim.vib_ana_setup
    assert eval_wf.artifacts.spec_window == sim.spec_eval_setup.ev_info.spectral_window
    assert eval_wf.artifacts.regions is None
    assert len(eval_wf.artifacts.terms) == 14
    assert eval_wf.artifacts.grid_manager is None

    import os
    os.remove('eval_wf.pkl')


def test_integration_evv_experiment_until_after_evaluation():

    from ....fixtures import evv_experiment
    from wilson_suite.wilson_utils.paths import SUITE_ROOT
    from wilson_suite.wilson_experiment.indep_vars_and_axes import SpectralAxisSet, IndependentVariableSet, \
        SignedPulseTuple, SpectralAxis

    evv_exp = evv_experiment()
    terms = ws.derive.derive.get_fully_enhanced_terms(experiment=evv_exp)
    #axes_choice = evv_exp.valid_axis_combs[0].valid_axis_combs[1]  # {'A': [(-1,)], 'B': [(2,)]}
    axes_choice = evv_exp.valid_axis_combs[0].valid_axis_combs[0] # {'A': [(-1,)], 'B': [(-1,), (2,)]}
    # axis_choice = SpectralAxisSet(
    #axes = (SpectralAxis(label='A', var_set=IndependentVariableSet(var_set=(SignedPulseTuple(pulse_refs=(1,)),))),
    #        SpectralAxis(label='B', var_set=IndependentVariableSet(
    #            var_set=(SignedPulseTuple(pulse_refs=(-1,)),
    #                     SignedPulseTuple(pulse_refs=(2,)))))))  # {'A': [(1,)], 'B': [(-1,), (2,)]}



    calc_setup = ws.main.abstractions.DataOriginInfo(source_type='gaussian',
                                                     lvl_theory='B3LYP',
                                                     basis_set='cc-pVQZ',
                                                     base_file_loc=SUITE_ROOT+'/../data_for_tests/g16_formaldehyde_B3LYPcc_pVQZ.out')

    sim = ws.main.workflow_abstractions.WilsonSimulation()
    sim.addExperiment(evv_exp)
    sim.addTerms(terms=terms)  # terms

    mol_system = ws.main.abstractions.MolecularSystem(name='formaldehyde', natoms=4)

    vib_ana = ws.main.abstractions.VibAnaSetup(system=mol_system, regime='GVPT2', vibana_own_analysis='none')

    sim.addSystem(mol_system)
    sim.addVibAnaSetup(vib_ana)
    sim.addPropEvalSetup(eval_uniform=calc_setup)

    sim.setPropsAndMaxStateLvl()  # setting up self.props/sim.props
    sim.dressPropsWithSetup()

    sim.setAxisChoiceAndTranslateTerms(axes_choice)

    from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import SpectralWindow, Box

    # These windows capture many of the same features for the respective axis set choices
    #bounds_dict = {'A': (-3000, -1200.), 'B': (1200., 6000.)} # {'A': [(-1,)], 'B': [(2,)]}
    bounds_dict = {'A': (-3000, -1200.), 'B': (500., 3000.)} # {'A': [(-1,)], 'B': [(-1,), (2,)]}

    # These bounds raise a "no features in window" error but I thought they would correspond to what I used for
    # {'A': [(-1,)], 'B': [(-1,), (2,)]} right above here
    # bounds_dict = {'A': (1200, 3000.), 'B': (500., 3000.)} # {'A': [(1,)], 'B': [(-1,), (2,)]}

    spectral_window = SpectralWindow(box=Box(bounds_dict))

    dynrange_log10 = 9 # 3 = dynamic range 1000

    evi = ws.main.spectrum_abstractions.EvaluationInfo(**{'spectral_window': spectral_window,
                                                          'Gamma': 4.7, 'Gamma_unit': 'cm-1',
                                                          'grid_resolution': {'A': 360, 'B': 600},
                                                          'dynamic_range': 10**dynrange_log10,
                                                          'box_range_safety_margin': 0.1,
                                                          'scale_wrt_max_intensity': True,
                                                          'minimum_box_padding': 30.0
                                                          }
                                                       )


    eval_setup = ws.main.spectrum_abstractions.SpecEvalSetup(ev_info=evi)

    sim.addSpecEvalSetup(eval_setup)

    from wilson_suite.wilson_utils.wilson_data_obtainer import wilson_data_obtainer
    sim.getResults(obtainer=wilson_data_obtainer)
    sim.vib_ana_setup.set_include_modes_list()

    sim.evaluate()

    import matplotlib.pyplot as plt


    Z = np.log(np.abs(sim.spec)**2)

    zmax = np.amax(Z)

    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            if Z[i, j] < (zmax - dynrange_log10):
                Z[i, j] = (zmax - dynrange_log10)


    x = np.unique(sim.spec_eval_setup.grid['A'])
    y = np.unique(sim.spec_eval_setup.grid['B'])

    # if Z.shape == (len(y), len(x)) -> no transpose; if Z.shape == (len(x), len(y)) -> transpose
    # matplotlib expects [y, x] ordering for images
    toplot = Z.T

    plt.pcolormesh(x, y, toplot, vmax=np.amax(Z), vmin=np.amax(Z)-dynrange_log10, shading="auto")
    plt.xlabel('A')
    plt.ylabel('B')
    plt.colorbar(label='log intensity')
    plt.show()

def test_full_integration_H2O_molecule():
    print()
    from ....fixtures import evv_experiment
    from wilson_suite.wilson_utils.paths import SUITE_ROOT

    evv_exp = evv_experiment()
    terms = ws.derive.derive.get_fully_enhanced_terms(experiment=evv_exp)
    axes_choice = evv_exp.valid_axis_combs[0].valid_axis_combs[3] # {'A': [(2,)], 'B': [(-1,), (2,)]}

    calc_setup = ws.main.abstractions.DataOriginInfo(source_type='gaussian',
                                                     lvl_theory='HF', 
                                                     basis_set='STO-3G', 
                                                     base_file_loc=SUITE_ROOT+'/../data_for_tests/g16_h2o_HF_STO3G.out')

    sim = ws.main.workflow_abstractions.WilsonSimulation()
    sim.addExperiment(evv_exp)
    sim.addTerms(terms=terms) # terms

    mol_system = ws.main.abstractions.MolecularSystem(name='h2o', natoms=3)

    vib_ana = ws.main.abstractions.VibAnaSetup(system=mol_system, regime='GVPT2', vibana_own_analysis='none')
    
    sim.addSystem(mol_system)
    sim.addVibAnaSetup(vib_ana)
    sim.addPropEvalSetup(eval_uniform=calc_setup)
    
    sim.setPropsAndMaxStateLvl() # setting up self.props/sim.props
    sim.dressPropsWithSetup()

    sim.setAxisChoiceAndTranslateTerms(axes_choice)

    from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import SpectralWindow, Box
    
    bounds_dict = {'A': (0., 5000.), 'B': (0., 5000.)}

    spectral_window = SpectralWindow(box=Box(bounds_dict))

    evi = ws.main.spectrum_abstractions.EvaluationInfo(**{'spectral_window': spectral_window,
                                                          'Gamma': 4.7, 'Gamma_unit': 'cm-1',
                                                          'grid_resolution': {'A': 10, 'B': 10}})
    
    eval_setup = ws.main.spectrum_abstractions.SpecEvalSetup(ev_info=evi)

    sim.addSpecEvalSetup(eval_setup)

    from wilson_suite.wilson_utils.wilson_data_obtainer import wilson_data_obtainer
    sim.getResults(obtainer=wilson_data_obtainer)
    sim.vib_ana_setup.set_include_modes_list()

    sim.evaluate()

    np.set_printoptions(linewidth=280, precision=1)

    import matplotlib.pyplot as plt

    Z = np.log(np.abs(sim.spec)**2)
    x = np.unique(sim.spec_eval_setup.grid['A'])
    y = np.unique(sim.spec_eval_setup.grid['B'])

    # if Z.shape == (len(y), len(x)) -> no transpose; if Z.shape == (len(x), len(y)) -> transpose
    # matplotlib expects [y, x] ordering for images
    toplot = Z.T

    plt.pcolormesh(x, y, toplot, shading="auto")
    plt.xlabel('A')
    plt.ylabel('B')
    plt.colorbar(label='log intensity')
    #plt.show()



def test_full_integration_EVV_axes_getResults_extra():
    print()
    from ....fixtures import evv_experiment
    from wilson_suite.wilson_utils.paths import SUITE_ROOT

    evv_exp = evv_experiment()
    terms = ws.derive.derive.get_fully_enhanced_terms(experiment=evv_exp)

    from wilson_suite.wilson_utils.builder_functions import make_SpectralAxisSet
    axes_choice: ws.main.spectrum_abstractions.SpectralAxisSet = make_SpectralAxisSet({'A': [1], 'B': [-1,2]}) # {'A': [(1,)], 'B': [(-1,), (2,)]}

    calc_setup = ws.main.abstractions.DataOriginInfo(source_type='gaussian', 
                                                     lvl_theory='B3LYP', 
                                                     basis_set='cc-pVQZ', 
                                                     base_file_loc=SUITE_ROOT+'/../data_for_tests/g16_formaldehyde_B3LYPcc_pVQZ.out')

    sim = ws.main.workflow_abstractions.WilsonSimulation()
    sim.addExperiment(evv_exp)
    sim.addTerms(terms=terms) # terms

    mol_system = ws.main.abstractions.MolecularSystem(name='FORM', natoms=4)

    vib_ana = ws.main.abstractions.VibAnaSetup(system=mol_system, regime='GVPT2', vibana_own_analysis='none')
    
    sim.addSystem(mol_system)
    sim.addVibAnaSetup(vib_ana)
    sim.addPropEvalSetup(eval_uniform=calc_setup)
    
    sim.setPropsAndMaxStateLvl() # setting up self.props/sim.props
    sim.dressPropsWithSetup()

    sim.setAxisChoiceAndTranslateTerms(axes_choice)

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>
    from wilson_suite.wilson_utils.wilson_data_obtainer import wilson_data_obtainer
    sim.getResults(obtainer=wilson_data_obtainer, get_geometry=True, get_displacements=True)
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>

    from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import SpectralWindow, Box
    
    bounds_dict = {'A': (0., 5000.), 'B': (0., 5000.)}

    spectral_window = SpectralWindow(box=Box(bounds_dict))

    evi = ws.main.spectrum_abstractions.EvaluationInfo(**{'spectral_window': spectral_window,
                                                          'Gamma': 4.7, 'Gamma_unit': 'cm-1',
                                                          'grid_resolution': {'A': 7, 'B': 10}})
    
    eval_setup = ws.main.spectrum_abstractions.SpecEvalSetup(ev_info=evi)

    sim.addSpecEvalSetup(eval_setup)

    sim.vib_ana_setup.set_include_modes_list()

    print('simulation.exp.polarization_avg_vector', sim.exp.polarization_avg_vector)
    sim.evaluate()

    print(len(sim._workflow.artifacts.features))
    assert len(sim._workflow.artifacts.features) == 9

    np.set_printoptions(linewidth=280, precision=1)

    import matplotlib.pyplot as plt

    Z = np.log(np.abs(sim.spec)**2)
    x = np.unique(sim.spec_eval_setup.grid['A'])
    y = np.unique(sim.spec_eval_setup.grid['B'])

    # if Z.shape == (len(y), len(x)) -> no transpose; if Z.shape == (len(x), len(y)) -> transpose
    # matplotlib expects [y, x] ordering for images
    toplot = Z.T

    plt.pcolormesh(x, y, toplot, shading="auto")
    plt.xlabel('A')
    plt.ylabel('B')
    plt.colorbar(label='log intensity')
    # plt.show()


'''
def test_debugging():
    print()
    from wilson_suite.wilson_utils.serialization import unpickle_smth_from
    wf = unpickle_smth_from('/home/vlev/monorepo/eval_wf.pkl')
    print(type(wf))
    from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import SpectralFeature
    SpectralFeature.print_list_features(wf.artifacts.features)
'''