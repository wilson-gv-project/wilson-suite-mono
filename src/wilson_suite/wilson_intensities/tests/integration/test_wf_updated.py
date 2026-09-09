import wilson_suite as ws
import numpy as np


# def test_full_integration_EVV_axes():
#     print()
#     from ....fixtures import evv_experiment
#     from wilson_suite.wilson_utils.paths import SUITE_ROOT

#     EVV_EXPERIMENT = evv_experiment()
#     DERIVED_EVV_TERMS = ws.derive.derive.get_fully_enhanced_terms(experiment=EVV_EXPERIMENT)

#     from wilson_suite.wilson_utils.builder_functions import make_SpectralAxisSet
#     axes_choice: ws.main.spectrum_abstractions.SpectralAxisSet = make_SpectralAxisSet({'A': [1], 'B': [-1,2]}) # {'A': [(1,)], 'B': [(-1,), (2,)]}

#     calc_setup = ws.main.abstractions.DataOriginInfo(source_type='gaussian', 
#                                                      lvl_theory='B3LYP', 
#                                                      basis_set='cc-pVQZ', 
#                                                      base_file_loc=SUITE_ROOT+'/../data_for_tests/g16_formaldehyde_B3LYPcc_pVQZ.out')

#     mol_system = ws.main.abstractions.MolecularSystem(name='FORM', natoms=4)

#     vib_ana = ws.main.abstractions.VibAnaSetup(system=mol_system, regime='GVPT2', vibana_own_analysis='none')
#     vib_ana.set_include_modes_list()


#     from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import SpectralWindow, Box
#     bounds_dict = {'A': (0., 5000.), 'B': (0., 5000.)}
#     spectral_window = SpectralWindow(box=Box(bounds_dict))

#     evi = ws.main.spectrum_abstractions.EvaluationInfo(**{'spectral_window': spectral_window,
#                                                           'Gamma': 4.7, 'Gamma_unit': 'cm-1',
#                                                           'grid_resolution': {'A': 7, 'B': 10}})
    
#     eval_setup = ws.main.spectrum_abstractions.SpecEvaxlSetup(ev_info=evi)

#     from wilson_suite.wilson_main.workflow_abstractions_updated import SimulationBuilder
#     sim = (SimulationBuilder(EVV_EXPERIMENT)
#                         .with_terms(DERIVED_EVV_TERMS)
#                         .with_system(mol_system)
#                         .with_vibana_setup(vib_ana)
#                         .with_calc_setup(calc_setup)
#                         .with_eval_setup(eval_setup)
#                         .build())

#     # only now does execution begin
#     data = sim.request_data()
#     results = sim.evaluate(data)
#     results.render()

#     sim.addSpecEvalSetup(eval_setup)
#     sim.vib_ana_setup.set_include_modes_list()


#     print(len(sim._workflow.artifacts.features))
#     assert len(sim._workflow.artifacts.features) == 9

#     np.set_printoptions(linewidth=280, precision=1)

#     import matplotlib.pyplot as plt

#     Z = np.log(np.abs(sim.spec)**2)
#     x = np.unique(sim.spec_eval_setup.grid['A'])
#     y = np.unique(sim.spec_eval_setup.grid['B'])

#     # if Z.shape == (len(y), len(x)) -> no transpose; if Z.shape == (len(x), len(y)) -> transpose
#     # matplotlib expects [y, x] ordering for images
#     toplot = Z.T

#     plt.pcolormesh(x, y, toplot, shading="auto")
#     plt.xlabel('A')
#     plt.ylabel('B')
#     plt.colorbar(label='log intensity')
#     # plt.show()


# def test_full_integration__EVV_axes_with_apply_exp_magn_conditions():
#     print()
#     from wilson_suite.wilson_utils.paths import SUITE_ROOT
#     from ....fixtures import evv_experiment

#     evv_exp = evv_experiment()
#     terms = ws.derive.derive.get_fully_enhanced_terms(experiment=evv_exp)

#     from wilson_suite.wilson_utils.builder_functions import make_SpectralAxisSet
#     axes_choice: ws.main.spectrum_abstractions.SpectralAxisSet = make_SpectralAxisSet({'A': [1], 'B': [-1,2]}) # {'A': [(1,)], 'B': [(-1,), (2,)]}

#     calc_setup = ws.main.abstractions.DataOriginInfo(source_type='gaussian', 
#                                                      lvl_theory='B3LYP', 
#                                                      basis_set='cc-pVQZ', 
#                                                      base_file_loc=SUITE_ROOT+'/../data_for_tests/g16_formaldehyde_B3LYPcc_pVQZ.out')

#     sim = ws.main.workflow_abstractions.WilsonSimulation()
#     sim.addExperiment(evv_exp)
#     sim.addTerms(terms=terms) # terms

#     mol_system = ws.main.abstractions.MolecularSystem(name='FORM', natoms=4)

#     vib_ana = ws.main.abstractions.VibAnaSetup(system=mol_system, regime='GVPT2', vibana_own_analysis='none')
    
#     sim.addSystem(mol_system)
#     sim.addVibAnaSetup(vib_ana)
#     sim.addPropEvalSetup(eval_uniform=calc_setup)
    
#     sim.setPropsAndMaxStateLvl() # setting up self.props/sim.props
#     sim.dressPropsWithSetup()

#     sim.setAxisChoiceAndTranslateTerms(axes_choice)

#     from wilson_suite.wilson_utils.wilson_data_obtainer import wilson_data_obtainer
#     sim.getResults(obtainer=wilson_data_obtainer)

#     from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import SpectralWindow, Box
    
#     bounds_dict = {'A': (0., 5000.), 'B': (0., 5000.)}

#     spectral_window = SpectralWindow(box=Box(bounds_dict))

#     evi = ws.main.spectrum_abstractions.EvaluationInfo(**{'spectral_window': spectral_window,
#                                                           'Gamma': 4.7, 'Gamma_unit': 'cm-1',
#                                                           'grid_resolution': {'A': 7, 'B': 10}})
    
#     eval_setup = ws.main.spectrum_abstractions.SpecEvalSetup(ev_info=evi)

#     sim.addSpecEvalSetup(eval_setup)

#     sim.vib_ana_setup.set_include_modes_list()

#     print('simulation.exp.polarization_avg_vector', sim.exp.polarization_avg_vector)
    
#     # this should filter features to draw
#     sim.apply_exp_magn_conditions(where='eval')
    
#     sim.evaluate()
    
#     print(len(sim._workflow.artifacts.features))
#     assert len(sim._workflow.artifacts.features) == 4

#     np.set_printoptions(linewidth=280, precision=1)

#     import matplotlib.pyplot as plt

#     Z = np.log(np.abs(sim.spec)**2)
#     x = np.unique(sim.spec_eval_setup.grid['A'])
#     y = np.unique(sim.spec_eval_setup.grid['B'])

#     # if Z.shape == (len(y), len(x)) -> no transpose; if Z.shape == (len(x), len(y)) -> transpose
#     # matplotlib expects [y, x] ordering for images
#     toplot = Z.T

#     plt.pcolormesh(x, y, toplot, shading="auto")
#     plt.xlabel('A')
#     plt.ylabel('B')
#     plt.colorbar(label='log intensity')
#     # plt.show()


def test_new():
    from wilson_suite.wilson_utils.paths import SUITE_ROOT
    from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import SpectralWindow, Box
    from wilson_suite.wilson_main.workflow_abstractions_updated import SimulationBuilder, SimulationRun
    from ....fixtures import evv_experiment
    from wilson_suite.wilson_utils.builder_functions import make_SpectralAxisSet

    evv_exp = evv_experiment()
    mol_system = ws.main.abstractions.MolecularSystem(name='FORM', natoms=4)
    calc_setup = ws.main.abstractions.DataOriginInfo(source_type='gaussian', 
                                                     lvl_theory='B3LYP', 
                                                     basis_set='cc-pVQZ', 
                                                     base_file_loc=SUITE_ROOT+'/../data_for_tests/g16_formaldehyde_B3LYPcc_pVQZ.out')

    vib_ana = ws.main.abstractions.VibAnaSetup(system=mol_system, regime='GVPT2', vibana_own_analysis='none')
    
    axes_choice: ws.main.spectrum_abstractions.SpectralAxisSet = make_SpectralAxisSet({'A': [1], 'B': [-1,2]}) # {'A': [(1,)], 'B': [(-1,), (2,)]}
    bounds_dict = {'A': (0., 5000.), 'B': (0., 5000.)}
    spectral_window = SpectralWindow(box=Box(bounds_dict))
    evi = ws.main.spectrum_abstractions.EvaluationInfo(**{'spectral_window': spectral_window,
                                                          'spectral_axes': axes_choice,
                                                          'Gamma': 4.7, 'Gamma_unit': 'cm-1',
                                                          'grid_resolution': {'A': 7, 'B': 10},
                                                          'minimum_box_padding': 10.,
                                                          'apply_magn_conditions': 'eval'})
    
    style_config = ws.main.spectrum_abstractions.PlotConfig(tick_step=50., 
                                                            colormap='magma_r',
                                                            show_top_ticks=False,
                                                            show_right_ticks=False
                                                            )
    rnd = ws.main.spectrum_abstractions.RenderingInfo(filename='fig_file.svg', 
                                                      reference_max=None,
                                                      style_config=style_config,
                                                      axes_labels={"x": "w1", "y": "w2-w1"})
    
    eval_setup = ws.main.spectrum_abstractions.SpecEvalSetup(ev_info=evi, rnd_info=rnd)

    builder = SimulationBuilder(experiment=evv_exp, system=mol_system)

    # validate these structures
    builder.set_eval_uniform(calc_setup)
    builder.set_spec_eval(eval_setup)
    builder.set_vib_ana(vib_ana)

    builder.dress_prop_order()
    # contains empty setup, without actual data for evaluation, with terms
    sealed_input = builder.seal()

    print(sealed_input.spec_eval.rnd_info)
    runner = SimulationRun(sealed=sealed_input)
    res = runner.execute()

    assert runner.is_data_filled
    assert runner.has_result

    assert res.artifacts.vibdiff_cache is not None
    assert res.artifacts.vibdiff_cache._cache != {}
    assert res.artifacts.motif_locs is not None
    assert res.artifacts.need_precalc is not None
    assert res.artifacts.precalculated is not None
    assert res.artifacts.coefficients is not None
    assert res.artifacts.features is not None
    assert res.artifacts.zero_feats is not None
    assert res.artifacts.grid_manager is not None
    assert res.artifacts.regions is not None
    assert res.artifacts.regions_results is not None

    runner.render()
    """
    """

def test_new_dryrun():
    from wilson_suite.wilson_main.workflow_abstractions_updated import SimulationBuilder
    from ....fixtures import evv_experiment
    from wilson_suite.wilson_utils.builder_functions import make_SpectralAxisSet

    evv_exp = evv_experiment()

    vib_ana = ws.main.abstractions.VibAnaSetup(regime='GVPT2', vibana_own_analysis='none')
    
    builder = SimulationBuilder(experiment=evv_exp)

    # validate these structures
    builder.set_vib_ana(vib_ana)

    axes_choice: ws.main.spectrum_abstractions.SpectralAxisSet = make_SpectralAxisSet({'A': [1], 'B': [-1,2]}) # {'A': [(1,)], 'B': [(-1,), (2,)]}
    builder.set_axis_choice(axes_choice)

    # contains empty setup, without actual data for evaluation, with terms
    sealed_input = builder.seal_dry_run()

    assert sealed_input.terms_in_axes.axes_choice == axes_choice
    assert sealed_input.terms_in_axes.magn_conditions == (('B',),)
    
