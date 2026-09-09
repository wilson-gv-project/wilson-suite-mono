import wilson_suite as ws
from wilson_suite.wilson_main.abstractions import DataOriginInfo
import numpy as np
from wilson_suite.wilson_utils.wilson_data_obtainer import wilson_data_obtainer

def test_saving_obtained_data():
    """
    complete_info_keys = ['cff', 'anharmonic_states', 'nc_sqrt_eigval', 'dipgrad', 'B', 'polgrad', 'coriolis', 'polhess', 'qff', 'diphess']
    compl_data = complete_info_keys + ['harmonic_states']
    """
    from ...fixtures import evv_experiment
    from wilson_suite.wilson_utils.paths import SUITE_ROOT

    evv_exp = evv_experiment()
    terms = ws.derive.derive.get_fully_enhanced_terms(experiment=evv_exp)

    sim = ws.main.workflow_abstractions.WilsonSimulation()
    sim.addExperiment(evv_exp)
    sim.addTerms(terms=terms) # terms

    mol_system = ws.main.abstractions.MolecularSystem(name='h2o', natoms=3)
    vib_ana = ws.main.abstractions.VibAnaSetup(system=mol_system, regime='GVPT2', vibana_own_analysis='none')
    
    sim.addSystem(mol_system)
    sim.addVibAnaSetup(vib_ana)

    calc_setup = ws.main.abstractions.DataOriginInfo(source_type='gaussian',
                                                     lvl_theory='HF', 
                                                     basis_set='STO-3G', 
                                                     base_file_loc=SUITE_ROOT+'/../data_for_tests/g16_h2o_HF_STO3G.out')
    calc_setup_blank = ws.main.abstractions.DataOriginInfo()
    sim.addPropEvalSetup(eval_uniform=calc_setup_blank)
    
    sim.setPropsAndMaxStateLvl() # setting up self.props/sim.props
    sim.dressPropsWithSetup()

    data_orig_g16 = DataOriginInfo(source_type='gaussian',
                                   base_file_loc='/home/vlev/monorepo/src/../data_for_tests/g16_h2o_HF_STO3G.out')
    rq_none_keys = ['dipgrad', 'polhess', 'polgrad', 'diphess', 'cff', 'nc_sqrt_eigval', 'anharmonic_states']
    request_dict = dict.fromkeys(rq_none_keys, data_orig_g16)

    # pickling data dict - testing save_obtained_data
    from wilson_suite.wilson_utils.paths import SUITE_ROOT
    filepath = SUITE_ROOT+'/wilson_suite/wilson_main/tests/test_compl_data.pkl'
    print(filepath, type(filepath))
    ws.utils.save_obtained_data(request_dict, format='pkl', filename=filepath)

    unpkl_compl_data = ws.utils.serialization.unpickle_smth_from(filepath)
    assert_equal(unpkl_compl_data, request_dict)

    # ---- workflow getResults
    sim.addPropEvalSetup(eval_uniform=calc_setup)
    sim.setPropsAndMaxStateLvl()
    sim.dressPropsWithSetup()

    # saving data
    filepath_wf = SUITE_ROOT+'/wilson_suite/wilson_main/tests/myfile.pkl'
    sim.getResults(obtainer=wilson_data_obtainer, save_to_filename=filepath_wf)

    # for requestData to have reset values for props and residual_vib_info
    sim.dressPropsWithSetup()
    
    compl_data_wf = wilson_data_obtainer(sim.requestData())
    unpkl_compl_data_wf = ws.utils.serialization.unpickle_smth_from(filepath_wf)
    assert_equal(unpkl_compl_data_wf, compl_data_wf)

def assert_equal(a, b):
    """
    assert equality of complete obtained data arrays

    some values are np.ndarrays, others are dicts with float values
    """

    assert type(a) is type(b)

    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for k in a:
            assert_equal(a[k], b[k])

    elif isinstance(a, np.ndarray):
        assert np.array_equal(a, b)

    else:
        # float comparison
        assert a == b


def test_save_wilsonsim():
    from ...fixtures import evv_experiment
    from wilson_suite.wilson_utils.paths import SUITE_ROOT

    evv_exp = evv_experiment()
    terms = ws.derive.derive.get_fully_enhanced_terms(experiment=evv_exp)

    sim = ws.main.workflow_abstractions.WilsonSimulation()
    sim.addExperiment(evv_exp)
    sim.addTerms(terms=terms)

    mol_system = ws.main.abstractions.MolecularSystem(name='h2o', natoms=3)
    vib_ana = ws.main.abstractions.VibAnaSetup(system=mol_system, regime='GVPT2', vibana_own_analysis='none')
    
    sim.addSystem(mol_system)
    sim.addVibAnaSetup(vib_ana)

    calc_setup = ws.main.abstractions.DataOriginInfo(source_type='gaussian',
                                                     lvl_theory='HF', 
                                                     basis_set='STO-3G', 
                                                     base_file_loc=SUITE_ROOT+'/../data_for_tests/g16_h2o_HF_STO3G.out')
    sim.addPropEvalSetup(eval_uniform=calc_setup)
    
    sim.setPropsAndMaxStateLvl() # setting up self.props/sim.props
    sim.dressPropsWithSetup()

    from wilson_suite.wilson_utils.paths import WORKFLOW_BASE_DIR
    print(WORKFLOW_BASE_DIR, type(WORKFLOW_BASE_DIR))

    # initialize WORKFLOW_BASE_DIR, if not there; initialize this sim workflow run directory
    sim.make_proj_dir()
    # save this sim to pickle in the sim workflow run directory
    sim.save_to_pkl(filename='sim0.pkl')

    sim.getResults(obtainer=wilson_data_obtainer, save_to_filename='data_file.pkl')

    vib_ana.set_include_modes_list()
    
    from wilson_suite.wilson_utils.builder_functions import make_SpectralAxisSet
    axes_choice = make_SpectralAxisSet({'A': [2], 'B': [-1, 2]})
    sim.setAxisChoiceAndTranslateTerms(axes_choice)

    from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import SpectralWindow, Box
    bounds_dict = {'A': (1000., 3100.), 'B': (-100., 2500.)}
    spectral_window = SpectralWindow(box=Box(bounds_dict))

    evi = ws.main.spectrum_abstractions.EvaluationInfo(**{'spectral_window': spectral_window,
                                                          'Gamma': 24.7, 'Gamma_unit': 'cm-1',
                                                          'dynamic_range': 1000,
                                                          'grid_resolution': {'A': 70, 'B': 100}})
    eval_setup = ws.main.spectrum_abstractions.SpecEvalSetup(ev_info=evi)

    sim.addSpecEvalSetup(eval_setup)
    sim.save_to_pkl(filename='sim1.pkl')

    # save EvaluationInputs for evaluation workflow
    sim.evaluate(save_evalinputs_pkl='EvalInputs.pkl')

    sim.save_to_pkl(filename='configs.pkl', configs_only=True)

