from wilson_suite.wilson_analysis.render.matplotlib_renderer import MatplotlibRendererGV
import wilson_suite as ws

from wilson_suite.wilson_utils.paths import SUITE_ROOT
from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import SpectralWindow, Box
from wilson_suite.wilson_main.workflow_abstractions_updated import SimulationBuilder, SimulationRun
from wilson_suite.fixtures import evv_experiment
from wilson_suite.wilson_utils.builder_functions import make_SpectralAxisSet

def get_eval_result_and_input():
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
                                                            'Gamma': 10., 'Gamma_unit': 'cm-1',
                                                            'grid_resolution': {'A': 700, 'B': 1000},
                                                            'minimum_box_padding': 250.,
                                                            'apply_magn_conditions': 'eval'})

    style_config = ws.main.spectrum_abstractions.PlotConfig(tick_step=50., 
                                                            cb_tick_label_fontsize=24,
                                                            colormap='magma_r',
                                                            show_top_ticks=False,
                                                            show_right_ticks=False
                                                            )
    rnd = ws.main.spectrum_abstractions.RenderingInfo(filename='fig_file.svg', 
                                                        reference_max=None,
                                                        style_config=style_config,
                                                        dynamic_range=1000,
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
    return res, sealed_input

def test_rndGV():
    res, sealed_input = get_eval_result_and_input()

    print(res.grid)
    renderer = MatplotlibRendererGV(eval_result=res, 
                                    setup_inputs=sealed_input)
    ax = renderer.contour()
    renderer.save(ax.figure, path='fig1.svg')