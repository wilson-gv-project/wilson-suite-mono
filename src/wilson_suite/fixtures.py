from . import wilson_experiment as ws_experiment
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from wilson_suite.wilson_derive.response_terms import VibPerturbedTerm


import logging
# wilson. - for hierarchy of loggers
logger = logging.getLogger("wilson.")


def evv_experiment() -> ws_experiment.experiment_abstractions.VibExperiment:
    """
    Returns VibExperiment instance for EVV experiment
    """
    import wilson_suite.wilson_experiment.experiment_abstractions as wexp

    pulse_ir_1 = wexp.make_impulsive_gaussian_pulse(tc=50.0, cf=0.0, cf_uv=0.0,
                                                           maxstr=1.0e-5, wv=(0.0, 0.0, 1.0), pol=(1.0, 0.0, 0.0), id=1)

    pulse_ir_2 = wexp.make_impulsive_gaussian_pulse(tc=100.0, cf=0.0, cf_uv=0.0,
                                                           maxstr=1.0e-5, wv=(0.0, 0.0, 1.0), pol=(1.0, 0.0, 0.0), id=2)

    pulse_uvvis_1 = wexp.make_impulsive_gaussian_pulse(tc=120.0, cf=0.0, cf_uv=0.072,
                                                           maxstr=1.0e-5, wv=(0.0, 0.0, 1.0), pol=(1.0, 0.0, 0.0), id=3)

    pulses = (pulse_ir_1, pulse_ir_2, pulse_uvvis_1)

    field_a = wexp.ElectricField(pulses)

    detector_a = wexp.SpecDetector(detection_method='freq',
                                                         detector_location=(0.0, 0.0, 1.0),
                                                         detection_polarization=(1.0, 0.0, 0.0),
                                                         detection_range=[0.003 + 0.0001 * i for i in range(101)],
                                                         wv_filter=[{1: -1, 2: 1, 3: 1}])

    # Push one carrier freq
    scan_obj_a = wexp.ScanObject('pulse', 'cf', id=1, coeff=1.0)
    scan_obj_b = wexp.ScanObject('detector', 'detection_range', id=0, coeff=1.0)
    scan_range_a = [0.0001 * i for i in range(101)]
    scan_a = wexp.SpecScan(scan_objs=(scan_obj_a, scan_obj_b), range=scan_range_a)

    return wexp.VibExperiment(field=field_a, detector=detector_a, scans=(scan_a,), magn_conditions=((-1, 2),),)

def evv_experiment_pulse_1_and_2_coincident() -> ws_experiment.experiment_abstractions.VibExperiment:
    """
    Returns VibExperiment instance for an EVV experiment variant with pulse 1 and 2 coincident in time
    """
    import wilson_suite.wilson_experiment.experiment_abstractions as wexp

    pulse_ir_1 = wexp.make_impulsive_gaussian_pulse(tc=100.0, cf=0.0, cf_uv=0.0,
                                                    maxstr=1.0e-5, wv=(0.0, 0.0, 1.0), pol=(1.0, 0.0, 0.0), id=1)

    pulse_ir_2 = wexp.make_impulsive_gaussian_pulse(tc=100.0, cf=0.0, cf_uv=0.0,
                                                    maxstr=1.0e-5, wv=(0.0, 0.0, 1.0), pol=(1.0, 0.0, 0.0), id=2)

    pulse_uvvis_1 = wexp.make_impulsive_gaussian_pulse(tc=120.0, cf=0.0, cf_uv=0.072,
                                                       maxstr=1.0e-5, wv=(0.0, 0.0, 1.0), pol=(1.0, 0.0, 0.0), id=3)

    pulses = (pulse_ir_1, pulse_ir_2, pulse_uvvis_1)

    field_a = wexp.ElectricField(pulses)

    detector_a = wexp.SpecDetector(detection_method='freq',
                                   detector_location=(0.0, 0.0, 1.0),
                                   detection_polarization=(1.0, 0.0, 0.0),
                                   detection_range=[0.003 + 0.0001 * i for i in range(101)],
                                   wv_filter=[{1: -1, 2: 1, 3: 1}])

    # Push one carrier freq
    scan_obj_a = wexp.ScanObject('pulse', 'cf', id=1, coeff=1.0)
    scan_obj_b = wexp.ScanObject('detector', 'detection_range', id=0, coeff=1.0)
    scan_range_a = [0.0001 * i for i in range(101)]
    scan_a = wexp.SpecScan(scan_objs=(scan_obj_a, scan_obj_b), range=scan_range_a)

    return wexp.VibExperiment(field=field_a, detector=detector_a, scans=(scan_a,), magn_conditions=((-1, 2),),)

def experiment_beta_alpha_cars() -> ws_experiment.experiment_abstractions.VibExperiment:
    """
    Returns VibExperiment instance for a CARS-like experiment where the (first) coherence is formed by a
    hyper-Raman-like process and the emitted light further results from a Raman-like process
    """
    import wilson_suite.wilson_experiment.experiment_abstractions as wexp


    pulse_uvvis_1 = wexp.make_impulsive_gaussian_pulse(maxstr=1.0e-5, tc = 10.0, cf=0.0,
                                                       cf_uv=0.072, wv=(0.0, 0.0, 1.0), pol=(1.0, 0.0, 0.0), id=1)

    pulse_uvvis_2 = wexp.make_impulsive_gaussian_pulse(maxstr=1.0e-5, tc = 10.0, cf=0.0,
                                                       cf_uv=0.072, wv=(0.0, 0.0, 1.0), pol=(1.0, 0.0, 0.0), id=2)

    pulse_uvvis_3 = wexp.make_impulsive_gaussian_pulse(maxstr=1.0e-5, tc = 10.0, cf=0.0,
                                                       cf_uv=0.144, wv=(0.0, 0.0, 1.0), pol=(1.0, 0.0, 0.0), id=3)

    pulse_uvvis_4 = wexp.make_impulsive_gaussian_pulse(maxstr=1.0e-5, tc = 120.0, cf=0.0,
                                                       cf_uv=0.072, wv=(0.0, 0.0, 1.0), pol=(1.0, 0.0, 0.0), id=4)


    pulses = (pulse_uvvis_1, pulse_uvvis_2, pulse_uvvis_3, pulse_uvvis_4)

    field_a = wexp.ElectricField(pulses)
    order = len(pulses)

    detector_a = wexp.SpecDetector(detection_method='freq',
                                   detector_location=(0.0, 0.0, 1.0),
                                   detection_polarization=(1.0, 0.0, 0.0),
                                   detection_range=[0.003 + 0.0001 * i for i in range(101)],
                                   wv_filter=[{1: 1, 2: 1, 3: -1, 4: 1}])

    # Push one carrier freq
    scan_obj_a = wexp.ScanObject('pulse', 'cf', id=1, coeff=1.0)
    scan_obj_b = wexp.ScanObject('detector', 'detection_range', id=0, coeff=1.0)
    scan_range_a = [0.0001 * i for i in range(101)]
    scan_a = wexp.SpecScan(scan_objs=(scan_obj_a, scan_obj_b), range=scan_range_a)

    return wexp.VibExperiment(field=field_a, detector=detector_a, scans=(scan_a,))

def evv_terms() -> list['VibPerturbedTerm']:
    """
    Returns EVV terms derived with wilson_derive: These terms are not translated to any specific axis choice and are
    instead expressed in terms of pulse ID interactions
    """
    import wilson_suite as ws

    return ws.derive.main.get_fully_enhanced_terms(experiment=evv_experiment())

def get_eval_ready_evv_terms():
    """
    Returns EVV terms derived with wilson_derive including translation to be expressed in terms of a choice of axes
    """
    import wilson_suite as ws

    experiment_a = evv_experiment()
    terms = ws.derive.derive.get_fully_enhanced_terms(experiment=experiment_a)

    print('valid ax combs', experiment_a.valid_axis_combs[0])
    print('len valid ax combs', len(experiment_a.valid_axis_combs[0].valid_axis_combs))
    return ws.derive.term_var_translate.translate_terms_to_axis_variables(terms, experiment_a.valid_axis_combs[0].valid_axis_combs[3])

def get_terms_from_json():
    from wilson_suite.wilson_derive.response_terms import VibPerturbedTerm
    from wilson_suite.wilson_utils.paths import SUITE_ROOT
    return VibPerturbedTerm.load_many_from_json(SUITE_ROOT+'/../terms_fuller_flat.json')

# # QC calculations/vibana parameters
# mol_system = abst_main.MolecularSystem(name='FORM', natoms=4)


# # spectrum eval/render parameters

# def spectral_grid():
#     axis1 = abst_main.SpectralAxis({1: 1})
#     axis2 = abst_main.SpectralAxis({1: 1, 2: -1})
#     start = {1: 250, 2: 100}
#     end = {1: 3850, 2: 7550}
#     spacer = {1: 230.8, 2: 230.8}
#     return abst_main.SpectralGrid({1: axis1, 2: axis2}, range_style='uniform',
#                                   start=start, end=end, spacer=spacer)

# def spec_eval_setup(spec_grid):
#     evi = {'dynrange': 500, 'Gamma': 4.7, 'diag_margin': 5., 'maxmax': None}
#     rndi = {'num_level_ticks': 15}
#     return abst_main.SpecEvalSetup(grid=spec_grid, ev_info=evi, rnd_info=rndi)


# def makeWilsonSimInstance(experiment, mol_system, vibanasetup, calc_setup, eval_prop_specify: dict, 
#                   spec_eval_setup: abst_main.SpecEvalSetup = None):
#     """
#     Should create an instance of a WilsonSimulation 
#     with a certain state based on the needs
    
#     """
#     # sim = abst_main.WilsonSimulation()

#     raise NotImplementedError('This fixture-making function is not yet implemented.')

# """
# 1. Ready for evaluation, WilsonSimulation instance should have:
#     a. self.system, self.terms, self.props, self.spec_eval_setup, self.vib_ana_setup
#     b. self.exp if evaluateFull()

# 2. Ready for rendering, WilsonSimulation instance should have:
#     a. self.spec, self.system, self.exp, self.diagn, self.name, self.spec_eval_setup

# 3. Ready to load numerical data (states, properties) (getResultsFromCalculationBatches()):
#     a. self.calc_batches, self.vib_ana_setup., self.props
#     b. [if self.vib_ana_setup.external_fill_from is not None: -- getResults() needs to call vibanalysis if no external fill]
#     c. getResultsFrom...() should not be needing any self. access? static methods?
#     d. getResults of CalculationBatch
     
# 4. Ready to makeCalculationBatches(): - this method is copying MolProps to batches where they are saved
#     a. self.props - props has all possible props. not getting calc batch for getting states? if requesting own vibana for stetes, how to get hessian?
# """