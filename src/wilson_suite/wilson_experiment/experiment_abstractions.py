from dataclasses import dataclass
from typing import Optional, Iterable
from operator import itemgetter
import copy
from math import inf as infinity

# TODO: Expand functionality according to below TODOs

@dataclass
class SpecDetector:
    """
    Class to represent a spectral detector

    ----
    detection_method: String: A detection method: "integrated", "time", "freq"
    If detection_method is "time" or "freq", then the detection data is one spectral dimension
    If detection_method is "integrated", then the detection data is a scalar
    Currently, only "freq" (frequency-range) detection is supported.

    detector_location: List of floats: Unit vector describing the direction in which the detector is facing.
    Default: [0.0, 0.0, 1.0]. Currently only used for defining polarization filtering.

    detection_polarization: List of floats: Detect only light with this specific polarization vector. Default:
    [1.0, 0.0, 0.0].

    detection_range: List of floats: For "time" or "freq" detection, tell over which points (the range)
    in either t/E space as relevant the data is collected

    wv_filter: List of dictionaries {pulse label: sign, ...}:  Detect only light along this/these particular
    phase-matching direction(s)

    ignore_collinear: If using a wavevector filter, ignore other effects collinear with this/these direction(s)?
    Currently not used.
    """
    detection_method: str
    
    detector_location: Optional[tuple[float]] = None
    
    detection_polarization: Optional[tuple[float]] = None
    
    # Comment: detection_range as None and detection_method as 'freq' is valid but results in no dimensionality
    detection_range: Optional[list[float]] = None
    wv_filter: Optional[list[dict]] = None
    ignore_collinear: bool = True

    overall_phase: complex = 1.0 + 0.0j

    def __post_init__(self):
        if self.detection_method not in {'time', 'freq', 'int'}:
            raise ValueError("The detection type must be either 'time', 'freq'(uency), or 'int'(egrated)")

        if not self.overall_phase == 1.0 + 0.0j:
            raise ValueError('Detector overall phase currently restricted to zero shift')

        # Will currently not be reach because of restriction to zero shift, but will be relevant when that is lifted
        if not(abs(self.overall_phase) - 1.0 < 1e-10):
            raise ValueError('Detector overall phase must be of unit length')

@dataclass
class ScanObject:
    """
    Class to represent some attribute of the experiment which could be scanned. Currently not in use.

    category: String: The (main) category of the object to be scanned, e.g. "pulse" or "detector"
    subcategory: String: The subcategory of the object to be scanned, e.g. "cf" (carrier frequency)
    id: Integer: An integer identifier for the specific instance of the main category, e.g. if 'category' is "pulse"
        (of which there are typically several, labelled by integer indices), then the present 'id' attribute
        tells which of these pulses are meant in this scan object
    coeff: Float: Coefficient: As a scan is performed across a range (see SpecScan), 'coeff' tells which scaling to use
    for the scan range increments when applying them to the present scan object. For example, one might want to have a
    scan that increases one parameter according to the scan range while concurrently decreasing another parameter twice
    as rapidly. This can be represented by associating the scan with two scan objects referring to the respective
    parameters, where the first has 'coeff' == 1.0 and the second has 'coeff' == -2.0. Default: 1.0.
    """

    category: str
    subcategory: str
    id: int = 0
    coeff: float = 1.0


    def __post_init__(self):

        # Testing against currently recognized scan categories/subcategories.
        # Presence in these lists is required but not sufficient to
        # indicate actual support for a specific scan.

        valid_scan_objs = ['pulse', 'detector']
        valid_scan_attributes = {'pulse': ['cf', 'tc', 'dev'], 'detector': ['detection_range']}

        if not self.category in valid_scan_objs:
            raise ValueError('Scan category not supported')
        if not self.subcategory in valid_scan_attributes[self.category]:
            raise ValueError('Scan subcategory not supported')


@dataclass
class SpecScan:
    """
    Class to represent a spectral scan (adding to the dimensionality of a spectrum)

    scan_objs: Tuple of ScanObject instances: Tells what this scan will vary

    range: Iterable over which scan objects are varied (scaled by their multipliers as represented by
    their 'coeff' attributes)
    """

    scan_objs: tuple[ScanObject]
    range: Iterable

    def __post_init__(self):

        from collections.abc import Iterable

        if not (isinstance(self.scan_objs, tuple)):
            raise TypeError('scan_objs must be a tuple of ScanObject instances')
        else:
            for i in self.scan_objs:
                if not (isinstance(i, ScanObject)):
                    raise TypeError('scan_objs must be a list of ScanObject instances')

        if not isinstance(self.range, Iterable):
            raise TypeError('range is not iterable')

@dataclass
class EmPulse:
    """
    Class to represent an electromagnetic pulse
    
    ----
    env: String: Pulse time-domain envelope type: The only current valid choice is "gaussian".
    maxstr: Float: Pulse amplitude at maximum of envelope. Default: Zero strength
    tc: Float: Point in time at which pulse envelope is at maximum
    cf: float: (Infrared-range) carrier frequency
    cf_uv: float: Designated "UV/VIS range" part of carrier frequency (for e.g. CARS-style cancellation).
        - cf_uv should be specified as a nonnegative value: Any cancellation should follow from the phase-matching
        wavevector--frequency combination in the experiment
        - For pulses where cf_uv != 0.0, then cf must be 0.0
    dev: float: Deviation parameter (e.g. broadness of Gaussian pulse). An impulsive-like pulse can be considered as the
    dev -> 0.0 limit of a Gaussian pulse. A "continuous-wave"-like pulse can be considered as the dev -> infty limit of
    a Gaussian pulse.
    wv: floats: Unit wavevector propagation direction with respect to laboratory axes
    pol: floats: Polarization: Unit vector describing polarization direction with respect to laboratory axes
        - Must be orthogonal to wavevector
        - Only linear polarization currently supported (no phase difference between orthogonal components of
        polarization vector in plane of polarization)
        - Default: (1.0, 0.0, 0.0)
    overall_phase: complex number defining a unit vector in the complex plane: Overall phase of pulse. Currently enforced as (1.0, 0.0)

    id: integer: Pulse ID label
    """
    env: str
    
    tc: float = None
    cf: float = None
    dev: float = None

    cf_uv: float = 0.0
    maxstr: float = 0.0

    wv: tuple[float] = (0.0, 0.0, 1.0)
    pol: tuple[float] = (1.0, 0.0, 0.0)
    overall_phase: complex = 1.0 + 0.0j

    id: int = None

    def __post_init__(self):
        
        # TODO: Generalize to arbitrary pulses and chirped pulses
        allowed_envelopes = ['gaussian']

        if self.env not in allowed_envelopes:
            raise ValueError('The only current recognized pulse envelope choices is "gaussian"')

        if not (self.cf_uv == 0.0):
            if not (self.cf == 0.0 or self.cf == None):
                raise ValueError('Pulses with non-zero UV/VIS carrier freqs. must have IR carrier freq part set to zero or None')

        if not (self.cf_uv >= 0.0):
            raise ValueError('A carrier frequency must be nonnegative')

        if self.cf is not None:
            if not (self.cf >= 0.0):
                raise ValueError('A carrier frequency must be nonnegative')

        if self.env == 'gaussian':
            if self.cf is None:
                if not self.dev == 0.0:
                    raise AssertionError('A non-impulsive-tending Gaussian pulse must have a carrier frequency')
            if self.tc is None:
                if not self.dev == infinity:
                    raise AssertionError('A non-continuous-wave-tending Gaussian pulse must have a parameter describing the time at which the temporal envelope is at its maximum')
            if self.dev is None:
                raise AssertionError('A Gaussian pulse must have a time deviation parameter')

        # Wavevector: In which unit vector direction is the pulse wave travelling
        if isinstance(self.wv, tuple):
            if len(self.wv) == 3:
                if all([isinstance(i, float) for i in self.wv]):
                    wv_len = (self.wv[0]**2.0 + self.wv[1]**2.0 + self.wv[2]**2.0)**0.5
                    if not wv_len == 1.0:
                        print('Wavevector was normalized')
                    self.wv = tuple([i/wv_len for i in self.wv])

                else:
                    raise AssertionError('The pulse wavevector must be a len 3 tuple of floats')
            else:
                raise AssertionError('The pulse wavevector must be a len 3 tuple of floats')
        else:
            raise AssertionError('The pulse wavevector must be a len 3 tuple of floats')

        # Polarization: Specify the polarization of the pulse
        # Currently supports unit linear polarization
        if isinstance(self.pol, tuple):
            if len(self.pol) == 3:
                if all([isinstance(i, float) for i in self.pol]):
                    pol_len = (self.pol[0]**2.0 + self.pol[1]**2.0 + self.pol[2]**2.0)**0.5
                    if not pol_len == 1.0:
                        print('Wavevector was normalized')
                    self.pol = tuple([i/pol_len for i in self.pol])

                    wv_pol_dot = self.pol[0] * self.wv[0] + self.pol[1] * self.wv[1] + self.pol[2] * self.wv[2]

                    if not(wv_pol_dot == 0.0):
                        raise AssertionError('Error: Wavevector of pulse not orthogonal to polarization vector')

                else:
                    raise AssertionError('The polarization vector must be a len 3 tuple of floats')
            else:
                raise AssertionError('The polarization vector must be a len 3 tuple of floats')
        else:
            raise AssertionError('The polarization must be a len 3 tuple of floats')

        if not self.overall_phase == 1.0 + 0.0j:
            raise AssertionError('Overall phase currently restricted to zero shift')

        if not((abs(self.overall_phase) - 1.0) < 1e-10):
            raise AssertionError('The overall phase must be of unit length')

    def tendsImpulsive(self):
        """
        Does this pulse tend to impulsive?
        """

        if self.env == 'gaussian':
            if self.dev == 0.0:
                return True
            return False
        else:
            raise ValueError('"Tends-impulsive" check currently not implemented for non-Gaussian pulses')

    def tendsContinuous(self):
        """
        Does this pulse tend to being continuous-wave-like?
        """

        if self.env == 'gaussian':
            if self.dev == infinity:
                return True
            return False
        else:
            raise ValueError('"Tends-impulsive" check currently not implemented for non-Gaussian pulses')


# FIXME: Here and next two fns: Change to be in terms of f(*, kw1=kw1, ...) style
def make_gaussian_pulse(tc: float, cf: float, dev: float, cf_uv: float = 0.0, maxstr: float=0.0,
                        wv: tuple[float] = (0.0, 0.0, 1.0), pol: tuple[float] = (1.0, 0.0, 0.0),
                        overall_phase: complex = 1.0 + 0.0j, id: int = None):
    """
    Helper function: Makes a Gaussian EmPulse instance with the selected parameters.
    See EmPulse for a definition of these parameters.

    Required arguments: tc, cf, dev.
    """

    return EmPulse(env = 'gaussian', tc = tc, cf = cf, dev = dev, cf_uv = cf_uv, maxstr = maxstr, wv = wv,
                   pol = pol, overall_phase = overall_phase, id = id)


def make_impulsive_gaussian_pulse(tc: float, cf: float = None, cf_uv: float = 0.0, maxstr: float=0.0,
                         wv: tuple[float] = (0.0, 0.0, 1.0), pol: tuple[float] = (1.0, 0.0, 0.0),
                         overall_phase: complex = 1.0 + 0.0j, id: int = None):
    """
    Helper function: Makes an impulsive-tending Gaussian EmPulse instance with the selected parameters.
    See EmPulse for definitions of these parameters. The present function does not accept a dev parameter
    because it is zero for an impulsive pulse and this is carried out explicitly in the call below. All
    other parameters are transferred directly to the EmPulse instance creation.

    Required arguments: tc. All other arguments optional.
    """

    return EmPulse(env = 'gaussian', tc = tc, cf = cf, dev = 0.0, cf_uv = cf_uv, maxstr = maxstr, wv = wv,
                   pol = pol, overall_phase = overall_phase, id = id)

def make_cw_gaussian_pulse(cf: float, tc: float = None, cf_uv: float = 0.0, maxstr: float=0.0,
                         wv: tuple[float] = (0.0, 0.0, 1.0), pol: tuple[float] = (1.0, 0.0, 0.0),
                         overall_phase: complex = 1.0 + 0.0j, id: int = None):
    """
    Helper function: Makes a continuous-wave-tending Gaussian EmPulse instance with the selected parameters.
    See EmPulse for definitions of these parameters. The present function does not accept a dev parameter
    because it is infinity for a cw pulse and this is carried out explicitly in the call below. All
    other parameters are transferred directly to the EmPulse instance creation.

    Required arguments: cf. All other arguments optional. NOTE changed ordering of cf and tc arguments compared to
    other make_..._pulse functions because of this optionality.
    """

    return EmPulse(env = 'gaussian', tc = tc, cf = cf, dev = infinity, cf_uv = cf_uv, maxstr = maxstr, wv = wv,
                   pol = pol, overall_phase = overall_phase, id = id)

# The field consists of a collection of pulses
@dataclass
class ElectricField:
    """
    Class to represent an electromagnetic field consisting of one or more pulses
    
    ----
    pulses: List of EmPulse instances: The pulses making up the field. While EmPulse itself does not require an ID to
    be specified, the use of pulses in ElectricField must have each pulse be assigned an ordinal integer ID starting from 1
    """

    pulses: tuple[EmPulse]

    def __post_init__(self):

        pulse_id_target = [i + 1 for i in range(len(self.pulses))]
        for i in self.pulses:
            if not i.id in pulse_id_target:
                raise ValueError('Pulse ID', i.id, ' of field not found in ordinal list')
            else:
                pulse_id_target.remove(i.id)

        # This condition should never be met but just to be safe
        if not pulse_id_target == []:
            raise AssertionError('Collection of pulse IDs do not correspond to ordinal list')

@dataclass
class VibExperiment:
    """
    Class to represent a vibrational wave-mixing experiment

    field: ElectricField instance: A "base" perturbing field (upon which scans may be imposed)

    detector: SpecDetector instance: The detector for this experiment

    scans: List of SpecScan instances: Tells which parameters will be scanned over (and how) in this experiment

    magn_conditions: Tuple of tuples: Magnitude conditions for use in identifying terms that will not become
    fully resononant in this experiment. Format: Outer tuple collects magnitude conditions. Each inner tuple is
    a magnitude condition and consists of signed pulse references (NOTE: Currently not using the SignedPulseTuple class)
    where the sum of the associated frequencies are understood to be significantly > 0, where "significantly > 0" means
    "never close to zero".
    Example: ( (-1, 2), (2, 3, -4) ) denotes two magnitude conditions:
        a) -w1 + w2 is always significantly > 0,
        b) w2 + w3 - w4 is always significantly > 0
    """

    field: ElectricField
    detector: SpecDetector
    scans: tuple[SpecScan] = None
    magn_conditions: tuple[tuple] = None

    def __post_init__(self):

        if not isinstance(self.field, ElectricField):
            raise TypeError('The field attribute must be an ElectricField instance')

        if not isinstance(self.detector, SpecDetector):
            raise TypeError('The detector attribute must be a SpecDetector instance')

        if self.scans is not None:
            if not isinstance(self.scans, tuple):
                raise TypeError('The scans attribute, if specified, must be a tuple of SpecScan instances')
            for i in self.scans:
                if not isinstance(i, SpecScan):
                    raise TypeError('The scans attribute, if specified, must be a tuple of SpecScan instances')

        if self.magn_conditions is not None:
            if not isinstance(self.magn_conditions, tuple):
                raise TypeError('The magn_conditions attribute, if specified, must be a tuple of tuples')
            for i in self.magn_conditions:
                if not isinstance(i, tuple):
                    raise TypeError('The magn_conditions attribute, if specified, must be a tuple of tuples')


        from wilson_suite.wilson_experiment.indep_vars_and_axes import (PhaseMatchingCondition, SignedPulseTuple,
                                                                        find_indep_exp_variables, find_valid_axes,
                                                                        find_canonical_axes)

        # Establishes an assumption: The order is the same as the number of pulses in the field
        # It is furthermore (but not strictly from this) assumed that in the experiment, the system will interact
        # once with each pulse
        self.order = len(self.field.pulses)

        # Determine which phase-matching condition(s) will come under consideration in this experiment
        relevant_phasematch = []

        # If no specified phase-matching (wavevector) filter, all are (potentially) relevant
        if self.detector.wv_filter is None:

            from itertools import product as iter_prod
            k = 0

            for i in iter_prod([1, -1], repeat=len(self.field.pulses)):

                new_phasematch = []

                for j in range(len(self.field.pulses)):
                    new_phasematch.append(self.field.pulses[j].id * i[j])

                relevant_phasematch.append(PhaseMatchingCondition(SignedPulseTuple(tuple(new_phasematch)), k))
                k += 1

        # Otherwise, registered only that/those specified for the detector
        else:

            k = 0

            for i in range(len(self.detector.wv_filter)):

                new_phasematch = []

                for j in self.detector.wv_filter[i]:
                    new_phasematch.append(j * self.detector.wv_filter[i][j])

                relevant_phasematch.append(PhaseMatchingCondition(SignedPulseTuple(tuple(new_phasematch)), k))
                k += 1

        self.relevant_phasematch = relevant_phasematch

        # FIXME: Replace with try...except in case not sufficient data specified
        try:
            self.dim = self.findDimensionality()
        except AssertionError:
            self.dim = None

        # Determine pulse "epochs" - i.e. disjoint (or taken to be disjoint) time partitions of the pulses
        self.epochs = find_epochs(self.field)

        # Determine which interaction orderings are causally possible
        self.int_sequences = self.findInteractionSequences()

        # Register UV/VIS range carrier frequencies of field for convenience
        self.cfuv = get_carrier_freqs_uv(self.field.pulses)

        # Find valid choices of independent variables
        self.indep_vars = find_indep_exp_variables(self.field.pulses, self.epochs, self.relevant_phasematch)

        # Find valid choices of spectral axes given the choices of independent variables determined above
        self.valid_axis_combs = find_valid_axes(self.indep_vars)

        # If no canonical axes can be determined, set to None
        try:
            self.canonical_axes = find_canonical_axes(self.indep_vars)
        except ValueError:
            self.canonical_axes = None

        # Register all polarization vectors (associated with the detector and pulses) for convenience
        # Here I establish a convention: Macroscopic ranks are with respect to pulse IDs but first rank refers to the
        # detected signal (so detected, pulse ID 1, pulse ID 2, ...)
        all_polarizations = [copy.deepcopy(self.detector.detection_polarization)]

        # Could probably be done more elegantly but works
        for i in range(len(self.field.pulses)):
            for j in self.field.pulses:
                if j.id == i + 1:
                    all_polarizations.append(copy.deepcopy(j.pol))

        self.all_polarizations = all_polarizations

        # Determine the macroscopic orientational average polarization vector
        from wilson_suite.wilson_intensities.amplitudes.averaging import get_pol_laser
        self.polarization_avg_vector = get_pol_laser(self.all_polarizations)

    def tell_axis_options(self):

        for i in range(len(self.valid_axis_combs)):
            self.valid_axis_combs[i].present_spectral_axis_choices(from_exp_index = i)

    def is_axis_set_by_ref_valid(self, ref: dict):
        """
        TODO: Take a reference describing an axis set choice in terms of
        the dict {label 1: ((indep var 1_1), (indep var 1_2 ) ...), label 2: ((indep var 2_1), ...) , ...}
        and return True if it's a valid choice for this experiment (TODO: ...phase-matching combination, indep var choice?)
        or False if it's not

        Sketch: Similar to choose_axis_set_by_ref

        """
        pass

    def choose_axis_set_by_ref(self, ref: dict):
        """
        TODO: Take a reference describing an axis set choice in terms of
        the dict {label 1: ((indep var 1_1), (indep var 1_2 ) ...), label 2: ((indep var 2_1), ...) , ...}
        and return the corresponding SpectralAxisSet instance

        Sketch: Look through self.valid_axis_choices and find out if any of the registered valid choices
        match: If yes, return that SpectralAxisSet instance - otherwise raise an error
        """
        pass

    def findDimensionality(self) -> int:
        """
        Using detector and scans information, determine the dimensionality of the spectral data
        that carrying out this experiment would produce

        Returns an integer d telling this dimensionality
        """

        d = 0
        d += len(self.scans)

        # Time or frequency series adds a dimension
        if self.detector.detection_method == 'time':
            if self.detector.detection_range is not None:
                return d + 1

        elif self.detector.detection_method == 'freq':
            if self.detector.detection_range is not None:
                return d + 1

        # Integrated intensity does not add a dimension
        elif self.detector.detection_method == 'int':
            return d

        else:
            raise AssertionError('Cannot determine dimensionality: Unrecognized detection method')

        return d
    
    def findInteractionSequences(self) -> list:
        """
        Based on the experiment information, find out if there must be a specific sequence/sequences of interactions with pulses

        Returns a list [[{sequence 1 interaction 1: pulse i}, {seq. 1 int. 2: pulse j}, ...],
                        [{seq. 2 int. 1: pulse k}, ... ], ...]
        """

        def interactionRecurse(res: list, curr_int: list, rem_wv: dict, curr_epoch: int, epochs: list):
            """
            Tail-recursive routine for finding interaction sequences

            res: List of lists: Results accumulator
            curr_int: List of dictionaries: Result currently being assembled
            rem_wv: Dictionary: One wavevector filter dictionary
            curr_epoch: Epoch counter
            epochs: List of epochs as determined by ElectricField.findEpochs
            """


            # Termination condition
            # If this interaction sequence satisfied the wavevector filter, append it
            if rem_wv == {}:
                res.append(tuple(curr_int))

            # Recursion
            else:
                # Can one or more pulses in requested wv be found at the current or later epoch?
                # If so, make all combinations, update rem wv and recurse further at same epoch

                for t in range(curr_epoch, len(epochs)):
                    for i in rem_wv:

                        if i in epochs[t]:

                            new_rem_wv = copy.deepcopy(rem_wv)
                            new_int = copy.deepcopy(curr_int)
                            new_int.append({i: new_rem_wv[i]})

                            del new_rem_wv[i]

                            interactionRecurse(res, new_int, new_rem_wv, t, epochs)



        if self.detector.wv_filter is None:
            raise AssertionError('Interaction sequence determination currently only implemented for wavevector filter detector')

        if len(self.detector.wv_filter) > 1:
            raise AssertionError('Interaction sequence determination currently not supported for more than one phase-matching direction')

        int_sequences = []
        int_seed = []

        for i in self.detector.wv_filter:

            interactionRecurse(int_sequences, int_seed, i, 0, find_epochs(self.field))

        return int_sequences
    

    def derive_terms(self) -> dict:
        """
        returns original derived terms with pulses (independent vars)
        """
        from wilson_suite.wilson_derive.derive import get_fully_enhanced_terms
        return get_fully_enhanced_terms(experiment=self)


def get_carrier_freqs_uv(pulses) -> dict:
    """
    Get dictionary of UV/VIS-range part of carrier frequencies

    Returns: Dictionary {pulse 1: UV/VIS carrier freq., ...}
    """
    cfuv_dict = {}
    for i in pulses:
        cfuv_dict[i.id] = i.cf_uv

    return cfuv_dict

def find_epochs(field, tol: float=0.0) -> list:
    """
    Divide field into epochs with either zero or finite tolerance
    Currently only supported for a field consisting of ideal or impulsive pulses

    tol: Float: Tolerance for non-temporal coincidence (currently not supported)
    # TODO: Add support for tolerance

    Returns: List of lists: [[epoch 1 pulse 1, epoch 1 pulse 2, ...], [epoch 2 pulse 1, ...], ...]
    """

    if not(tol == 0.0):
        raise ValueError('Non-zero tolerance not yet supported in find_epochs')

    for i in field.pulses:
        if not i.tendsImpulsive():
            raise AssertionError('Can currently only determine epochs for fields with impulsive-tending pulses')
        if i.id is None:
            raise AssertionError('All pulses must have IDs for valid epoch determination')

    times_ids = sorted([(i.tc, i.id) for i in field.pulses], key=itemgetter(0))
    epochs = [[]]
    epoch = 0
    curr_time = times_ids[0][0]

    for i in times_ids:
        if not(i[0] == curr_time):
            epochs.append([])
            epoch += 1
            curr_time = i[0]
        epochs[epoch].append(i[1])

    sorted_epochs = []

    for i in epochs:
        sorted_epochs.append(sorted(i))

    return sorted_epochs

def uv_cancels(coll: tuple, cfs_uv: dict, tol: float=1e-10) -> bool:
    """
    Do the UV/VIS frequency components of this collection of pulses cancel?

    coll: tuple: (Signed) pulse ID references
    cfs_uv: dictionary {non-signed pulse ID: non-signed UV/VIS frequency component, ...}
    tol: Tolerance: If the magnitude of the requested combination of freq components is beneath tol,
    then accept as cancelling

    return: True if cancellation was determined, False otherwise

    # FIXME: This routine is a bit confusingly made but should
    """


    if tol < 0.0:
        raise ValueError('The tolerance must be a nonnegative number')

    acc = 0.0

    for i in coll:

        sgn = (i > 0) - (i < 0)

        acc += sgn * cfs_uv[sgn * i]

    sgnacc = (acc > 0) - (acc < 0)

    return ((sgnacc * acc) <= tol)