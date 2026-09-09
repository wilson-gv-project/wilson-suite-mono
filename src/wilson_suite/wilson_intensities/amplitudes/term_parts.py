from wilson_suite.wilson_derive.abstractions import ResonanceCondition, HarmOscStateSymbolic, PolProp, VibDiffTerm
from dataclasses import dataclass, field
from wilson_suite.wilson_intensities.amplitudes.numerical_abstractions import NumericalResonanceMotif
from wilson_suite.wilson_utils.prop_trivname import prop_trivname
from ...wilson_main.abstractions import MolPropsCollection
from typing import Callable
from wilson_suite.wilson_main.abstractions import VibState
from wilson_suite.wilson_utils.unit_convertor import convNu2Ene
# from collections.abc import Mapping


from typing import Mapping, Iterator
from types import MappingProxyType


import numpy as np
from typing import Tuple

import copy
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from ..amplitudes.vibene_differences import VibDiffCache, VibDiff
    from .spectrum_composition import ResLocGeoObject

@dataclass
class PropsCollection:
    """
    Collects several PolProp class instances and enables some group operations:
        get_avegaded_props - extract averaged properties
        get_non_avegaded_props - extract non-averaged properties
        get_cart_axes
        get_mode_indices
        get_total_difforder
    """
    props: list[PolProp]

    def __post_init__(self):
        self.props = tuple(self.props)

    def __iter__(self):
        yield from self.props

    def __hash__(self):
        # return hash(tuple([tuple(self.get_cart_axes()), self.get_total_difforder()]))
        return hash(tuple([tuple(self.get_cart_axes()), tuple(self.get_mode_indices())]))
    
    def __eq__(self, other):
        """
        Now depends on comparison of PolProp instances.
        Now PolProp instances are considered equal if the have the same lists of operators (ops)
            (further, equality of QOperator instances) and same differentiation order (dord)
        """
        if isinstance(other, PropsCollection):
            return all([p in other.props for p in self.props])
        return False
    
    def get_averaged_props(self):
        return PropsCollection(props=[p for p in self.props if p.ops])
    def get_non_averaged_props(self):
        return PropsCollection(props=[p for p in self.props if not p.ops])
    
    def get_cart_axes(self):
        return [op.o for p in self.props for op in p.ops]
    def get_mode_indices(self):
        groups = [p.inds if p.inds is not None else [] for p in self.props]
        return [idx for p_inds in groups for idx in p_inds]
    
    # UNUSED
    def get_mode_indices_grouped(self):
        return [p.inds if p.inds is not None else [] for p in self.props]
    
    def get_mode_indices_group_template(self):
        return [len(p.inds) if p.inds is not None else [] for p in self.props]
    
    # UNUSED
    def get_total_difforder(self):
        return sum([p.dord for p in self.props])
    
    def _set_attr_for_all_props(self, attr, value):
        for prop in self.props:
            prop.__setattr__(attr, value)
    
    def identify_avrg_motif(self):
        """
        Averaged properties motif/ID - inds will be set to None.
        Indices will be added later, when several terms with with avrg motifs are concidered
        """
        averaged = copy.deepcopy(self.get_averaged_props())
        if averaged.props:
            averaged._set_attr_for_all_props('inds', None)
            return averaged
    
    def sort(self):
        """
        non-averaged props will be in the end of the tuple
        """
        self.props = sorted(self.props, key=lambda j: j.ops[0].o if j.ops else float('inf'))
        return self

    def __repr__(self):
        inds_all = [len(p.inds) if p.inds else 0 for p in self.props]
        full_string = [f'{prop_trivname(ord_geo=inds_all[i], ord_el=len(p.ops))}{p.inds}{[i.o for i in p.ops]}_d{p.dord}' for i, p in enumerate(self.props)]
        return ' * '.join(full_string)

@dataclass
class FreqTermsCollection:
    freqterms: list[VibDiffTerm]
    
    def __post_init__(self):
        self.freqterms = tuple(self.freqterms)

    def __iter__(self):
        for freqt in self.freqterms:
            yield freqt

    def __hash__(self):
        return hash(self.freqterms)
    
    def __eq__(self, other):
        """
        Now depends on comparison of PolProp instances.
        Now PolProp instances are considered equal if the have the same lists of operators (ops)
            (further, equality of QOperator instances) and same differentiation order (dord)
        """
        if isinstance(other, FreqTermsCollection):
            return all([ft in other.freqterms for ft in self.freqterms])
        return False

    def get_vibenedenom(self):
        return FreqTermsCollection(freqterms=[ft for ft in self.freqterms if not ft.is_pert_wf_diff])
    
    def get_pert_wf_diff(self):
        return FreqTermsCollection(freqterms=[ft for ft in self.freqterms if ft.is_pert_wf_diff])
    
    def get_num_indices_vibenedenom(self):
        """
        """
        # these vibdiffterms have only sl, sr is zero
        return tuple(sorted(set(i for vd in self.get_vibenedenom() for i in vd.sl.q)))

@dataclass
class ResonanceMotif:
    """
    Collects ResonanceCondition instances into a "resonance motif"
        get_vibdiffs
        get_freq_axes
    """
    resonance_conditions: list[ResonanceCondition]
    
    def __iter__(self):
        for condition in self.resonance_conditions:
            yield condition

    def __eq__(self, other):
        if isinstance(other, ResonanceMotif):
            return self._tuplify() == other._tuplify()
        return False
    def __hash__(self):
        return hash(self._tuplify())
    
    def _tuplify(self):
        conditions = []
        for cond in self.resonance_conditions:
            new_pf = tuple(cond.pf)
            new_diff = tuple([tuple(cond.diff.sl.q), tuple(cond.diff.sr.q)])

            conditions.append(tuple([new_diff, new_pf]))
        return tuple(conditions)
    
    def to_str(self):
        """
        EVV / paper1 spectific here
        """
        strings = []
        for cond in self.resonance_conditions:
            if 'B' in cond.pf[0]:
                state = []
                for i in [cond.diff.sl.q, cond.diff.sr.q]:
                    if len(i)!=0:
                        state.append('+'.join(i))
                    else:
                        state.append('.')
                strings.append(f'{','.join(state)}')
        return ' x '.join(strings)

    # UNUSED?
    @classmethod
    def from_tuples(cls, tupleOfTuples):
        """
        motif1 = (((('a', 'b'), ('a',)), ('A',)), ((('b',), ('a',)), ('B',)))
        motif2 = (((('a', 'b'), ('a',)), ('A',)),)
        motif3 = (((('',), ('a',)), ('B',)), ((('',), ('a',)), ('A', '-B')))
        motif4 = (((('',), ('a',)), ('B',)), ((('b',), ('a',)), ('B',)))
        """
        r_conditions = []
        for rc_tuple in tupleOfTuples:
            print('rc tuple', rc_tuple)
            rc = ResonanceCondition(diff=VibDiffTerm(sl=HarmOscStateSymbolic(q=rc_tuple[0][0]),
                                                     sr=HarmOscStateSymbolic(q=rc_tuple[0][1])), pf=rc_tuple[1])
            r_conditions.append(rc)
        return cls(r_conditions)

    # UNUSED
    @classmethod
    def from_dicts(cls, res_conds_listdict: list[dict]):
        """
        motif3 = (
                  ((left-('',), right-('a',)), pert_freqs-('B',)), 
                  ((left-('',), right-('a',)), pert_freqs-('A', '-B')))

        res_conds_dict = [{'left': tuple, 'right': tuple, 'pert_freqs': tuple},
                          {'left': tuple, 'right': tuple, 'pert_freqs': tuple}]
        """
        r_conditions = []

        for rc_dict in res_conds_listdict:
            rc = ResonanceCondition(diff=VibDiffTerm(sl=HarmOscStateSymbolic(q=rc_dict['left']),
                                                     sr=HarmOscStateSymbolic(q=rc_dict['right'])), 
                                                     pf=rc_dict['pert_freqs'])
            r_conditions.append(rc)
        return cls(r_conditions)

    def __repr__(self):
        return f'{self.resonance_conditions}'
    
    def __len__(self):
        """
        Returns the number of elements in the container.
        """
        return len(self.resonance_conditions)
    
    # UNUSED
    @property
    def resonance_location_class(self, total_num_axes):
        return total_num_axes - len(self.resonance_conditions)
    
    # UNUSED
    def get_vibdiffs(self):
        return {i: cond.diff for i, cond in enumerate(self.resonance_conditions)}
    # UNUSED
    def get_freq_axes(self):
        return {i: tuple(cond.pf) for i, cond in enumerate(self.resonance_conditions)}
    
    def get_max_different_freq_axes(self):
       return set([i.strip('-') for cond in self.resonance_conditions for i in cond.pf])
    
    def get_nm_indices(self):
        return set([label for cond in self.resonance_conditions for i in cond.diff for label in i.q])


@dataclass(frozen=True)
class ParameterSet(Mapping[str, int]):
    """
    Immutable mapping of parameter label -> index value.
    """
    _parameters: Mapping[str, int]

    def __init__(self, parameters: Mapping[str, int]):
        if not isinstance(parameters, Mapping):
            raise TypeError("ParameterSet must be initialized with a mapping.")

        params = dict(parameters)

        if 'zero' not in params:
            params['zero'] = 'zero'

        object.__setattr__(
            self,
            "_parameters",
            MappingProxyType(params),
        )

    # --- Mapping interface ---

    def __getitem__(self, key: str):
        if key == '':
            key = 'zero'
        return self._parameters[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._parameters)

    def __len__(self) -> int:
        return len(self._parameters)

    # --- Equality & hashing ---

    def __eq__(self, other):
        if not isinstance(other, ParameterSet):
            return NotImplemented
        return self._parameters == other._parameters

    def __hash__(self):
        # Order-independent, value-based hash
        return hash(frozenset(self._parameters.items()))

    # def __lt__(self, other):
    #     if not isinstance(other, ParameterSet):
    #         return NotImplemented
    #     # Sort keys to ensure we compare 'a', then 'b', then 'c' 
    #     # regardless of insertion order.
    #     self_values = tuple(self[k] for k in sorted(self.keys()))
    #     other_values = tuple(other[k] for k in sorted(other.keys()))
        
    #     return self_values < other_values

    def __lt__(self, other):
        if not isinstance(other, ParameterSet):
            return NotImplemented
        
        sort_keys = ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h')
        relevant_keys = [k for k in sort_keys if k in self or k in other]

        # Sort keys to ensure we compare 'a', then 'b', then 'c' 
        # regardless of insertion order.
        self_vals = tuple(self.get(k, 0) for k in relevant_keys)
        other_vals = tuple(other.get(k, 0) for k in relevant_keys)
        
        return self_vals < other_vals

    # --- Convenience ---

    def parameter_labels(self):
        return [k for k in self._parameters if k != 'zero']

    def indices(self):
        return [v for v in self._parameters.values() if v != 'zero']

    def to_dict(self):
        return dict(self._parameters)

    def __repr__(self):
        repr_d = {k: v for k, v in self._parameters.items() if k != 'zero'}
        return f"{self.__class__.__name__}({repr_d})"

    # --- Pickle support ---
    def __getstate__(self):
        # Return plain dict instead of mappingproxy
        return {'_parameters': dict(self._parameters)}

    def __setstate__(self, state):
        object.__setattr__(self, "_parameters", MappingProxyType(state['_parameters']))

@dataclass
class VibStatesData:
    """
    Holds vib states data and can compute vib states energy differences
    """
    allstates: tuple[VibState]
    harmonic_osc_states_labels: tuple[int] = None
    number_of_nmodes: int = None
    
    def __post_init__(self):
        tmp_allstates = list(self.allstates)
        tmp_allstates.append(VibState(harm_quanta_coeffs={}, state_label='zero', energy=0.))
        self.allstates = tuple(tmp_allstates)
        
        self.allenergies_map = {i.state_label: i.energy for i in self.allstates}
        self.allstates_map = {i.state_label: i for i in self.allstates}
        self._storage = dict()

    # UNUSED
    def _fill_storage(self):
        for vlabel_a, energy_a in self.allenergies_map:
            for vlabel_b, energy_b in self.allenergies_map:
                self._storage[(vlabel_a, vlabel_b)] = convNu2Ene(energy_a - energy_b)


    def get_harmonic_osc_states(self):
        """
        i.state_label - TODO: make a convention, rules how to describe vibstates
        now i.state_label is str
        """
        harm_states_str = [i for i in self.allstates if ',' not in i.state_label and i.state_label!='zero']
        harm_states = {int(i.state_label): i.energy for i in harm_states_str if int(i.state_label) in self.harmonic_osc_states_labels}

        return dict(sorted(harm_states.items()))
    
    def get_state_by_label(self, state_label):
        if state_label in self.allstates_map:
            return self.allstates_map.get(state_label)
        else:
            raise ValueError(f'Requested state label - {state_label} - is not in VibStatesData')
    
    # UNUSED
    def get_energy_by_label(self, state_label):
        if state_label in self.allstates_map:
            return self.allenergies_map.get(state_label)
        else:
            raise ValueError(f'Requested state label - {state_label} - is not in VibStatesData')

@dataclass(frozen=True)
class EvaluationDataAndConfigs:
    props_data: MolPropsCollection = None
    vibstates_data: 'VibStatesData' = None
    number_of_nmodes: int = None
    nm_inds_choices: list[int] = None
    pulse_polarization_vector: list = None
    nc_sqrt_eigval: dict = None


@dataclass()
class PrecalculatedData:
    vibdiff_cache: 'VibDiffCache' = None
    avrg_tensors: dict = None
    avrg_expr_tensor_mapping: dict = None
    vibenedenoms_tensors: dict = None

@dataclass(frozen=True)
class TermParametersChoice:
    """
    Minimal representation of a group of terms that share a resonance motif.
    Each term is identified only by its integer term_id.

    in compile_feature() 
    """
    res_motif: "ResonanceMotif"
    states_parameters: Tuple["ParameterSet", ...]
    term_ids: tuple[int|str, ...] = field(default_factory=tuple)

    def sort_parameters(self) -> "TermParametersChoice":
        """
        Returns a new TermParametersChoice where the states_parameters 
        are sorted according to the ParameterSet comparison logic.
        """
        new_params = tuple(sorted(self.states_parameters))
        
        # Replace the current tuple with the sorted one
        # Using dataclasses.replace for frozen instances
        from dataclasses import replace
        return replace(self, states_parameters=new_params)


    def __hash__(self):
        return hash((self.term_ids, self.states_parameters))

    def __eq__(self, other):
        return (
            isinstance(other, TermParametersChoice)
            and self.term_ids == other.term_ids
            and self.states_parameters == other.states_parameters
        )

    # def __lt__(self, other):
    #     if len(self.term_ids) != len(other.term_ids):
    #         return len(self.term_ids) < len(other.term_ids)
    #     else:
    #         self.states_parameters
    #     return

    def __lt__(self, other):
        if not isinstance(other, TermParametersChoice):
            return NotImplemented
        
        # 1. Compare lengths of term_ids
        if len(self.term_ids) != len(other.term_ids):
            return len(self.term_ids) < len(other.term_ids)
        
        # 2. Compare term_ids values (lexicographical)
        if self.term_ids != other.term_ids:
            return self.term_ids < other.term_ids
            
        # 3. Compare the sequences of ParameterSets
        # Python will compare self.states_parameters[0] < other.states_parameters[0], etc.
        return self.states_parameters < other.states_parameters
# -------------------------------------------------------


# UNUSED
def is_tuple_of_tuples(my_variable):
    if not isinstance(my_variable, tuple):
        return False

    for item in my_variable:
        if not isinstance(item, tuple):
            return False

    return True



def safe_arange_inclusive_scaled(start, stop, step):
    """
    google ai

    Generates a range by scaling to integers, including the stop value. 
    Best for decimal steps like 0.1, 0.01 etc.
    """
    # Scale inputs to integers
    start_int = round(start / step)
    stop_int = round(stop / step)
    
    # Create the integer range, going one step further than 'stop_int' to guarantee inclusion
    int_range = np.arange(start_int, stop_int + 1) 
    
    # Scale back down
    return int_range * step

def linspace_with_step(start, stop, step):
    """
    Generates a range using np.linspace that respects exact start/stop points 
    while accepting a step size input. Includes the stop point.
    """
    # Calculate the number of intervals required using integer math logic 
    # to avoid floating point accumulation errors
    num_intervals = round((stop - start) / step)
    
    # The number of points is the number of intervals + 1 (fencepost error principle)
    num_points = int(num_intervals) + 1
    
    # Use linspace, which is precise with boundaries, now that we have the exact num_points
    return np.linspace(start, stop, num=num_points, endpoint=True)



@dataclass
class EvalTerm:
    """
    Parametrized term - VibPertTerm but with abc parameters specified.

    resonance_function - resonance part parametrized for this term
    """
    amplitude_coeff: float
    res_loc: 'ResLocGeoObject'
    # function of grid → value
    compiled_res_motif: NumericalResonanceMotif # needs states data to be compiled; can be evaluated with meshgrids
    parameters: dict

    @classmethod
    def from_VibPertTerm():
        
        return


@dataclass
class EvalFeature:
    location: dict  # {'A':1100,'B':2300}
    terms: list[EvalTerm]
    amplitude: float
    lineshape_param: float





def make_resonance_function(res_motif: ResonanceMotif, 
                            meshgrids: dict[str, np.ndarray]) -> Callable:
    """
    take symbolic ResonanceMotif and make a function for 
    """
    pfreqs_for_res_motif = {}
    
    for rc in res_motif:
        pfreqs_for_res_motif[rc] = sum([meshgrids[ax.strip('-')] * rc.pf_dict[ax.strip('-')] for ax in rc.pf])
    
    def param_func(param_set: ParameterSet, vibstates_data: VibStatesData, 
                   vibdiff_cache: 'VibDiffCache', lineshape_parameter: float):
        rs_difs_num = []
        
        for rc in res_motif:
            # for each resonace condition (rc) compute vibdiff part
            vd = VibDiff.from_symbolic(rc.diff, param_set, vibstates_data)
            vd.cache_it(vibdiff_cache)
            # for this resonance condition get grids for axes

            # in reciprocal centimeters
            rc_num = vd.energy_difference(au=False) - pfreqs_for_res_motif[rc] - 1j*lineshape_parameter
            
            # in au
            rs_difs_num.append(convNu2Ene(rc_num))
        return 
    
    return param_func