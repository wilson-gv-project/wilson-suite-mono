"""
===
purification of abstractions
===

rsp_evaluator


"""

import copy
from collections.abc import Callable, Iterator, Mapping, Sequence
from dataclasses import InitVar, dataclass, field
from types import MappingProxyType
from typing import TYPE_CHECKING, Any

import numpy as np

from wilson_suite.wilson_derive.abstractions import (
    HarmOscStateSymbolic,
    PolProp,
    ResonanceCondition,
    VibDiffTerm,
)
from wilson_suite.wilson_utils.prop_trivname import prop_trivname
from wilson_suite.wilson_utils.unit_convertor import convNu2Ene

if TYPE_CHECKING:
    from wilson_suite.wilson_derive.response_terms import VibPerturbedTerm

# symbolic
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
    props: Sequence[PolProp]

    def __post_init__(self):
        self.props = tuple(self.props)

    def __iter__(self):
        yield from self.props

    def __hash__(self):
        return hash( (self.get_cart_axes(), self.get_mode_indices()) )
    
    def __eq__(self, other):
        """
        Now depends on comparison of PolProp instances.
        Now PolProp instances are considered equal if the have the same lists of operators (ops)
            (further, equality of QOperator instances) and same differentiation order (dord)
        """
        if isinstance(other, PropsCollection):
            return all(p in other.props for p in self.props)
        return False
    
    def get_averaged_props(self):
        return PropsCollection(props=[p for p in self.props if p.ops])
    def get_non_averaged_props(self):
        return PropsCollection(props=[p for p in self.props if not p.ops])
    
    def get_cart_axes(self):
        return tuple(op.o for p in self.props for op in p.ops)
    def get_mode_indices(self):
        groups = [p.inds if p.inds is not None else [] for p in self.props]
        return tuple(idx for p_inds in groups for idx in p_inds)
    
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

# symbolic
@dataclass
class FreqTermsCollection:
    """
    VibDiffTerm's states are HarmOscStateSymbolic instances
    """
    freqterms: Sequence[VibDiffTerm]
    
    def __post_init__(self):
        self.freqterms = tuple(self.freqterms)

    def __iter__(self):
        yield from self.freqterms

    def __hash__(self):
        return hash(self.freqterms)
    
    def __eq__(self, other):
        """
        Now depends on comparison of PolProp instances.
        Now PolProp instances are considered equal if the have the same lists of operators (ops)
            (further, equality of QOperator instances) and same differentiation order (dord)
        """
        if isinstance(other, FreqTermsCollection):
            return all(ft in other.freqterms for ft in self.freqterms)
        return False

    def get_vibenedenom(self):
        """
        freqterm.sl. and freqterm.sr. should be HarmOscStateSymbolic instances
        """
        def sr_or_sl_only(freqterm: VibDiffTerm):
            return (freqterm.sl.q == []) or (freqterm.sr.q == [])

        return FreqTermsCollection(freqterms=[ft for ft in self.freqterms if not ft.is_pert_wf_diff or sr_or_sl_only(ft)])
    
    def get_pert_wf_diff(self):
        return FreqTermsCollection(freqterms=[ft for ft in self.freqterms if ft.is_pert_wf_diff])
    
    def get_num_indices_vibenedenom(self):
        """
        these vibdiffterms have only sl, sr is zero
        """
        return tuple(sorted({i for vd in self.get_vibenedenom() for i in vd.sl.q}))


# symbolic
@dataclass
class ResonanceMotif:
    """
    Collects ResonanceCondition instances into a "resonance motif"
        get_vibdiffs
        get_freq_axes
    """
    resonance_conditions: Sequence[ResonanceCondition]
    
    def __iter__(self):
        yield from self.resonance_conditions

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
            new_diff = (tuple(cond.diff.sl.q), tuple(cond.diff.sr.q))

            conditions.append((new_diff, new_pf))
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


# numerical
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


## ------------------------------------------------------------------
# numerical
@dataclass(frozen=True)
class MolSystemData:
    """Everything obtained externally. No configuration."""

    name: str
    states: tuple
    eigenvals: np.ndarray
    eigenvecs: np.ndarray
    mol_props: 'MolPropsCollection'
    natoms: int | None = None
    geo: Any = None
    geo_extra: Any = None
    linear: bool = False
    conformer: str = "conf1"


@dataclass(frozen=True)
class DataOriginInfo:
	"""
	Class to represent computational setups for properties obtained external to Wilson
	Does not need to pertain to an actual program and could also be used for "get from no specific calculation"/
	"get from file"

	----
	source_type: String: Options: gaussian, cfour, wilson
	lvl_theory: String: Level of theory
	basis_set: String: Basis set
	base_file_loc: String: path to the base file
	"""
	# Strings
	source_type: str = ''
	
	lvl_theory: str = ''
	basis_set: str = ''

	base_file_loc: str = ''


	def __hash__(self):
		def to_tuple(x):
			return tuple(sorted(x.items())) if isinstance(x, dict) else x
		return hash((self.source_type, self.lvl_theory, self.basis_set, to_tuple(self.base_file_loc)))

	def __eq__(self, other):
		if not isinstance(other, DataOriginInfo):
			return False
		
		return (
            self.source_type == other.source_type and
            self.lvl_theory == other.lvl_theory and
            self.basis_set == other.basis_set and
			self.base_file_loc == other.base_file_loc
		)


# numerical
@dataclass
class MolecularProperty:
	"""
	Class to represent a molecular (energy derivative or similar) property
	Can both be used "head only" (only prop_spec, target_basis, target_units) to specify only the concept of a property
    and "full" (system, calc_setup) for a particular realization (optional with/without values)
	
	----
	prop_spec: Dictionary {'attr name': val, ...}: Info like perturbing operators, frequencies etc. (all values must be hashable)
	triv_name: String: Trivial name For simplified reference
	vals: Form not specified: Values of properties - could be array or dictionary
	system: MolecularSystem instance: For which system?
	calc_setup: DataOriginInfo instance: For which calculation setup?

	see more in test_main_dataclasses.py::test_MolecularProperty
	"""
	# FIXME: Improve on prop_spec name; settle more consistently what the attributes will be and what must be default
	prop_spec: dict
	trivial_name: str | None = None
	vals: InitVar[Any] = field(default=None, repr=False)
	calc_setup: DataOriginInfo | None = None
	extra_data: dict | None = None

	def to_dict(self):
		return {
			"prop_spec": self.prop_spec,
			"trivial_name": self.trivial_name,
		}

	# FIXME: Complete/update this
	def h(self, htype: int) -> int:
		"""
		Hashing function with four hash types

		htype: integer: Hash type: Valid choices are

		1: "head only" information (only hash(prop_spec))
		2: hash involves attributes from 1) but also tgt basis, tgt units
		3: hash involves attributes from 2) but also system, calc_setup
		4: hash involves attributes from 3) but also in_basis, in_units
		# TODO: Check for adequate property specification and values format when known
		# TODO: Consider enforcing specification of units and basis when values are provided

		Returns an integer hash value
		"""

		hlist = []

		if (htype < 1) or (htype > 4):

			raise AssertionError('Property hash must be requested with type argument (1-4)')

		if htype >= 1:

			for i in self.prop_spec:

				hlist.append(i)
				hlist.append(self.prop_spec[i])

		if htype >= 3:

			hlist.append(self.system.h())
			hlist.append(self.calc_setup.h())


		return hash(tuple(hlist))

	def addSystem(self, system: MolSystemData):
		"""
		Associate a MolecularSystem instance

		system:	MolecularSystem instance: The system to be attached
		"""

		self.system = system

	def addCalcSetup(self, calc_setup):
		"""
		Associate an DataOriginInfo instance

		calc_setup: DataOriginInfo instance: The setup to be attached
		"""

		self.calc_setup = calc_setup
	
	# Add values (usually scalars or a numPy array)
	def addValues(self, values: Any):
		"""
		Associate values to this property

		values: Undetermined form: The values to be added
		"""

		self.vals = values


# numerical
@dataclass
class MolPropsCollection:
	properties: list[MolecularProperty]

	def get(self, trivial_name: str):
		d = {prop.trivial_name: prop for prop in self.properties}
		if trivial_name not in d:
			raise ValueError(f'trivial_name {trivial_name} is not in MolPropsCollection')
		return d.get(trivial_name)

	def __getitem__(self, trivial_name):
		"""Allow coll[name] syntax."""
		return self.get(trivial_name)

	def __contains__(self, trivial_name: str) -> bool:
		"""Allow `name in coll` syntax."""
		return any(p.trivial_name == trivial_name for p in self.properties)

	def __iter__(self):
		"""Allow `for p in coll` syntax."""
		return iter(self.properties)

	def __len__(self) -> int:
		return len(self.properties)

	def names(self) -> list[str]:
		"""All trivial names."""
		return [p.trivial_name for p in self.properties if p.trivial_name is not None]

	def filter(self, predicate: Callable[[MolecularProperty], bool]) -> 'MolPropsCollection':
		return MolPropsCollection([p for p in self.properties if predicate(p)])

	def without_values(self) -> 'MolPropsCollection':
		"""Properties still awaiting data — useful for finding what's missing."""
		return self.filter(lambda p: p.vals is None)

	def by_calc_setup(self, origin: DataOriginInfo) -> 'MolPropsCollection':
		"""All properties computed with a given setup."""
		return self.filter(lambda p: p.calc_setup == origin)

	def of_order(self, order: int) -> 'MolPropsCollection':
		"""Properties of a specific differentiation order, e.g. 1 for dipole, 2 for polarizability."""
		return self.filter(lambda p: p.order == order)

	def group_by_calc_setup(self) -> dict[DataOriginInfo, 'MolPropsCollection']:
		"""Bucket properties by which setup they use. For batching QC jobs."""
		from collections import defaultdict
		groups = defaultdict(list)
		for p in self.properties:
			groups[p.calc_setup].append(p)
		return {k: MolPropsCollection(v) for k, v in groups.items()}

	def dress(self, uniform: DataOriginInfo | None = None, 
			by_name: dict[str, DataOriginInfo] | None = None):
		"""Attach DataOriginInfo to each property. 
		by_name takes precedence; uniform is the fallback."""
		if uniform is None and by_name is None:
			raise ValueError("Provide `uniform` or `by_name` (or both).")
		
		for p in self.properties:
			if by_name and p.trivial_name in by_name:
				p.addCalcSetup(by_name[p.trivial_name])
			elif uniform is not None:
				p.addCalcSetup(uniform)
			else:
				raise ValueError(f"No setup for property {p}")

	@property
	def are_dressed(self) -> bool:
		return all(isinstance(p.calc_setup, DataOriginInfo) for p in self.properties)

	def build_request_dict(self) -> dict[str, DataOriginInfo]:
		"""Build a {name: DataOriginInfo} shopping list."""
		if not self.are_dressed:
			raise RuntimeError("Collection must be dressed before requesting data.")
		return {p.trivial_name: p.calc_setup for p in self.properties}

	def fill_from(self, data_dict: dict):
		"""Load obtained data into each property's .vals."""
		for p in self.properties:
			if p.trivial_name in data_dict:
				p.vals = data_dict[p.trivial_name]

	@property
	def is_filled(self) -> bool:
		return all(p.vals is not None for p in self.properties)


# numerical
@dataclass
class VibState:
	"""
	Class to represent a vibrational state.
	This is for a "concrete" vibrational state and not the same as its symbolic namesake in wilson-derive.

	----
	s: dictionary {(harm. quanta): coeff, (harm. quanta): coeff, ...}: Specify the state in terms of harm. osc. WFs
	e: float: State energy level
	d: type not specified: Should be some form of vector to represent displacement in terms of atomic coordinates

	UPD:
	dictionary self.s is not JSON-serializable (tuples can't be keys), but self.serial_s is.
	self.serial_s is set up in post_init; deserialize_state_dict will return original self.s based on self.serial_s.

	Notes:
	s: InitVar[dict] = field(repr=False) - means that this atribute will not be in repr() of the class instance
	InitVar - is an init-only variable
	This seems to be okay for now, but should mind this feature
	"""
	harm_quanta_coeffs: dict[tuple[int, ...], float]
	energy: float = 0.0
	displacement: Any = None
	serial_harm_quanta_coeffs: dict[str, float] = field(init=False)
	state_label: str = None
	harmonic_WF: bool = None

	def __post_init__(self) -> None:
		"""Convert tuple keys to comma-separated strings for JSON serialization."""
		self.serial_harm_quanta_coeffs = {
			",".join(map(str, k)): v
			for k, v in self.harm_quanta_coeffs.items()
		}

	def deserialize_state_dict(self) -> dict[tuple[int, ...], float]:
		"""Convert serialized dictionary back to original format with tuple keys."""
		return {
			tuple(int(x) for x in k.split(",")): v
			for k, v in self.serial_harm_quanta_coeffs.items()
		}

	def __eq__(self, other: 'VibState') -> bool:
		if not isinstance(other, VibState):
			return NotImplemented
		return self.state_label == other.state_label and np.isclose(self.energy, other.energy)

	def __lt__(self, other: 'VibState') -> bool:
		if not isinstance(other, VibState):
			return NotImplemented
		return self.state_label < other.state_label
	
	@classmethod
	def get_1q_states(cls, states: list['VibState']):
		return [s for s in states if len(s.state_label.split(','))==1]


# numerical
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


def make_vibdiff_key(vibdiff_term: VibDiffTerm, index_dict: dict) -> tuple[str, str]:
    """
    Non-sorted key for VibDiffBank_cache

    returns keys for vibdiff bank for vib states expression and choice of indices
    """
    left_state_symb = vibdiff_term.sl.q
    right_state_symb = vibdiff_term.sr.q


    left_state_label = ','.join([str(i) for i in sorted([index_dict[i] for i in left_state_symb])])
    right_state_label = ','.join([str(i) for i in sorted([index_dict[i] for i in right_state_symb])])
    
    if left_state_label == '':
        left_state_label = 'zero'
    if right_state_label == '':
        right_state_label = 'zero'
    
    return (left_state_label, right_state_label)


# numerical
@dataclass
class VibDiff:
    """
    Represents difference between two vibrational states.
    Numerical representation that holds values, as opposed to VibDiffTerm which is symbolic.
    Handles special case of zero states (ground state) in comparisons.
    """
    left: VibState
    right: VibState
    
    def is_zero_state(self, state: VibState) -> bool:
        """
        Check if state is a zero (ground) state.

        #TODO more criteria?
        """
        return state.state_label == 'zero'
    
    def normalized(self) -> 'VibDiff':
        """
        Return normalized form where left <= right.
        Zero states are considered smaller than any other state.
        """
        left_is_zero = self.is_zero_state(self.left)
        right_is_zero = self.is_zero_state(self.right)
        
        # If both are zero states or neither is zero, use standard comparison
        if left_is_zero == right_is_zero:
            if self.left < self.right:
                return VibDiff(self.left, self.right)
            return VibDiff(self.right, self.left)
            
        # Zero state should always be on the left
        if left_is_zero:
            return VibDiff(self.left, self.right)
        return VibDiff(self.right, self.left)
    
    def energy_difference(self, *, au=False) -> float:
        """
        Calculate energy difference between states.
        For zero states, energy is considered to be 0.0
        """
        left_energy = 0.0 if self.is_zero_state(self.left) else self.left.energy
        right_energy = 0.0 if self.is_zero_state(self.right) else self.right.energy
        if au:
            return convNu2Ene(left_energy - right_energy)
        else:
            return left_energy - right_energy

    @classmethod
    def from_symbolic(cls, 
                    vibdiff_term_symb: VibDiffTerm,
                    index_dict: dict,
                    vibstates_data: 'VibStatesData') -> 'VibDiff':
        """Construct VibDiff from symbolic representation."""
        # Get state labels from symbolic term
        left_label, right_label = make_vibdiff_key(vibdiff_term_symb, index_dict)
        # Look up states in vibstates_data
        left_state = (
            VibState(harm_quanta_coeffs={}, state_label='zero', energy=0.0)
            if left_label == 'zero'
            else vibstates_data.get_state_by_label(left_label)
        )
        
        right_state = (
            VibState(harm_quanta_coeffs={}, state_label='zero', energy=0.0)
            if right_label == 'zero'
            else vibstates_data.get_state_by_label(right_label)
        )

        return cls(left=left_state, right=right_state)

    def cache_it(self, vibdiff_cache: 'VibDiffCache'):
        """Ensure this VibDiff's energy is cached."""
        if vibdiff_cache.get(self) is None:
            energy = self.energy_difference()
            vibdiff_cache.add(self, energy)


# numerical
@dataclass
class VibDiffCache:
    """
    bank keys
    ('0', '2')   -> sorted version is the key --- ('0', '2')
    ('0,2', '2') -> sorted version is the key --- ('2', '0,2')
    ('1,2', '3') -> sorted version is the key --- ('3', '1,2')
    ('1,2,4', '3,1') -> sorted version is the key --- ('1,3', '1,2,4')
    ('4,1,2', '3,1') -> sorted version is the key --- ('1,3', '1,2,4')

    """
    def __init__(self):
        self._cache: dict[tuple[str, str], float] = {}
    
    def __repr__(self):
        return str(self._cache)
    
    def get(self, vib_diff: VibDiff) -> float | None:
        """Get cached energy difference"""
        key = (vib_diff.left.state_label, vib_diff.right.state_label)
        norm_diff = vib_diff.normalized()
        norm_key = (norm_diff.left.state_label, norm_diff.right.state_label)
        
        if key in self._cache:
            return self._cache[key]
        if norm_key in self._cache:
            return -self._cache[norm_key] if key != norm_key else self._cache[norm_key]
        return None
        
    def add(self, vib_diff: VibDiff, energy: float):
        """Cache energy difference"""
        norm_diff = vib_diff.normalized()
        self._cache[(norm_diff.left.state_label, norm_diff.right.state_label)] = energy

## ------------------------------------------------------------------


@dataclass(frozen=True)
class EvaluationDataAndConfigs:
    """
    data holding abstractions are used here
    """
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


@dataclass
class RspEvalTerm:

    term_id: int | str          # provenance back to the symbolic term

    rot_avrg_props: PropsCollection | None
    rot_invr_props: PropsCollection | None

    res_conds: ResonanceMotif | None

    vibdiffs: FreqTermsCollection | None
    ene_prefac: FreqTermsCollection | None
    num_coeff: FreqTermsCollection | None

    summation_indices: tuple[str, ...] | None  # from tellNonSummSummIndices
    non_summation_indices: tuple[str, ...] | None
    all_indices: tuple[str, ...] | None       # sorted union



def parse_vibpert_term(term: 'VibPerturbedTerm') -> RspEvalTerm:
    
    # extract AVRG and NON_AVRG expressions
    avrg_expr = PropsCollection(props=term.props).get_averaged_props().sort()
    non_avrg_expr = PropsCollection(props=term.props).get_non_averaged_props()

    res_conds = ResonanceMotif(term.res)
    
    # extract frequency terms and their differences
    freqterms_all = FreqTermsCollection(freqterms=term.freqterms)
    extra_freqterms = freqterms_all.get_pert_wf_diff() # fixme

    ene_prefac = freqterms_all.get_vibenedenom() # fixme
    """
    precalculated_data.vibenedenoms_tensors[freqterms.get_num_indices_vibenedenom()]
    """

    num_coeff = float(term.coeff)

    idx_summ, idx_nonsumm = term.tellNonSummSummIndices()

    return RspEvalTerm(rot_avrg_props=avrg_expr,
                       rot_invr_props=non_avrg_expr,
                       res_conds=res_conds,
                       vibdiffs=extra_freqterms,
                       ene_prefac=ene_prefac,
                       num_coeff=num_coeff,
                       summation_indices=idx_summ,
                       non_summation_indices=idx_nonsumm)


def evaluate_term_coeffs(term: 'VibPerturbedTerm', 
                         relevant_indices: list[dict], 
                         necessary_data: tuple['EvaluationDataAndConfigs', 'PrecalculatedData'], 
                         zero_tol: float = 1e-18) -> dict['ParameterSet', float]:
    """
    Evaluate the coefficient part of the term 'term' for all of the indices in 'relevant_indices',
    handling hierarchical summation over indices.

    example_relevant_indices = [
        {'e': 0},  # Sum over 'a', 'b', 'c', 'd'
        {'e': 1},  # Sum over 'a', 'b', 'c', 'd'
        {'e': 2},  # Sum over 'a', 'b', 'c', 'd'
    ]

    Parameters:
        term: VibPerturbedTerm
            The term to evaluate, containing properties and frequency terms.
        relevant_indices: List[Dict]
            List of dictionaries specifying the relevant indices for evaluation.
        necessary_data: Tuple[EvaluationDataAndConfigs, PrecalculatedData]
            Tuple containing data and configurations required for evaluation.
        zero_tol: float
            Tolerance for considering a value as zero.
    Returns:
        Dict[ParameterSet, float]: A dictionary mapping ParameterSet to computed coefficients.
    """
    results = {}
    data_and_configs, precalculated_data = necessary_data
    
    # extract AVRG and NON_AVRG expressions
    avrg_expr = PropsCollection(props=term.props).get_averaged_props().sort()
    non_avrg_expr = PropsCollection(props=term.props).get_non_averaged_props()
    
    # extract frequency terms and their differences
    freqterms = FreqTermsCollection(freqterms=term.freqterms)
    extra_freqterms = freqterms.get_pert_wf_diff()
    
    # Get all indices
    idx_summ, idx_nonsumm = term.tellNonSummSummIndices()
    term_idx_all = sorted(idx_summ + idx_nonsumm)
    
    def hierarchical_sum(index_dict: dict, remaining_indices: list[str], dict_of_sum: dict) -> float:
        """
        Perform hierarchical summation over the remaining indices.
        Parameters:
            index_dict: Dict
                The current index dictionary with some indices fixed.
            remaining_indices: List[str]
                The list of indices that still need to be summed over.
        Returns:
            float: The result of the summation for the given index dictionary.
            
            dict_of_sum: top level: 
                    {ParameterSet(index_dict): {}}
        """
        # Base case: no remaining indices to sum over
        if not remaining_indices:
            value, contribs = evaluate_single_index_dict(term, index_dict, avrg_expr, non_avrg_expr, extra_freqterms, freqterms, data_and_configs, precalculated_data, zero_tol)
            dict_of_sum[ParameterSet(index_dict)] = contribs
            # returns coef and dict with param contribs

            return value, dict_of_sum
        
        # Get the next index to sum over
        current_index = remaining_indices[0]
        remaining = remaining_indices[1:]
        
        # Perform summation over all possible values for the current index
        n_modes = data_and_configs.number_of_nmodes
        total_sum = 0.0

        for value in range(n_modes):
            # Update the index dictionary with the current value
            new_index_dict = index_dict.copy()
            new_index_dict[current_index] = value
            
            parent_key = ParameterSet(index_dict)
            child_key = ParameterSet(new_index_dict)
            # if parent_key not in dict_of_sum:
            #     dict_of_sum[parent_key] = {}
            # dict_of_sum[parent_key][child_key] = {}
            
            # Recursively compute the sum for the remaining indices
            recurse_res = hierarchical_sum(new_index_dict, remaining, dict_of_sum[ParameterSet(index_dict)])
            total_sum += recurse_res[0]
        
        return total_sum, dict_of_sum
    
    # Iterate over the relevant indices -- 
    for index_dict in relevant_indices:
        # Identify missing indices
        missing_indices = [index for index in term_idx_all if index not in index_dict]
        
        dict_of_sum = {ParameterSet(index_dict): {}}
        
        # Perform hierarchical summation for the current index_dict
        result = hierarchical_sum(index_dict, missing_indices, dict_of_sum)
        results[ParameterSet(index_dict)] = result
    
    return results


def evaluate_single_index_dict(term: 'VibPerturbedTerm', 
                               index_dict: dict, 
                               avrg_expr, 
                               non_avrg_expr, 
                               extra_freqterms, 
                               freqterms, 
                               data_and_configs, 
                               precalculated_data, 
                               zero_tol: float) -> tuple[float, dict]:
    """
    Evaluate the term for a single index dictionary.
    Parameters:
        (same as evaluate_term_coeffs)
    
    Returns:
        float: The computed coefficient for the given index dictionary.
    """
    # Evaluate NON_AVRG
    NON_AVRG = eval_non_avrg_per_indexdict(non_avrg_expr, index_dict, data_and_configs, zero_tol)
    if NON_AVRG == 0.0:
        # print('\nNON_AVRG zero - ', non_avrg_expr, index_dict, '\n\n')
        AVRG = eval_avrg_per_indexdict(avrg_expr, index_dict, precalculated_data, zero_tol)
        return 0.0, {'NON_AVRG': NON_AVRG, 'AVRG': AVRG}
    # Evaluate AVRG
    AVRG = eval_avrg_per_indexdict(avrg_expr, index_dict, precalculated_data, zero_tol)
    if AVRG == 0.0:
        # print('\nAVRG zero - ', avrg_expr, index_dict, '\n\n')
        return 0.0, {'AVRG': AVRG}
    # Evaluate VIBDIFF_TERMS
    VIBDIFF_TERMS = eval_vibdiff_pert_wf_diff(extra_freqterms, index_dict, precalculated_data, data_and_configs)
    if VIBDIFF_TERMS == 0.0:
        # print('\nVIBDIFF_TERMS zero - ', extra_freqterms, index_dict, '\n\n')
        return 0.0, {'VIBDIFF_TERMS': VIBDIFF_TERMS}
    # Evaluate VIBENE_DENOM
    VIBENE_DENOM = eval_vibenedenom(freqterms, index_dict, precalculated_data)
    if VIBENE_DENOM == 0.0:
        # print('\nVIBENE_DENOM zero - ', freqterms, index_dict, '\n\n')
        return 0.0, {'VIBENE_DENOM': VIBENE_DENOM}
    # Compute the product
    product_all = NON_AVRG * AVRG * VIBDIFF_TERMS * VIBENE_DENOM

    # print('\nindex_dict', index_dict)
    # print('NON_AVRG', NON_AVRG)
    # print('AVRG', AVRG)
    # print('VIBDIFF_TERMS', VIBDIFF_TERMS)
    # print('VIBENE_DENOM', VIBENE_DENOM, '\n')

    dict_contribs = {'NON_AVRG': NON_AVRG, 'AVRG': AVRG, 'VIBDIFF_TERMS': VIBDIFF_TERMS, 'VIBENE_DENOM': VIBENE_DENOM}

    return float(term.coeff) * float(product_all), dict_contribs


# TODO: error handling for missing or invalid data for all functions below
def eval_non_avrg_per_indexdict(non_avrg_expr: PropsCollection, 
                                index_dict: dict, 
                                data_and_configs: EvaluationDataAndConfigs, 
                                zero_tol: float = 1e-18):
    """
    non_avrg_expr - extracted part of VibPerturbed term 

    order of indices generally: a,b,c,... 
    E.g. in CFF tensor index tuple is (a,b,c)
    """
    product_all = 1.
    
    for non_avrg_prop in non_avrg_expr:
        # accessing values for non-averaged properties from data
        na_prop_inds = tuple([index_dict[i] for i in non_avrg_prop.inds])
        triv_name = prop_trivname(ord_el=len(non_avrg_prop.ops), ord_geo=non_avrg_prop.dord)

        NON_AVRG = data_and_configs.props_data.get(triv_name).vals[na_prop_inds]

        if np.isclose(NON_AVRG, zero_tol):
            return 0.
        else:
            product_all *= NON_AVRG
    
    return product_all


def get_ind_tuple_from_base(expr: PropsCollection, base_expr: PropsCollection, index_dict: dict):
    """Map expr to indices according to base expression's unique symbols."""
    base_unique = sorted(list(set(base_expr.get_mode_indices())))
    expr_inds = expr.get_mode_indices()

    if len(base_unique) < len(expr_inds):
        # walk through only base_unique ind labels
        return tuple(index_dict[sym] for sym in base_unique)
    elif len(base_unique) == len(expr_inds):
        # walk through all expr_inds labels, there are repeated labels
        return tuple(index_dict[sym] for sym in expr_inds)
    else:
        # this should not be possible in the worflow
        raise ValueError('This base_expr cannot be a base expression for this expr')


def eval_avrg_per_indexdict(avrg_expr: PropsCollection, 
                            index_dict: dict, 
                            precalculated_data: PrecalculatedData,
                            zero_tol: float = 1e-18):
    avrg_tensor_expr = precalculated_data.avrg_expr_tensor_mapping[avrg_expr]
    
    avrg_tensor = precalculated_data.avrg_tensors[avrg_tensor_expr]
    
    avrg_index_tuple = get_ind_tuple_from_base(expr=avrg_expr, 
                                                         base_expr=avrg_tensor_expr, 
                                                         index_dict=index_dict)
    if np.isclose(avrg_tensor[avrg_index_tuple], zero_tol):
        return 0.
    return avrg_tensor[avrg_index_tuple]

def eval_vibdiff_pert_wf_diff(extra_freqterms: FreqTermsCollection,
                              index_dict: dict,
                              precalculated_data: PrecalculatedData,
                              data_and_configs: EvaluationDataAndConfigs):
    product_all = 1.

    for vibdiff in extra_freqterms:
        vib_diff_w_value = VibDiff.from_symbolic(vibdiff, index_dict, 
                                                        data_and_configs.vibstates_data)
        vib_diff_w_value.cache_it(vibdiff_cache=precalculated_data.vibdiff_cache)

        product_all *= 1./ vib_diff_w_value.energy_difference(au=True)
    
    return product_all

def eval_vibenedenom(freqterms: FreqTermsCollection,
                     index_dict: dict,
                     precalculated_data: PrecalculatedData):
    
    vibenedenoms_tensor = precalculated_data.vibenedenoms_tensors[freqterms.get_num_indices_vibenedenom()]

    vibeneden_index_tuple = tuple([index_dict[i] for i in freqterms.get_num_indices_vibenedenom()])

    return vibenedenoms_tensor[vibeneden_index_tuple]

