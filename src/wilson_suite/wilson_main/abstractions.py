from dataclasses import dataclass, field, asdict, is_dataclass, InitVar
from typing import Callable, Any, Optional


import logging

import numpy as np
logger = logging.getLogger("wilson")

# A system is here only the system name, molecular geometry and atoms (masses for isotopes?)
@dataclass
class MolecularSystem:
	"""
	name: String: System name

	geo: Currently not fully fixed form but must specify system geometry in a form
	[[atom name, [x, y, z](Ångström)], ...], where atom_name must be and element name and where
	default isotopes are assumed

	natoms: integer: Number of atoms

	geo_extra: Currently not fixed form but reserved for situations where extra information is needed
	Speculated example format: # [[atom name, [x, y, z] (Ångström), (# protons, # neutrons)], ...]
	and possibly more attributes if needed
	"""
	name: str
	natoms: int = None
	geo: Any = None
	geo_extra: Any = None
	linear: bool = False
	conformer: str = 'conf1'

	@property
	def Nnmodes(self):
		if self.linear:
			return 3*self.natoms-5
		else:
			return 3*self.natoms-6

	def __post_init__(self):
		if (self.geo is not None) and (self.geo_extra is not None):
			raise AssertionError('Ambigious definition: Both default form geometry geo and extended form geometry geo_extra was defined')
		if (self.geo is None) and (self.geo_extra is None):
			logger.info('Note: Molecular system was instantiated without geometry information')
		if self.natoms is None and self.geo is None:
			raise ValueError('Incomplete definition: Either the number of atoms (natoms) or the geometry (geo) of the MolecularSystem is required')
		
		if self.natoms < 3:
			self.linear = True

	def __hash__(self):
		# FIXME: Influence hash with other attributes as well
		return hash(self.name)



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
	trivial_name: str=None
	vals: InitVar[Any] = field(default=None, repr=False)
	calc_setup: DataOriginInfo = None
	extra_data: dict = None

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

	def addSystem(self, system: MolecularSystem):
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


# FIXMEs: Improved handling of mode exclusion; possibly methods changes
@dataclass
class VibAnaSetup:
	"""
    Class for setup for vibrational analysis and storage of the resulting information
	
	----
	regime: string: Vibrational analysis regime (e.g. "harmonic", "GVPT2", "VPT2")
	system: MolecularSystem instance: System to which this instance pertains
	regime_subinfo: dictionary: Extra configuration info for vibrational regime (e.g. skip rotational effects)
	max_state_lvl: integer: Maximum number of vibrational quanta per harmonic state involved in states
	states: List of VibState instances: Specification of each vibrational state in scope
	nc_sqrt_eigval: dictionary {mode index: value}: Harmonic vibrational energy levels
	nc_eigvec: dictionary {mode index: [values]}: Normal mode displacements (canonically in Cartesian basis)
	vibana_own_analysis: String: Which kinds of properties will I need to actually carry out the vibrational analysis? 
		Choices: 
			"full": I need properties for both harmonic and (if chosen) anharmonic analysis,
			"anharm": I only need properties to carry out an anharmonic analysis [I will or have already gotten the harmonic data], 
			"none": I don't need any properties [I will or have already gotten both harmonic and anharmonic data]
	exclude_modes: list: Tells which modes (if any) to exclude in this vibrational analysis
	"""
	regime: str=None
	system: MolecularSystem=None
	regime_subinfo: dict=None

	# 'full': Will need properties for harmonic (and, if requested, anharmonic) analysis
	# 'anharm': Will only need props. for anharmonic analysis (harmonic results will be provided from external source)
	# 'none': All results will be provided by external source
	vibana_own_analysis: str='full'

	max_state_lvl: int=None
	states: list[VibState]=field(default=None, repr=False)
	
	number_of_modes: int = None
	
	# Dictionary: {nm index: w}
	nc_sqrt_eigval: dict=None
	nc_eigvec: dict=None

	# TODO: MODE EXCLUSION, REGISTERING OF FERMI RESONANCES (TO BE PASSED TO EVALUATOR)
	exclude_modes: list = None
	diagn: dict = None

	def __post_init__(self):
		if self.number_of_modes is None:
			if self.system is not None:
				self.number_of_modes = self.system.Nnmodes
		
		if self.exclude_modes is None:
			if self.system is not None:
				self.exclude_modes = []
		else:
			if self.system is None:
				logger.warning('VibAnaSetup().exclude_modes attribute is not meaningfull without having set the VibAnaSetup().system attribute')

	@property
	def modes_indices(self):
		"""
		Automatically set up based on number of modes in the system (3*N-5 or 3*N-6) and exclude_modes list, 
			which is an empty list by default
		Returns empty list if no system attribute
		"""
		if self.system is None:
			raise AttributeError('VibAnaSetup().modes_indices attribute cannot be created without having set the VibAnaSetup().system attribute')

		import numpy as np
		return [int(i) for i in np.arange(self.system.Nnmodes) if i not in self.exclude_modes]
	
	# FIXME: Possibly return to this later
	@property
	def serial_states(self):
		return [{'s': vibst.serial_s, 
		   'e': vibst.e, 'd': vibst.d} for vibst in getattr(self, 'states')]

	@property
	def has_all_states(self):
		"""
		#FIXME: now state level (number of quanta) is obtained from the label!
		"""
		if self.states is None:
			return False
		states_lvls = [len(next(iter(st.harm_quanta_coeffs))) for st in self.states] # takes first key in harm_quanta_coeffs dict
		if self.max_state_lvl in states_lvls:
			return True
		print('states_lvls', states_lvls)
		print('self.max_state_lvl', self.max_state_lvl)
		return False

	@property
	def isAllSet(self):
		"""
		Checking status of VibAna data.
		"""
		if self.nc_sqrt_eigval is not None and self.has_all_states:
			return True
		print('\nself.nc_sqrt_eigval is not None:', self.nc_sqrt_eigval is not None)
		print('self.has_all_states:', self.has_all_states, '\n')
		return False
	
	def setStates(self, states: list[VibState]):
		"""
		Set vibrational states

		states: List of VibState instances: The states to be set
		"""

		self.states = states

	def upd_exclude_modes(self, upd_exclude_modes: list = None):
		"overwrites self.exclude_modes"
		
		if self.exclude_modes is None:
			if self.system is not None:
				self.exclude_modes = []
		elif upd_exclude_modes is not None:
			self.exclude_modes = upd_exclude_modes
		else:
			if self.system is None:
				logger.warning('VibAnaSetup().exclude_modes attribute is not meaningfull without having set the VibAnaSetup().system attribute')


	def set_include_modes_list(self):
		if self.nc_sqrt_eigval is None:
			raise ValueError("Check nc_sqrt_eigval before setting up a modes inclusion list (nc_sqrt_eigval needs to be set).")

		self.include_list = tuple([v for v in list(self.nc_sqrt_eigval.keys()) if v not in self.exclude_modes])
		if self.include_list == tuple():
			raise ValueError("include_list of included normal modes labels is empty")
