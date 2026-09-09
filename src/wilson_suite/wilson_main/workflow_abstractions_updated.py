"""
1. VibExperiment object  -> 
2. derive Terms -> 
3. configure necessary settings (eval/render) -> 
4. request Data for evaluation [necessary for intensities eval, vib info - max state level only, nothing about vib analysis] -> 
5. get DataResults

---
terms [props, max_state_lvl] -- from get_fully_enhanced_terms(VibExperiment, ...)
magn_conds - from VibExperiment
pulse_polarization_vector -- terms should have that as a parameter???
gamma
axis_choice

vib_ana_configs (harm/ahnarm - anharm_inhouse or not)

dyn_range
spec_window

data:
number of normal modes (total and choices here)
states and energies (all)


get_fully_enhanced_terms(VibExperiment, ...)


"""
import copy
from pathlib import Path
from typing import Any, Callable
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
	from .spectrum_abstractions import SpecEvalSetup
	from wilson_suite.wilson_analysis.render.spectrum_renderer import SpectrumRenderer
	from .abstractions import VibAnaSetup

from .main_functions import find_props_and_max_state_lvl
from .abstractions import MolPropsCollection
from .abstractions import MolecularSystem, DataOriginInfo
from ..wilson_derive.response_terms import VibPerturbedTerm, VibPertTermsCollection
from wilson_suite.wilson_experiment.experiment_abstractions import VibExperiment
from ..wilson_experiment.indep_vars_and_axes import SpectralAxisSet
from wilson_suite.wilson_utils.wilson_data_obtainer import wilson_data_obtainer
from wilson_suite.wilson_analysis.render.matplotlib_renderer import MatplotlibRendererGV
from wilson_suite.wilson_derive.term_var_translate import translate_magn_conditions_to_axisvars
import numpy as np

import logging
logger = logging.getLogger("wilson")


@dataclass
class SimulationSetup:
    experiment: VibExperiment
    system: MolecularSystem
    terms: list[VibPerturbedTerm]
    vib_ana: 'VibAnaSetup'
    spec_eval: 'SpecEvalSetup'
    axis_choice: SpectralAxisSet | None = None
    eval_uniform: DataOriginInfo | None = None
    eval_by_prop: dict[str, DataOriginInfo] | None = None

# derived from setup + terms
@dataclass
class PropertyRequest:
	props_coll: MolPropsCollection
	residual_vib_info: dict
	max_state_lvl: int

	@property
	def is_dressed(self):
		return self.props_coll.are_dressed


@dataclass
class SimulationResult:
    spec: np.ndarray | None = None
    grid: dict | None = None
    rendering: Any = None
    diagnostics: dict = field(default_factory=dict)

class RunDirectory:
	def __init__(self, base_dir: Path):
		self.base_dir = base_dir

	def save_pickle(self, obj, name: str):
		pass

	def save_json(self, obj, name: str):
		pass
	
	@property
	def figures(self) -> Path:
		pass


@dataclass(frozen=True)
class SealedSetup:
    """Phase A output, Phase B input. Immutable by construction."""
    experiment: VibExperiment
    system: MolecularSystem
    terms_in_axes: VibPertTermsCollection # list[VibPerturbedTerm]
    axis_choice: SpectralAxisSet
    vib_ana: 'VibAnaSetup'
    spec_eval: 'SpecEvalSetup'
    prop_order: PropertyRequest

    @property
    def is_dressed(self) -> bool:
        return self.prop_order.is_dressed and self.vib_ana.isAllSet


class SimulationBuilder:
	"""
	Part 1 - prepare setup and configs.

	The builder takes ownership of all inputs via deep copy. 
	Modifications to the caller's objects after setting will not affect the builder.

	experiment - defines derived terms and possible axes choices


	"""
	def __init__(self, 
			  		experiment: VibExperiment = None, 
					system: MolecularSystem = None,
					vib_ana: 'VibAnaSetup' = None,
					spec_eval: 'SpecEvalSetup' = None,
					eval_uniform: DataOriginInfo = None,
					eval_by_prop: dict[str, DataOriginInfo] = None):
		
		self._experiment = copy.deepcopy(experiment)

		self._explicit_axes = None

		self._explicit_terms: dict | None = None   # user-supplied via set_terms
		self._cached_derived_terms: dict | None = None   # cache of self._experiment.derive_terms()

		self._system = copy.deepcopy(system)
		self._vib_ana = copy.deepcopy(vib_ana)
		self._spec_eval = copy.deepcopy(spec_eval)
		self._eval_uniform = copy.deepcopy(eval_uniform)
		self._eval_by_prop = copy.deepcopy(eval_by_prop)

		self._prop_order: PropertyRequest | None = None
		self._dressed: bool = False

	@property
	def terms(self) -> dict:
		if self._explicit_terms is not None:
			return self._explicit_terms
		if self._cached_derived_terms is None:
			if self._experiment is None:
				raise RuntimeError(
					"Terms not set and no experiment to derive from. "
					"Provide terms directly via set_terms() or set an experiment."
				)
			self._cached_derived_terms = self._experiment.derive_terms()
		return self._cached_derived_terms

	@property
	def prop_order(self) -> PropertyRequest:
		"""
		Resolved property order. Cached after first access.
		freqs='static' only now
		"""
		if self._prop_order is None:
			self._require('_vib_ana')
			props, resid, max_lvl = find_props_and_max_state_lvl(
				self.terms, self._vib_ana, freqs='static'
			)
			self._prop_order = PropertyRequest(
				props_coll=MolPropsCollection(props), 
				residual_vib_info=resid, 
				max_state_lvl=max_lvl
			)
		return self._prop_order

	@property
	def axis_choice(self) -> SpectralAxisSet:
		"""Resolved axis choice: explicit > spec_eval.ev_info > experiment.canonical_axes."""
		if self._explicit_axes is not None:
			return self._explicit_axes
		if self._spec_eval and self._spec_eval.ev_info.spectral_axes is not None:
			return self._spec_eval.ev_info.spectral_axes
		if self._experiment is not None:
			return self._experiment.canonical_axes
		raise RuntimeError("No axis_choice set and no defaults available")

	# --- setters ---
	def set_experiment(self, experiment: VibExperiment): 
		"""
		Setting an experiment will remove terms, prop_order if previously set, 
		also _explicit_axes choice could not be appropriate for new experiment
		"""
		self._experiment = copy.deepcopy(experiment)

		self._explicit_axes = None
		self._cached_derived_terms = None
		self._prop_order = None
		self._dressed = False

	def set_system(self, system: MolecularSystem):
		# FIXME - how cheap is it to copy MolecularSystem
		# FIXME - upd so system initiates with natoms and linear or not - for Nnmodes
		self._system = copy.deepcopy(system)

	def set_terms(self, terms: dict[int, dict[tuple, VibPerturbedTerm]]):
		self._explicit_terms = copy.deepcopy(terms)

	def set_axis_choice(self, axes: SpectralAxisSet):
		"""
		updated axes choice for builder but not in spec eval
		"""
		# TODO: validate axes against experiment - as what are possible axes
		self._explicit_axes = copy.deepcopy(axes)

	def set_vib_ana(self, vib_ana: 'VibAnaSetup'):
		"""
		making a copy of va.

		Setting a VibAnaSetup will remove prop_order if previously set
		"""
		# FIXME - how cheap is it to copy VibAnaSetup
		self._vib_ana = copy.deepcopy(vib_ana)

		self._prop_order = None
		self._dressed = False

	def set_spec_eval(self, spec_eval: 'SpecEvalSetup'):
		"""
		making a copy of se
		"""
		self._spec_eval = copy.deepcopy(spec_eval)

	def set_eval_uniform(self, origin: DataOriginInfo): 
		"""
		self will be set as undressed, assuming this is a new DataOriginInfo
		"""
		self._eval_uniform = copy.deepcopy(origin)
		self._dressed = False

	def set_eval_by_prop(self, mapping: dict[str, DataOriginInfo]):
		"""
		self will be set as undressed, assuming this is a new DataOriginInfo mapping
		"""
		self._eval_by_prop = copy.deepcopy(mapping)
		self._dressed = False


	@property
	def _terms_collection(self) -> VibPertTermsCollection:
		"""Fresh wrapper around terms; built on each access."""
		return VibPertTermsCollection(
			term_dict=self.terms, experiment=self._experiment
		)

	def _dress_residual(self, uniform, by_name):
		for key in self.prop_order.residual_vib_info:
			if by_name and key in by_name:
				self.prop_order.residual_vib_info[key] = by_name[key]
			elif uniform is not None:
				self.prop_order.residual_vib_info[key] = uniform
			else:
				raise ValueError(f"No setup for residual vib info item {key!r}")
	
	def _dress_props(self, uniform, by_name):
		self.prop_order.props_coll.dress(
			uniform=uniform, 
			by_name=by_name
		)

	def dress_prop_order(self):
		"""
		TODO: check this
		"""
		if self._eval_uniform is None and self._eval_by_prop is None:
			raise RuntimeError("Provide eval_uniform or eval_by_prop first.")

		self._dress_props(
			uniform=self._eval_uniform, 
			by_name=self._eval_by_prop
		)

		self._dress_residual(
			uniform=self._eval_uniform, 
			by_name=self._eval_by_prop
		)
		
		self._dressed = True

	def fill_defaults(self):
		"""Explicit soft defaults. Call before seal()"""
		# TODO: spectral window, resolution, damping defaults go here
		pass


	@property
	def missing(self) -> list[str]:
		required = ['_experiment', '_system',
					'_vib_ana', '_spec_eval']
		return [name[1:] for name in required if getattr(self, name) is None]

	def _build_sealed(self):

		axes = self.axis_choice
		
		experiment=copy.deepcopy(self._experiment)
    	
		magn_conditions_translated = None
		if experiment is not None and experiment.magn_conditions is not None:
			magn_conditions_translated = translate_magn_conditions_to_axisvars(
				experiment.magn_conditions, axes
			)

		ev_info = self._spec_eval.ev_info if self._spec_eval is not None else None
		if ev_info is not None:
			# legacy field: read only by EvaluationWorkflow. Normally identical to what we
			# just translated (same source, same function) - only flag a real mismatch.
			if (ev_info.exp_magn_conditions is not None
					and ev_info.exp_magn_conditions != magn_conditions_translated):
				logger.warning(
					"spec_eval.ev_info.exp_magn_conditions (%r) disagrees with the conditions "
					"translated from experiment.magn_conditions (%r). The sealed-setup path uses "
					"the latter; the ev_info field is read only by the legacy EvaluationWorkflow. "
					"This usually means it was set by hand, or translated under a different "
					"axis_choice than the current one.",
					ev_info.exp_magn_conditions, magn_conditions_translated,
				)
			if ev_info.apply_exp_magn_conditions_eval and magn_conditions_translated is None:
				raise RuntimeError(
					"apply_exp_magn_conditions_eval is True but no magnitude conditions are "
					"available: experiment.magn_conditions is None. Set them on the experiment, "
					"or switch the eval filter off."
				)


		terms_in_axes = self._terms_collection.translate_to_ax_choice(axes, magn_conditions_translated)
		
		return SealedSetup(
			experiment=experiment,
			system=copy.deepcopy(self._system),
			terms_in_axes=terms_in_axes,
			axis_choice=copy.deepcopy(axes),
			vib_ana=copy.deepcopy(self._vib_ana),
			spec_eval=copy.deepcopy(self._spec_eval),
			prop_order=copy.deepcopy(self.prop_order),
		)

	def seal(self) -> SealedSetup:
		"""
		returns full setup ready to be passed to evaluation step.
		"""
		if self.missing:
			raise RuntimeError(f"Cannot seal: missing {self.missing}")		
		
		if not self._dressed:
			raise RuntimeError("Cannot seal: call dress_prop_order() first")
		
		return self._build_sealed()

	@property
	def missing_for_dry_run(self) -> list[str]:
		"""What's needed to produce a dry-run seal."""
		required = ['_experiment', '_vib_ana']  # subset
		return [n[1:] for n in required if getattr(self, n) is None]

	def seal_dry_run(self) -> SealedSetup:
		"""Seal without dressing properties — useful for inspecting the shopping list
		before committing to calc setups. The resulting SealedSetup has is_dressed=False
		and cannot be executed by SimulationRun."""
		if self.missing_for_dry_run:
			raise RuntimeError(f"Cannot dry-run seal: missing {self.missing_for_dry_run}")
		
		return self._build_sealed()


	def _require(self, *attrs):
		absent = [a for a in attrs if getattr(self, a) is None]
		if absent:
			raise RuntimeError(f"Need to set: {absent}")


@dataclass
class EvaluatedResult:
    spec: np.ndarray
    grid: dict
    artifacts: dict


def default_renderer_factory(result: EvaluatedResult, 
							 sealed: SealedSetup) -> 'SpectrumRenderer':
	"""Default renderer: matplotlib, using whatever rnd_info the sealed setup specifies."""
    
	pass
	backend = sealed.spec_eval.rnd_info.backend
	if backend == 'matplotlib':
		return MatplotlibRendererGV(result, sealed)
	else:
		raise NotImplementedError(f"Backend {backend!r} not supported")


def default_evaluator(sealed: SealedSetup, run_dir: RunDirectory,
			 save_evalinputs_pkl: str = None, verbose: bool = False) -> EvaluatedResult:
	"""
	Evaluating method, using EvaluationWorkflow
	"""
	from wilson_suite.wilson_intensities.amplitudes.evaluation_wf import EvaluationWorkflow_NEW
	
	# save EvaluationInputs data optionally to a pickle file
	if save_evalinputs_pkl is not None:
		
		from wilson_suite.wilson_utils.serialization import pickle_this_to
		pickle_this_to(sealed, filenamepkl='EvaluationInputs.pkl', save_to=run_dir)
	
	# sealed.spec_eval doesn't need to have rnd_info
	workflow = EvaluationWorkflow_NEW(setup_inputs=sealed)
	wf_result = workflow.run(verbose=verbose)
	
	result = EvaluatedResult(spec=wf_result['result'], 
						  	 grid={k:v for k,v in wf_result.items() if k!='result'},
							 artifacts=workflow.artifacts)
	return result

class SimulationRun:
	"""Phase B: execute a sealed setup."""
	def __init__(self, 
					sealed: SealedSetup, 
					obtainer: Callable = wilson_data_obtainer, 
					evaluator: Callable = default_evaluator,
					renderer_factory: Callable = default_renderer_factory, 
					run_dir: RunDirectory | None = None):
		"""
		obtainer=wilson_data_obtainer -- dafault obtainer function
		renderer_factory=default_renderer_factory, 

		"""
		self._sealed = sealed
		self._obtainer = obtainer
		self._evaluator = evaluator
		self._renderer_factory = renderer_factory
		self._run_dir = run_dir
		self._result: EvaluatedResult | None = None

		self._data_filled = False
	

	def _build_request_dict(self) -> dict:
	
		req = self._sealed.prop_order.props_coll.build_request_dict()
		req.update(self._sealed.prop_order.residual_vib_info)
		return req

	def _fill_results(self, data_dict: dict):
		"""
		loading data into self.props (and optionally to self.vib_ana_setup)
		
		data_dict: dict - {data_name: values}

		"""
		from .main_functions import fill_residual_vib_info_results
		self._sealed.prop_order.props_coll.fill_from(data_dict)
		fill_residual_vib_info_results(self._sealed.vib_ana, self._sealed.prop_order.residual_vib_info, data_dict)

	def obtain_data(self, exclude_modes: tuple = tuple()):
		"""Request data via obtainer, fill into prop_order."""
		request = self._build_request_dict() 
		data_dict = self._obtainer(request)
		self._fill_results(data_dict)

		if exclude_modes:
			self._sealed.vib_ana.exclude_modes = exclude_modes
		self._sealed.vib_ana.set_include_modes_list()

		self._data_filled = True

	def evaluate(self):
		if not self._data_filled:
			raise RuntimeError("Data not obtained yet. Call obtain_data() first.")
		self._result = self._evaluator(self._sealed, self._run_dir)


	def get_renderer(self):
		"""For interactive use: get a renderer object and call its methods directly."""
		if self._result is None:
			raise RuntimeError("No result to render. Call evaluate() or run() first.")
		return self._renderer_factory(self._result, self._sealed)

	def render(self, mode='contour', **kwargs):
		renderer = self.get_renderer()
		if not hasattr(renderer, mode):
			available = [m for m in dir(renderer) if not m.startswith('_')]
			raise ValueError(f"Unknown render mode {mode!r}. Available: {available}")
		return getattr(renderer, mode)(**kwargs)

	@property
	def result(self) -> EvaluatedResult:
		"""
		exposing result without possibility to write directly
		"""
		if self._result is None:
			raise RuntimeError("No result yet. Call evaluate() or run() first.")
		return self._result

	def execute(self) -> EvaluatedResult:
		if not self.is_data_filled:
			self.obtain_data()
		self.evaluate()
		return self._result

	@property
	def is_data_filled(self) -> bool:
		return self._data_filled

	@property
	def has_result(self) -> bool:
		return self._result is not None

