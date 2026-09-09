from wilson_suite.wilson_experiment.indep_vars_and_axes import IndependentVariableSet, SignedPulseTuple, SpectralAxis, SpectralAxisSet


def make_IndependentVariableSet(pulse_refs_tuples: list) -> IndependentVariableSet:
    return IndependentVariableSet(var_set=tuple(SignedPulseTuple(pulse_refs=(t,)) for t in pulse_refs_tuples))


def make_SpectralAxisSet(axes_set_dict: dict[str,list[tuple]]) -> SpectralAxisSet:
    """
    axes_set_dict: {'A': [-1], 'B': [-1, 2]} - lists are sufficient here because each element of the list will be in a SignedPulseTuple

    """
    axes = tuple(SpectralAxis(label=label, var_set=make_IndependentVariableSet(vars)) for label, vars in axes_set_dict.items())
    return SpectralAxisSet(axes=axes)