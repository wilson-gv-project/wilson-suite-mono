from wilson_suite.wilson_experiment.indep_vars_and_axes import SpectralAxis, SpectralAxisChoices
from wilson_suite.wilson_derive.response_terms import VibPerturbedTerm

# --------- SpectralAxis and related
def show_valid_axis_combs(valid_axis_combs: list[SpectralAxisChoices], output=False):
    text = []
    
    for i, phasematch_axchoices in enumerate(valid_axis_combs):
        text.append(f'{i} -- axes choices for phase matching condition {phasematch_axchoices.phasematch_cond.pulses.pulse_refs}')
        text.append(f'Independent variables are: {phasematch_axchoices.ind_vars.var_set}')
        
        for i, axchoice in enumerate(phasematch_axchoices.valid_axis_combs):
            text.append(f"\nSpectralAxisSet {i} with axes:")
            for ax in axchoice.axes:
                text.append(f'---- {show_axis(ax, output=True)}')
    if output:
        return '\n'.join(text)
    print('\n'.join(text))

def show_axis(ax: SpectralAxis, output: False):
    text = f"SpectralAxis with label {ax.label} and var_set {ax.var_set.var_set}"
    if output:
        return text
    print(text)

# --------- VibPerturbedTerm and related
def show_term_latex(term: VibPerturbedTerm, part_of_term: str = None, output=False):
    """ 
    part_of_term: 'res', 'coeff', 'props', 'freqterms'
    """
    res_conditions_denom = ''.join([rc.to_latex() for rc in term.res])
    if res_conditions_denom == '':
        res_conditions_str = ''
    else:
        res_conditions_str = rf'\frac{{1}}{{{res_conditions_denom}}}'

    coefficients_str = rf'\frac{{{term.coeff.numerator}}}{{{term.coeff.denominator}}}'
    properties_str = ''.join([p.to_latex() for p in term.props])

    freqterms_denom = ''.join([rf'\omega_{{{vd.to_latex()}}}' for vd in term.freqterms])
    if freqterms_denom == '':
        freqterms_str = ''
    else:
        freqterms_str = rf'\frac{{{1}}}{{{freqterms_denom}}}'

    if part_of_term is not None:
        if part_of_term == 'res':
            return res_conditions_str
        elif part_of_term == 'coeff':
            return coefficients_str
        elif part_of_term == 'props':
            return properties_str
        elif part_of_term == 'freqterms':
            return freqterms_str
    else:
        return coefficients_str + freqterms_str + properties_str + res_conditions_str
    
