from ..builder_functions import make_SpectralAxisSet

from ..some_reprs import show_valid_axis_combs, show_term_latex
import wilson_suite as ws

evv_exp = ws.fixtures.evv_experiment()
terms = ws.derive.derive.get_fully_enhanced_terms(experiment=evv_exp)

axis_choice = make_SpectralAxisSet({'A': [-1], 'B': [-1, 2]})
translated_terms = ws.derive.term_var_translate.translate_terms_to_axis_variables(terms, axis_choice)

flat_dict_orig = ws.utils.termdict_from_symb_term.derived_terms_flat(terms, tolistonly=False)
flat_dict = ws.utils.termdict_from_symb_term.derived_terms_flat(translated_terms, tolistonly=False)
flat_list = ws.utils.termdict_from_symb_term.derived_terms_flat(translated_terms, tolistonly=True)


def test_show_valid_axis_combs():
    print('\n')
    show_valid_axis_combs(evv_exp.valid_axis_combs)

def test_make_SpectralAxisSet():
    axis_choice = make_SpectralAxisSet({'A': [-1], 'B': [-1, 2]})
    print(axis_choice)

def test_show_term_latex():
    for k,v in flat_dict.items():
        print('\n', k)
        print(show_term_latex(v))

    for k,v in flat_dict_orig.items():
        print('\n', k)
        print(show_term_latex(v))
