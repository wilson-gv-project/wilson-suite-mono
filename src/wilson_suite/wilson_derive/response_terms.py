import copy
import json
from fractions import Fraction
from dataclasses import dataclass

from wilson_suite.wilson_derive.abstractions import PolProp, VibDiffTerm, ResonanceCondition, QOperator, \
    HarmOscStateSymbolic, TransitionIntegral

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from wilson_suite.wilson_experiment.experiment_abstractions import VibExperiment
    from wilson_suite.wilson_main.abstractions import MolecularProperty
    from wilson_suite.wilson_derive.term_var_translate import SpectralAxisSet

class VibPerturbedTerm:
    """
    Class to represent post-Hermite treatment term: Instances of this class are the main end result of a typical
    wilson-derive run.
    """

    def __init__(self, coeff: Fraction, props: list[PolProp], freqterms: list[VibDiffTerm], res: list[ResonanceCondition]):
        """
        A term consists of the following attributes:

        coeff: Coefficient of term (Fraction)

        props: Molecular properties: At the stage of derivation where the present class is relevant, all of these
        would involve at least one order of geometrical differentiation (list of PolProp instances)

        freqterms: Prefactors of inverse differences between vibrational energy levels as identified
        by the derivation (list of VibDiffTerm instances)

        res: Resonance conditions as identified by the derivation (list of ResonanceCondition instances)
        """

        if not (isinstance(coeff, Fraction) ):
            raise TypeError('VibPerturbedTerm coeff must be a Fraction instance')
        self.coeff = coeff

        if not (all([isinstance(i, PolProp) for i in props])):
            raise TypeError('All transition integral arguments in VibPerturbedTerm must be transitionIntegral instances')
        self.props = props

        if not (all([isinstance(i, VibDiffTerm) for i in freqterms])):
            raise TypeError('All frequency arguments in VibPerturbedTerm must be vibDiffTerm instances')
        self.freqterms = freqterms

        if not (all([isinstance(i, ResonanceCondition) for i in res])):
            raise TypeError('All resonance condition arguments in VibPerturbedTerm must be resonanceCondition instances')
        self.res = res

        # Boolean: Has this term been sorted according to the below sort method?
        self.was_sorted = False

        # Hash (currently indeterminate)
        self.hsh = None

        # to have info about the term in it
        self.anharmonicity = None

        self.note = None

    def __repr__(self):
        return f"VibPerturbedTerm(coeff = {self.coeff}, props = {self.props}, freqterms = {self.freqterms}, res = {self.res})"

    def tellNonSummSummIndices(self):

        summation_indices = []
        non_summation_indices = []

        for i in self.res:

            for j in i.diff.sl.q:
                if not j in non_summation_indices:
                    non_summation_indices.append(j)

            for j in i.diff.sr.q:
                if not j in non_summation_indices:
                    non_summation_indices.append(j)

        candidate_summation_indices = []

        for i in self.freqterms:

            for j in i.sl.q:
                if not j in candidate_summation_indices:
                    candidate_summation_indices.append(j)

            for j in i.sr.q:
                if not j in candidate_summation_indices:
                    candidate_summation_indices.append(j)

        for i in candidate_summation_indices:

            if not i in non_summation_indices:
                summation_indices.append(i)


        return non_summation_indices, summation_indices


    def nmRenameAndInternalResort(self, mask: dict):
        """
        Take a mask (dictionary of single-character key: value pairs) in an ordering to be replaced by the canonical normal mode
        index list: Example: mask is ['b': 'a', 'a': 'b', 'c': 'c'] -> Replace every reference to 'b' in self with 'a',
        replace every (original) 'a' with 'b' and leave 'c' unchanged
        """

        # In properties
        for i in range(len(self.props)):
            for j in range(len(self.props[i].inds)):
                self.props[i].inds[j] = mask[self.props[i].inds[j]]
            self.props[i].inds = sorted(self.props[i].inds)

        # In (non-resonance condition) frequency terms
        for i in range(len(self.freqterms)):
            for j in range(len(self.freqterms[i].sl.q)):
                self.freqterms[i].sl.q[j] = mask[self.freqterms[i].sl.q[j]]
            self.freqterms[i].sl.q = sorted(self.freqterms[i].sl.q)
            for j in range(len(self.freqterms[i].sr.q)):
                self.freqterms[i].sr.q[j] = mask[self.freqterms[i].sr.q[j]]
            self.freqterms[i].sr.q = sorted(self.freqterms[i].sr.q)

        # In resonance conditions
        for i in range(len(self.res)):
            for j in range(len(self.res[i].diff.sl.q)):
                self.res[i].diff.sl.q[j] = mask[self.res[i].diff.sl.q[j]]
            self.res[i].diff.sl.q = sorted(self.res[i].diff.sl.q)
            for j in range(len(self.res[i].diff.sr.q)):
                self.res[i].diff.sr.q[j] = mask[self.res[i].diff.sr.q[j]]
            self.res[i].diff.sr.q = sorted(self.res[i].diff.sr.q)

    def sort(self, nm_inds):
        """
        Sort and possibly rename indices in term to put in canonical term:
        - Sort resonance conditions in increasing order of number of perturbing frequencies
        - Rename normal mode indices according to encountered state labels in resonance conditions (currently leaving
        further indices in arbitrary order)
        - Sort operator references in properties in increasing numerolexical order
        - Sort properties among each order in increasing order of geometric differentiation
        - Then sort properties at same order of differentiation according to operator (numerolexical) ordering
        - Sort frequency difference terms internally to have greatest number of quanta in bra state
        - Sort terms tied wrt. previous sorting according to lexical order of state tuples
        - Sort freq diff terms wrt. each other according to lexical ordering of bras (sort tied terms by ket lex. ordering)

        nm_inds: List of canonically ordered normal mode indices (['a', 'b', ...])
        """

        # Sort resonance conditions in increasing order of number of perturbing frequencies
        self.res = sorted(self.res, key=lambda j:len(j.pf))

        # Make normal mode index sorting mask:
        # First gather according to encountered state labels in resonance conditions
        # Then get remaining indices as encountered in derivatives (FIXME: Likely degree of freedom)

        nm_labels = []

        for i in self.res:

            for j in i.diff.sl.q:
                if not j in nm_labels:
                    nm_labels.append(j)

            for j in i.diff.sr.q:
                if not j in nm_labels:
                    nm_labels.append(j)

        for i in self.props:
            for j in i.inds:
                if not j in nm_labels:
                    nm_labels.append(j)

        # If the label progression doesn't correspond to the canonical progression, then rename to make canonical
        # and re-sort (FIXME: Could be remaining sorting "slack" after this)
        if not (nm_labels == nm_inds[:len(nm_labels)]):
            nm_mask = {nm_labels[i]: nm_inds[i] for i in range(len(nm_labels))}
            self.nmRenameAndInternalResort(nm_mask)

        # Internal sorting: Operators in each property in increasing order
        for i in self.props:
            i.ops = sorted(i.ops, key=lambda j: j.o)

        # Inter-term sorting

        # First sort properties according to increasing order of differentiation
        self.props = sorted(self.props, key=lambda j: j.dord)

        # Then sort properties at same order of differentiation according to operators
        # NOTE: I believe this should take care of both number and lexical ordering

        m = 0
        dord_starts = [0]

        for i in range(len(self.props)):
            if self.props[i].dord > self.props[dord_starts[m]].dord:
                dord_starts.append(i)
                m += 1

        dord_starts.append(len(self.props) - 1)

        for i in range(len(dord_starts)):
            # No need to sort if at the last entry
            if not dord_starts[i] == len(self.props) - 1:
                # No need to sort a length 1 category
                if not dord_starts[i + 1] == dord_starts[i] + 1:
                    self.props[dord_starts[i]:dord_starts[i+1]] = sorted(self.props[dord_starts[i]:dord_starts[i+1]],
                                                                         key=lambda j: [k.o for k in j.ops])

        # Sort all freq term bras and kets

        # Internal sorting (assumes all individual nm tuples already sorted)

        # Make all freq diff terms have the greatest number of quanta in bra - flip and change sign if necessary
        for i in range(len(self.freqterms)):
            if len(self.freqterms[i].sl.q) < len(self.freqterms[i].sr.q):
                self.coeff *= Fraction(-1)
                tmp = copy.deepcopy(self.freqterms[i].sl)
                self.freqterms[i].sl = copy.deepcopy(self.freqterms[i].sr)
                self.freqterms[i].sr = tmp

            # Sort tied terms according to nm indices
            elif len(self.freqterms[i].sl.q) == len(self.freqterms[i].sr.q):
                if not (sorted([self.freqterms[i].sl.q, self.freqterms[i].sr.q]) ==
                               [self.freqterms[i].sl.q, self.freqterms[i].sr.q]):
                    self.coeff *= Fraction(-1)
                    tmp = copy.deepcopy(self.freqterms[i].sl)
                    self.freqterms[i].sl = copy.deepcopy(self.freqterms[i].sr)
                    self.freqterms[i].sr = tmp

        # Inter-term sorting

        # Sort freq diff terms first according to ket indices
        self.freqterms = sorted(self.freqterms, key=lambda j:j.sr.q)

        # Then sort tied freq diff terms according to bra indices
        ketq_starts = [0]

        m = 0

        for i in range(len(self.freqterms)):
            if not (self.freqterms[i].sr.q == self.freqterms[ketq_starts[m]].sr.q):
                ketq_starts.append(i)
                m += 1

        ketq_starts.append(len(self.freqterms))

        for i in range(len(ketq_starts)):
            if not ketq_starts[i] == len(self.freqterms):
                self.freqterms[ketq_starts[i]:ketq_starts[i + 1]] = \
                    sorted(self.freqterms[ketq_starts[i]:ketq_starts[i + 1]], key=lambda j: j.sl.q)

        # Mark term as having been sorted
        # WARNING: This flag doesn't update if any subsequent changes are made
        self.was_sorted = True

    def h(self, also_sort: bool=False, nm_inds: list=None) -> int:
        """
        Hashing function

        also_sort: Also sort term before hash is calculated?
        nm_inds: List of normal mode indices if sorting
        """

        if not(self.was_sorted) and not(also_sort):
            raise AssertionError('Term for which hash was requested has not been sorted')

        if also_sort:
            self.sort(nm_inds)

        # Getting hashes of constituent parts
        props_h = tuple([hash(i) for i in self.props])
        ft_h = tuple([i.h() for i in self.freqterms])
        res_h = tuple([i.h() for i in self.res])

        # Combine constituent hashes for collective hash for this term
        self.hsh = hash((props_h, ft_h, res_h))

        return self.hsh

    def full_enhancement_possible(self, magn_conditions=None) -> bool:
        """
        Determine: Given the setup/requested frequency ranges, is it possible for this term to become resonant within these
        ranges with respect to all of its resonance conditions?

        FIXME: Currently only for sole harmonic oscillator state, may later need extension for variationally resolved states

        magn_conditions: List [[A = signed freq i, B = signed freq j], ...] signifying that
        B - A > 0
        """

        # Discernment by frequency ranges implementation not ready
        #
        #if field_freq_ranges is not None:
        #
        #    for i in self.res:
        #        if not(i.couldBeResonantWithFieldByRanges(field_freq_ranges)):
        #            return False
        #
        #    return True

        if magn_conditions is not None:

            prev_res = None
            for i in self.res:
                if not (i.couldBeResonantWithFieldByConditions(magn_conditions, given_prev_res=prev_res)):
                    return False
                if not (i.couldBeResonantWithFieldByConditions(magn_conditions, given_prev_res=None)):
                    return False
                prev_res = copy.deepcopy(i)

            return True

        # If no information given, return True to keep term
        else:
            return True

    def present(self):
        """
        Formatted printing of own attributes
        """
        print(' >> VibPerturbedTerm presents:')

        print('Coefficient:', self.coeff)
        print('Properties:')
        for i in self.props:
            i.present()
        print('Freqterms:')
        for i in self.freqterms:
            i.present()
        print('Resonances:')
        for i in self.res:
            i.present()

        print('\nHas attributes: coeff, freqterms, res, was_sorted, props, hsh'
              '\nHas methods: nmRenameAndResort, sort, h, full_enhancement_possible, present'
              '\nPresenting also elements of: self.props, self.freqterms, self.res')

    def present_better(self):
        """
        Formatted printing of own attributes
        """
        print('\n >> VibPerturbedTerm presents:')

        print('Coefficient:', self.coeff)
        print('Properties:')
        for i in self.props:
            print('    ', i)
        print('Freqterms:')
        for i in self.freqterms:
            print('    ', i)
        print('Resonances:')
        for i in self.res:
            print('    ', i)

        print('\nHas attributes: coeff, freqterms, res, was_sorted, props, hsh'
              '\nHas methods: nmRenameAndResort, sort, h, full_enhancement_possible, present'
              '\nPresenting also elements of: self.props, self.freqterms, self.res')

    def to_latex(self, part=None):
        """
        """
        res_conditions_denom = ''.join([rc.to_latex() for rc in self.res])
        if res_conditions_denom == '':
            res_conditions_str = ''
        else:
            res_conditions_str = rf'\frac{{1}}{{{res_conditions_denom}}}'

        coefficients_str = rf'\frac{{{self.coeff.numerator}}}{{{self.coeff.denominator}}}'
        properties_str = ''.join([p.to_latex() for p in self.props])

        freqterms_denom = ''.join([rf'\omega_{{{vd.to_latex()}}}' for vd in self.freqterms])
        if freqterms_denom == '':
            freqterms_str = ''
        else:
            freqterms_str = rf'\frac{{{1}}}{{{freqterms_denom}}}'

        if part is not None:
            if part=='res':
                return res_conditions_str
            elif part=='coeff':
                return coefficients_str
            elif part=='props':
                return properties_str
            elif part=='freqterms':
                return freqterms_str
        else:
            return coefficients_str + freqterms_str + properties_str + res_conditions_str

    def to_dict(self) -> dict:
        """Convert VibPerturbedTerm to a dictionary for JSON serialization"""
        return {
            "coeff": [self.coeff.numerator, self.coeff.denominator],
            "props": [
                {
                    "ops": [{"o": op.o, "op_type": op.op_type, "ax": op.ax} for op in p.ops],
                    "dord": p.dord,
                    "inds": p.inds
                } for p in self.props
            ],
            "freqterms": [
                {
                    "sl": {"q": ft.sl.q},
                    "sr": {"q": ft.sr.q},
                    "is_pert_wf_diff": ft.is_pert_wf_diff
                } for ft in self.freqterms
            ],
            "res": [
                {
                    "diff": {
                        "sl": {"q": r.diff.sl.q},
                        "sr": {"q": r.diff.sr.q}
                    },
                    "pf": r.pf,
                    "id": r.id
                } for r in self.res
            ],
            "was_sorted": self.was_sorted
        }

    
    def to_str(self, nm_inds_for_sort=None):
        """
        Make a string identifyer - insteead of hash()

        """
        if not self.was_sorted:
            self.sort(nm_inds_for_sort)

        coeff = f'{self.coeff.numerator}/{self.coeff.denominator}'
        
        from wilson_suite.wilson_utils.prop_trivname import prop_trivname
        props = [f'{prop_trivname(ord_el=len(prop.ops), ord_geo=prop.dord)}_({','.join(prop.inds)})' for prop in self.props]

        freqterms = [f"<{ft.sl.q},{ft.sr.q}>" for ft in self.freqterms]

        res = [f'(<{r.diff.sl.q}{r.diff.sr.q}> - {r.pf} -iG)' for r in self.res]

        return coeff + ' * ' + ' * '.join(props) + ' / ' + ' / '.join(freqterms) + ' / ' + ' / '.join(res)


    @classmethod
    def from_dict(cls, data: dict) -> 'VibPerturbedTerm':
        """Create VibPerturbedTerm from a dictionary"""
        # Reconstruct coefficient
        coeff = Fraction(data["coeff"][0], data["coeff"][1])

        # Reconstruct properties
        props = []
        for p in data["props"]:
            ops = [QOperator(op["o"], op["op_type"], op["ax"]) for op in p["ops"]]
            prop = PolProp(ops, p["dord"])
            prop.inds = p["inds"]
            props.append(prop)

        # Reconstruct frequency terms
        freqterms = []
        for ft in data["freqterms"]:
            sl = HarmOscStateSymbolic(ft["sl"]["q"])
            sr = HarmOscStateSymbolic(ft["sr"]["q"])
            freqterms.append(VibDiffTerm(sl, sr, ft["is_pert_wf_diff"]))

        # Reconstruct resonance conditions
        res = []
        for r in data["res"]:
            sl = HarmOscStateSymbolic(r["diff"]["sl"]["q"])
            sr = HarmOscStateSymbolic(r["diff"]["sr"]["q"])
            diff = VibDiffTerm(sl, sr)
            res.append(ResonanceCondition(diff, r["pf"], r["id"]))

        # Create instance
        instance = cls(coeff, props, freqterms, res)
        instance.was_sorted = data["was_sorted"]
        return instance

    def to_json(self, filepath: str = None) -> str:
        """Convert to JSON string or save to file"""
        json_str = json.dumps(self.to_dict(), indent=2)
        if filepath:
            with open(filepath, 'w') as f:
                f.write(json_str)
        return json_str

    @classmethod
    def from_json(cls, json_str: str = None, filepath: str = None) -> 'VibPerturbedTerm':
        """Create instance from JSON string or file"""
        if filepath:
            with open(filepath) as f:
                data = json.load(f)
        else:
            data = json.loads(json_str)
        return cls.from_dict(data)

    @classmethod
    def save_many_to_json(cls, terms: list['VibPerturbedTerm'], filepath: str):
        """Save multiple terms to a single JSON file"""
        terms_data = [term.to_dict() for term in terms]
        with open(filepath, 'w') as f:
            json.dump(terms_data, f, indent=2)

    @classmethod
    def load_many_from_json(cls, filepath: str) -> list['VibPerturbedTerm']:
        """Load multiple terms from a JSON file"""
        with open(filepath) as f:
            terms_data = json.load(f)
        return [cls.from_dict(term_data) for term_data in terms_data]


class VibContribTerm:
    """
    Class to represent a vibrational term in SOS expression up to but not including Hermite treatment
    """

    def __init__(self, coeff: Fraction, ints: list[TransitionIntegral], res: list[ResonanceCondition]):
        """
        A term consists of the following attributes:

        coeff: Coefficient of term (Fraction)

        ints: Transition integrals involving a bra state, one or more operators, and a ket state

        res: Resonance conditions as identified by the derivation (list of ResonanceCondition instances)
        """

        if not (isinstance(coeff, Fraction) ):
            raise TypeError('vibContribTerm coeff must be a Fraction instance')
        self.coeff = coeff

        if not (all([isinstance(i, TransitionIntegral) for i in ints])):
            raise TypeError('All transition integral arguments in VibContribTerm must be TransitionIntegral instances')
        self.ints = ints

        if not (all([isinstance(i, ResonanceCondition) for i in res])):
            raise TypeError('All resonance condition arguments in VibContribTerm must be ResonanceCondition instances')
        self.res = res

        # Initialize frequency difference terms in anticipation of possible addition
        # from mechanical anharmonicity handling
        self.freqdiff = []

    def addFreqTerm(self, new_term: VibDiffTerm):
        """
        Add (inverted) frequency difference term from mechanical anharmonicity handling
        """

        if not (isinstance(new_term, VibDiffTerm)):
            raise TypeError('VibContribTerm frequency difference terms must be VibDiffTerm instances')

        self.freqdiff.append(new_term)

    def dressWithPulseInteractions(self, int_seq: tuple[dict]):
        """
        Substitute interaction dummy indices in my resonance conditions
        and in my polarization property axes with those of a specific pulse interaction sequence

        int_seq: tuple of dictionaries ({interaction #1 pulse label: sign}, {int. #2 pulse label: sign}, ...)
        """
        # FIXME: Must get tests as part of averaging work

        int_ids = [list(i.keys())[0] for i in int_seq]

        # Dress perturbing frequencies of resonance conditions
        for i in self.res:
            for j in range(len(i.pf)):
                i.pf[j] = int_ids[i.pf[j] - 1] * int_seq[j][int_ids[i.pf[j]  - 1]]

        # Dress operators of polarization properties
        for i in self.ints:
            for j in i.prop.ops:
                if not j.o == 0:
                    j.o = int_ids[j.o - 1]


    def allUVCancels(self, cfs_uv: dict, tol: float=0.0) -> bool:
        """
        Boolean function: Do all UV/VIS carrier frequency components sum to zero
        (leaving only potentially vibrationally resonant values) for my resonance conditions?

        cfs_uv: Dictionary {pulse label: UV/VIS freq. component, ...}
        tol: Tolerance (allow a verdict of "cancels" if sum is within tolerance)
        """

        all_cancel = True

        for i in self.res:
            all_cancel = all_cancel and i.uvCancels(cfs_uv, tol)

        return all_cancel


    def allElRspEpochContained(self, epochs, op_ind_omega) -> bool:
        """
        Boolean function: Do all the individual operator tuples of my (electronic) response functions
        belong to pulses corresponding to the same epoch?

        epochs: List of epochs and which pulses (pulse IDs) are in each
        op_ind_omega: Which operator is the "omega" (interaction w/detected field) operator?
        """

        epoch_contained = True

        for i in self.ints:
            epoch_contained = epoch_contained and i.prop.epochContained(epochs, op_ind_omega)

        return epoch_contained

    def present(self):
        """
        Formatted printing of own attributes
        """

        print('Coefficient:', self.coeff)

        for i in self.ints:
            i.present()

        for i in self.res:
            i.present()

        for i in self.freqdiff:
            i.present()

        print('')



@dataclass
class VibPertTermsCollection:
    """
    term_dict -- as orriginally derived terms - from ws.derive.derive.get_fully_enhanced_terms() or a list already
    experiment -- provenance of these terms
    axes_choice -- tracking if translated to axes choice, and which ones
    """
    term_dict: dict[int, dict[tuple, VibPerturbedTerm]]
    experiment: 'VibExperiment'
    axes_choice: 'SpectralAxisSet' = None
    magn_conditions: tuple = None

    def show_as(self, format: str='latex', part: str = None):
        """
        format - latex, str
        """
        if isinstance(self.term_dict, dict):
            flat_terms_dict = self.make_flat(self.term_dict)
        elif isinstance(self.term_dict, list):
            flat_terms_dict = self.term_dict

        if format=='latex':

            for id, term in flat_terms_dict.items():
                print('&'+term.to_latex(part=part) + rf' \\ % term_id {id}')
        
        elif format =='str':
        
            for id, term in flat_terms_dict.items():
                print(term.to_str())

        else:
            raise NotImplementedError('This format is no implemented for VibPertTermsCollection.show_as()')


    def make_flat(self, as_list=False) -> dict[str, VibPerturbedTerm] | list[VibPerturbedTerm]:
        """
        returns flat either dict or list, not a VibPertTermsCollection
        """
        from wilson_suite.wilson_utils.termdict_from_symb_term import derived_terms_flat
        return derived_terms_flat(self.term_dict, as_list)


    def required_data(self, freqs: str='static') -> tuple[int, list['MolecularProperty']]:
        """
        find_props here returns a list of MolecularProperty objects - is it too much at this point?
        """
        from wilson_suite.wilson_main.main_functions import find_max_state_lvl, find_props

        return find_props(self.term_dict, freqs), find_max_state_lvl(self.term_dict)
    

    def available_axes_choices(self):
        self.experiment.tell_axis_options()

    
    def translate_to_ax_choice(self, axes_choice, magn_conditions_translated=None) -> 'VibPertTermsCollection':
        """
        returns a new VibPertTermsCollection object

        should have some default axes_choice?
        """

        from wilson_suite.wilson_derive.term_var_translate import translate_terms_to_axis_variables
        return VibPertTermsCollection(translate_terms_to_axis_variables(self.term_dict, axes_choice), 
                                      experiment=self.experiment,
                                      axes_choice=axes_choice,
                                      magn_conditions=magn_conditions_translated)
