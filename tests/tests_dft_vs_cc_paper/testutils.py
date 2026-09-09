from wilson_suite.wilson_intensities.amplitudes.spectrum_composition import SpectralFeature, ResLocGeoObject, SpectralWindow
from wilson_suite.wilson_intensities.amplitudes.term_parts import TermParametersChoice, ResonanceMotif, ResonanceCondition, ParameterSet
import wilson_suite as ws
from wilson_suite.wilson_main import abstractions as wm_abst

class MakeObjects:
    
    @staticmethod
    def mk_feature_single_onetermid() -> SpectralFeature:
        """
        with one random term_contributions for now
        """
        rcs = [ResonanceCondition.make_from_tuples(left_state=(), right_state=('a',), pert_freqs=('-A',)),
               ResonanceCondition.make_from_tuples(left_state=('a', 'b'), right_state=('b',), pert_freqs=('B',))]
        res_motif = ResonanceMotif(rcs)
        
        param_set = ParameterSet({'a': 0, 'b': 1})

        term_contrib = TermParametersChoice(res_motif=res_motif,
                                            states_parameters=(param_set,),
                                            term_ids=tuple([hash(res_motif)]))
        
        feat = SpectralFeature(location=ResLocGeoObject({'A': 1119.5, 'B': 2921.}), 
                               lineshape_parameter=4.5, 
                               amplitude_coeff=-1.12e-06, 
                               term_contributions=(term_contrib,))
        features = SpectralFeature.dress_these_with_boxes([feat],
                                                          max_intensity=feat.get_intensity()+1., 
                                                          min_intensity=feat.get_intensity()-1.0e3,
                                                          box_range_safety_margin=0.1,
                                                          scale_wrt_max_intensity=True,
                                                          minimum_box_padding=30.0,
                                                          )
        return features[0]

    @staticmethod
    def mk_feature_single_multitermids() -> SpectralFeature:
        """
        with random 1 term_contributions with several term_ids with same res_motif (that is how features are).
        Number of term_contributions is usually one (they differ by res motif and parameters choices).

        """
        rcs = [ResonanceCondition.make_from_tuples(left_state=(), right_state=('a',), pert_freqs=('-A',)),
               ResonanceCondition.make_from_tuples(left_state=('a', 'b'), right_state=('b',), pert_freqs=('B',))]
        res_motif = ResonanceMotif(rcs)
        
        param_set = ParameterSet({'a': 0, 'b': 1})

        term_contrib = TermParametersChoice(res_motif=res_motif,
                                            states_parameters=(param_set,),
                                            term_ids=tuple([111, 112, 113, 114]))
        
        feat = SpectralFeature(location=ResLocGeoObject({'A': 1119.5, 'B': 2921.}), 
                               lineshape_parameter=4.5, 
                               amplitude_coeff=-1.12e-06, 
                               term_contributions=(term_contrib,))
        features = SpectralFeature.dress_these_with_boxes([feat],
                                                          max_intensity=feat.get_intensity()+1., 
                                                          min_intensity=feat.get_intensity()-1.0e3,
                                                          box_range_safety_margin=0.1,
                                                          scale_wrt_max_intensity=True,
                                                          minimum_box_padding=30.0,
                                                          )
        return features[0]

    @staticmethod
    def mk_features_non_ovrl() -> list[SpectralFeature]:

        rcs1 = [ResonanceCondition.make_from_tuples(left_state=(), right_state=('a',), pert_freqs=('-A',)),
               ResonanceCondition.make_from_tuples(left_state=('a', 'b'), right_state=('b',), pert_freqs=('B',))]
        res_motif1 = ResonanceMotif(rcs1)

        param_set1 = ParameterSet({'a': 0, 'b': 1})
        term_contrib1 = TermParametersChoice(res_motif=res_motif1,
                                            states_parameters=(param_set1,),
                                            term_ids=tuple([hash(res_motif1)]))

        rcs2 = [ResonanceCondition.make_from_tuples(left_state=(), right_state=('a',), pert_freqs=('-A',)),
               ResonanceCondition.make_from_tuples(left_state=('b',), right_state=('a',), pert_freqs=('B',))]
        res_motif2 = ResonanceMotif(rcs2)

        param_set2 = ParameterSet({'a': 1, 'b': 2})
        term_contrib2 = TermParametersChoice(res_motif=res_motif2,
                                            states_parameters=(param_set2,),
                                            term_ids=tuple([hash(res_motif2)]))
        
        feat1 = SpectralFeature(location=ResLocGeoObject({'A': 1119.5, 'B': 2921.}), 
                                lineshape_parameter=4.5, 
                                amplitude_coeff=-4.32e-06, 
                                term_contributions=(term_contrib1,))
        feat2 = SpectralFeature(location=ResLocGeoObject({'A': 964., 'B': 270.}), 
                                lineshape_parameter=4.5, 
                                amplitude_coeff=-1.12e-05, 
                                term_contributions=(term_contrib2,))
        
        features = SpectralFeature.dress_these_with_boxes([feat1, feat2],
                                                          max_intensity=feat2.get_intensity()+1., 
                                                          min_intensity=feat1.get_intensity()-1.0e3,
                                                          box_range_safety_margin=0.1,
                                                          scale_wrt_max_intensity=True,
                                                          minimum_box_padding=30.0,
                                                          )
        return features

    @staticmethod
    def mk_features_ovrl() -> list[SpectralFeature]:

        rcs1 = [ResonanceCondition.make_from_tuples(left_state=(), right_state=('a',), pert_freqs=('-A',)),
               ResonanceCondition.make_from_tuples(left_state=('a', 'b'), right_state=('b',), pert_freqs=('B',))]
        res_motif1 = ResonanceMotif(rcs1)

        param_set1 = ParameterSet({'a': 0, 'b': 1})
        term_contrib1 = TermParametersChoice(res_motif=res_motif1,
                                            states_parameters=(param_set1,),
                                            term_ids=tuple([hash(res_motif1)]))

        rcs2 = [ResonanceCondition.make_from_tuples(left_state=(), right_state=('a',), pert_freqs=('-A',)),
                ResonanceCondition.make_from_tuples(left_state=('b',), right_state=('a',), pert_freqs=('B',))]
        res_motif2 = ResonanceMotif(rcs2)

        param_set2 = ParameterSet({'a': 0, 'b': 1})
        term_contrib2 = TermParametersChoice(res_motif=res_motif2,
                                            states_parameters=(param_set2,),
                                            term_ids=tuple([hash(res_motif2)]))
        
        feat1 = SpectralFeature(location=ResLocGeoObject({'A': 1000., 'B': 520.}), 
                                lineshape_parameter=4.5, 
                                amplitude_coeff=-4.32e-06, 
                                term_contributions=(term_contrib1,))
        feat2 = SpectralFeature(location=ResLocGeoObject({'A': 1000., 'B': 500.}), 
                                lineshape_parameter=4.5, 
                                amplitude_coeff=-1.12e-05, 
                                term_contributions=(term_contrib2,))
        features = SpectralFeature.dress_these_with_boxes([feat1, feat2],
                                                          max_intensity=feat2.get_intensity()+1., 
                                                          min_intensity=feat1.get_intensity()-1.0e3,
                                                          box_range_safety_margin=0.1,
                                                          scale_wrt_max_intensity=True,
                                                          minimum_box_padding=30.0,
                                                          )
        assert features[0].feat_box.overlaps(features[1].feat_box)
        return features
    
    @staticmethod
    def mk_gridregion_for_specwindow(specwindow: SpectralWindow, coords: dict):
        """

        coords is like, e.g.:
        {'A': array([[1089.5, 1089.5, 1089.5, 1089.5, 1089.5],
                     [1103.9, 1103.9, 1103.9, 1103.9, 1103.9],
                     [1118.3, 1118.3, 1118.3, 1118.3, 1118.3],
                     [1132.7, 1132.7, 1132.7, 1132.7, 1132.7],
                     [1147.1, 1147.1, 1147.1, 1147.1, 1147.1]]), 
        'B': array([[2891. , 2905.4, 2919.8, 2934.2, 2948.6],
                    [2891. , 2905.4, 2919.8, 2934.2, 2948.6],
                    [2891. , 2905.4, 2919.8, 2934.2, 2948.6],
                    [2891. , 2905.4, 2919.8, 2934.2, 2948.6]])}
        """
        from wilson_suite.wilson_intensities.amplitudes.grid_manager_evaluator import GridRegion, RectangularDomain
        domain = RectangularDomain(box=specwindow.box,
                                   full_features=specwindow.full_features,
                                   contrib_features=specwindow.contrib_features)

        return GridRegion(domain=domain,
                          coords=coords)
    

    @staticmethod
    def mk_gridregion_for_feature(feature, coords):
        return


    @staticmethod
    def mk_vibstates_states(closer=False):
        """
        rcs = [ResonanceCondition.make_from_tuples(left_state=('',), right_state=('a',), pert_freqs=('-A',)),
               ResonanceCondition.make_from_tuples(left_state=('a', 'b'), right_state=('b',), pert_freqs=('-B',))]
        res_motif = ResonanceMotif(rcs)
        
        param_set = ParameterSet({'a': 0, 'b': 1})

            'A': 1119.5, 'B': 2921.
        1119.5 = a 0
        2921 = 0,1-1
        """
        if not closer:
            states = (
                wm_abst.VibState(harm_quanta_coeffs={(0,):1.}, state_label='0', energy=1119.5, harmonic_WF=True),
                wm_abst.VibState(harm_quanta_coeffs={(1,):1.}, state_label='1', energy=964., harmonic_WF=True),
                wm_abst.VibState(harm_quanta_coeffs={(2,):1.}, state_label='2', energy=1234., harmonic_WF=True),
                wm_abst.VibState(harm_quanta_coeffs={(0, 0):1.}, state_label='0,0', energy=1864., harmonic_WF=False),
                wm_abst.VibState(harm_quanta_coeffs={(0, 1):1.}, state_label='0,1', energy=3885., harmonic_WF=False),
                wm_abst.VibState(harm_quanta_coeffs={(1, 1):1.}, state_label='1,1', energy=2368., harmonic_WF=False),
                wm_abst.VibState(harm_quanta_coeffs={(1, 2):1.}, state_label='1,2', energy=2360., harmonic_WF=False),
                wm_abst.VibState(harm_quanta_coeffs={(2, 2):1.}, state_label='2,2', energy=2362., harmonic_WF=False),
                wm_abst.VibState(harm_quanta_coeffs={(0, 2):1.}, state_label='0,2', energy=2274., harmonic_WF=False),

                wm_abst.VibState(harm_quanta_coeffs={(0, 0, 0):1.}, state_label='0,0,0', energy=2685., harmonic_WF=False),
                wm_abst.VibState(harm_quanta_coeffs={(1, 1, 1):1.}, state_label='1,1,1', energy=3581., harmonic_WF=False),
                wm_abst.VibState(harm_quanta_coeffs={(2, 2, 2):1.}, state_label='2,2,2', energy=3690., harmonic_WF=False),
                wm_abst.VibState(harm_quanta_coeffs={(0, 1, 2):1.}, state_label='0,1,2', energy=3742., harmonic_WF=False),
                wm_abst.VibState(harm_quanta_coeffs={(0, 1, 1):1.}, state_label='0,1,1', energy=3318., harmonic_WF=False),
                wm_abst.VibState(harm_quanta_coeffs={(0, 2, 2):1.}, state_label='0,2,2', energy=3680., harmonic_WF=False),
                wm_abst.VibState(harm_quanta_coeffs={(0, 0, 1):1.}, state_label='0,0,1', energy=3155., harmonic_WF=False),
                wm_abst.VibState(harm_quanta_coeffs={(0, 0, 2):1.}, state_label='0,0,2', energy=3498., harmonic_WF=False),
                wm_abst.VibState(harm_quanta_coeffs={(1, 1, 2):1.}, state_label='1,1,2', energy=3594., harmonic_WF=False),
                wm_abst.VibState(harm_quanta_coeffs={(1, 2, 2):1.}, state_label='1,2,2', energy=3642., harmonic_WF=False),
            )
        else:
            states = (
                wm_abst.VibState(harm_quanta_coeffs={(0,):1.}, state_label='0', energy=1000., harmonic_WF=True),
                wm_abst.VibState(harm_quanta_coeffs={(1,):1.}, state_label='1', energy=1500., harmonic_WF=True),
                wm_abst.VibState(harm_quanta_coeffs={(0, 1):1.}, state_label='0,1', energy=2020., harmonic_WF=False),
            )

        return states

    @staticmethod
    def mk_vibstates_data():
        return 
    
    @staticmethod
    def mk_data_for_eval(list_of_states: list[wm_abst.VibState], 
                         include_states_list,
                         list_of_props, 
                         pulse_polarization_vector):
        """
        data:
        - pulse_polarization_vector (experiment property)
        - vib_data
        - vibdiff_cache
        - data_configs

        raw data:
        - vib_ana_setup - with [states, include_list, number_of_modes]
        - list[MolecularProperty] - with data
        - pulse_polarization_vector (experiment property)

        pkl_file contains a dict with keys:
            ['anharmonic_states', 'harmonic_states', 
            'nc_sqrt_eigval', 'B', 'coriolis', 
            'dipgrad', 'diphess', 'polgrad', 'polhess', 
            'cff', 'cff_rc', 'qff', 'qff_rc', 
            'MolecularSystem', 'DataOriginInfo']
        ---
        vibstates_data = VibStatesData(allstates=tuple(vib_ana_setup.states), 
                                harmonic_osc_states_labels=vib_ana_setup.include_list,
                                number_of_nmodes=vib_ana_setup.number_of_modes)
        vibdiff_cache = VibDiffCache()
        props = MolPropsCollection(properties=props)
        data_and_configs = EvaluationDataAndConfigs(props_data=props,
                                                    vibstates_data=vibstates_data,
                                                    number_of_nmodes=vib_ana_setup.number_of_modes,
                                                    nm_inds_choices=include_list,
                                                    pulse_polarization_vector=pulse_polarization_vector)
        
        """
        """
        data_and_configs = EvaluationDataAndConfigs(props_data=props,
                                            vibstates_data=vibstates_data,
                                            number_of_nmodes=vib_ana_setup.number_of_modes,
                                            nm_inds_choices=include_list,
                                            pulse_polarization_vector=pulse_polarization_vector)
        """
        vib_ana_setup = ws.main.abstractions.VibAnaSetup()
        vib_ana_setup.states = list_of_states
        vib_ana_setup.include_list = include_states_list
        vib_ana_setup.number_of_modes = 3
        vib_ana_setup.nc_sqrt_eigval = {int(i.state_label): i.energy+15.  for i in list_of_states if len(i.state_label.split(','))==1}

        from wilson_suite.wilson_intensities.amplitudes.term_parts import VibStatesData, EvaluationDataAndConfigs
        vibstates_data = VibStatesData(allstates=tuple(vib_ana_setup.states), 
                                harmonic_osc_states_labels=vib_ana_setup.include_list,
                                number_of_nmodes=vib_ana_setup.number_of_modes)
        
        """
        from wilson_suite.wilson_intensities.amplitudes.vibene_differences import VibDiffCache
        vibdiff_cache = VibDiffCache()
        """

        props = ws.main.abstractions.MolPropsCollection(properties=list_of_props)
        data_and_configs = EvaluationDataAndConfigs(props_data=props,
                                                    vibstates_data=vibstates_data,
                                                    number_of_nmodes=vib_ana_setup.number_of_modes,
                                                    nm_inds_choices=include_states_list,
                                                    pulse_polarization_vector=pulse_polarization_vector,
                                                    nc_sqrt_eigval=vib_ana_setup.nc_sqrt_eigval)

        return data_and_configs

    @staticmethod
    def mk_empty_prop_by_trivname(trivial_name):
        """
        electric props are static: 'freq': (0.0, 0.0, 0.0)
        """

        if trivial_name == 'cff':
            return ws.main.abstractions.MolecularProperty(
				{'ops': tuple(['g', 'g', 'g']), 'freq': (0.0, 0.0, 0.0)},
				trivial_name=trivial_name)
        elif trivial_name == 'qff':
            return ws.main.abstractions.MolecularProperty(
				{'ops': tuple(['g', 'g', 'g', 'g']), 'freq': (0.0, 0.0, 0.0, 0.0)},
				trivial_name=trivial_name)
        elif trivial_name == 'B':
            return ws.main.abstractions.MolecularProperty(
				{'ops': tuple(['r']), 'freq': (0.0)},
				trivial_name=trivial_name)
        elif trivial_name == 'coriolis':
            return ws.main.abstractions.MolecularProperty(
				{'ops': tuple(['g', 'g', 'r']), 'freq': (0.0, 0.0, 0.0)},
				trivial_name=trivial_name)
        elif trivial_name in ['dipgrad', 'diphess', 'polgrad', 'polhess']:
            if 'grad' in trivial_name:
                ops_g = ['g']
            elif 'hess' in trivial_name:
                ops_g = ['g', 'g']
            
            if 'dip' in trivial_name:
                ops_f = ['f']
            elif 'pol' in trivial_name:
                ops_f = ['f', 'f']
            
            return ws.main.abstractions.MolecularProperty(
				{'ops': tuple(ops_g+ops_f), 'freq': (0.0, 0.0, 0.0)},
				trivial_name=trivial_name)
        else:
            raise NotImplementedError("hess, dipgrad, diphess, polgrad, polhess, B, coriolis")
    
    @staticmethod
    def mk_prop_with_vals(trivial_name):
        """
        first geo, then electric - dimensions of tensors.
        """
        import numpy as np

        prop = MakeObjects.mk_empty_prop_by_trivname(trivial_name)

        if trivial_name == 'dipgrad':
            # gf
            vals = np.zeros((3, 3))
            vals[0, 0] = 0.1
            vals[1, 2] = -0.1
            vals[2, 1] = -0.1
        
        elif trivial_name == 'diphess':
            # ggf
            vals = np.zeros((3, 3, 3))
            vals[0, 0, 0] = -0.1
            # vals[1, 1, 1] = 3.4
            # vals[0, 1, 1] = vals[1, 0, 1] = 1.2
            vals[2, 0, 1] = vals[0, 2, 1] = -0.1

        elif trivial_name == 'polgrad':
            # gff
            vals = np.zeros((3, 3, 3))
            vals[0, 0, 0] = -0.1
            vals[1, 1, 2] = vals[1, 2, 1] = 0.1
            # vals[2, 1, 0] = vals[2, 0, 1] = 0.1

        elif trivial_name == 'polhess':
            # ggff
            vals = np.zeros((3, 3, 3, 3))
            vals[0, 0, 0, 0] = 0.1
            # vals[0, 1, 0, 1] = vals[1, 0, 0, 1] = vals[0, 1, 1, 0] = vals[1, 0, 1, 0] = -3.
            # vals[1, 1, 1, 1] = 0.5
            # vals[1, 2, 1, 2] = vals[2, 1, 1, 2] = vals[1, 2, 2, 1] = vals[2, 1, 2, 1] = -0.7
            vals[2, 2, 1, 2] = vals[2, 2, 2, 1] = 0.1

        elif trivial_name == 'cff':
            # ggg
            vals = np.zeros((3, 3, 3))
            vals[0, 0, 0] = 0.1
            vals[1, 1, 2] = vals[1, 2, 1] = vals[2, 1, 1] = -0.1
            # vals[0, 0, 1] = vals[0, 1, 0] = vals[1, 0, 0] = -1.3
            # vals[0, 1, 1] = vals[1, 0, 1] = vals[1, 1, 0] = 2.4
            # vals[1, 2, 2] = vals[2, 1, 2] = vals[2, 2, 1] = 2.4
            # vals[0, 1, 2] = vals[2, 0, 1] = vals[1, 2, 0] = vals[0, 2, 1] = vals[2, 1, 0] = vals[1, 0, 2] = -0.1

        """
        if trivial_name == 'dipgrad':
            # gf
            vals = np.zeros((3, 3))
            vals[0, 0] = 1.5
            vals[1, 2] = -0.5
            vals[2, 1] = -2.2

        elif trivial_name == 'diphess':
            # ggf
            vals = np.zeros((3, 3, 3))
            vals[0, 0, 0] = -0.2
            vals[1, 1, 1] = 3.4
            vals[0, 1, 1] = vals[1, 0, 1] = 1.2
            vals[2, 0, 1] = vals[0, 2, 1] = -0.8          # fixed symmetry (was vals[0,2,2])
            vals[1, 2, 0] = vals[2, 1, 0] = 0.3            # new — fixes [1,2,:] and [2,1,:]
            vals[2, 2, 2] = -1.1                            # new — fixes [2,2,:]
        
        elif trivial_name == 'polgrad':
            # gff
            vals = np.zeros((3, 3, 3))
            vals[0, 0, 0] = -0.2
            vals[1, 1, 2] = vals[1, 2, 1] = 1.2
            vals[2, 1, 0] = vals[2, 0, 1] = 0.1
            vals[0, 1, 1] = 0.8                             # new — fixes [0,1,:]
            vals[0, 2, 2] = -0.5                            # new — fixes [0,2,:]
            vals[1, 0, 0] = 1.4                             # new — fixes [1,0,:]
            vals[2, 2, 2] = -0.6                            # new — fixes [2,2,:]
        
        elif trivial_name == 'polhess':
            # ggff
            vals = np.zeros((3, 3, 3, 3))
            vals[0, 0, 0, 0] = 0.1
            vals[0, 1, 0, 1] = vals[1, 0, 0, 1] = vals[0, 1, 1, 0] = vals[1, 0, 1, 0] = -3.
            vals[1, 1, 1, 1] = 0.5
            vals[1, 2, 1, 2] = vals[2, 1, 1, 2] = vals[1, 2, 2, 1] = vals[2, 1, 2, 1] = -0.7
            vals[2, 2, 1, 2] = vals[2, 2, 2, 1] = 1.8
            # new entries:
            vals[0, 0, 1, 1] = 0.9                          # fixes [0,0,1,:]
            vals[0, 0, 2, 2] = -0.3                         # fixes [0,0,2,:]
            vals[0, 1, 2, 2] = vals[1, 0, 2, 2] = 0.6      # fixes [0,1,2,:] and [1,0,2,:]
            vals[0, 2, 0, 2] = vals[2, 0, 0, 2] = vals[0, 2, 2, 0] = vals[2, 0, 2, 0] = 0.4  # fixes [0,2,0,:], [2,0,0,:], [0,2,2,:], [2,0,2,:]
            vals[0, 2, 1, 1] = vals[2, 0, 1, 1] = -0.5     # fixes [0,2,1,:] and [2,0,1,:]
            vals[1, 1, 0, 0] = 1.3                          # fixes [1,1,0,:]
            vals[1, 1, 2, 2] = -0.2                         # fixes [1,1,2,:]
            vals[1, 2, 0, 0] = vals[2, 1, 0, 0] = 0.7      # fixes [1,2,0,:] and [2,1,0,:]
            vals[2, 2, 0, 0] = 1.5                          # fixes [2,2,0,:]
        
        elif trivial_name == 'cff':
            # ggg
            vals = np.zeros((3, 3, 3))
            vals[0, 0, 0] = 2.8
            vals[0, 0, 1] = vals[0, 1, 0] = vals[1, 0, 0] = -1.3
            vals[1, 1, 2] = vals[1, 2, 1] = vals[2, 1, 1] = -0.1
            vals[0, 1, 1] = vals[1, 0, 1] = vals[1, 1, 0] = 2.4
            vals[1, 2, 2] = vals[2, 1, 2] = vals[2, 2, 1] = 2.4
            vals[0, 1, 2] = vals[2, 0, 1] = vals[1, 2, 0] = vals[0, 2, 1] = vals[2, 1, 0] = vals[1, 0, 2] = -0.4

        # elif trivial_name == 'qff':
        #     # gggg
        #     vals = np.zeros((3, 3, 3, 3))
        #     vals[0, 0, 0, 0] = -0.8
        #     vals[0, 0, 1, 1] = vals[0, 1, 0, 1] = vals[1, 0, 0, 1] = vals[1, 1, 0, 0] = vals[0, 1, 1, 0] = vals[1, 0, 1, 0] = 2.4
        """
        prop.addValues(vals)
        
        return prop

def fillStatesData(data_dict):
    """
    data_dict should have key anharmonic_states, and value is a dictionary like:
    {('2',): 1497.353, ('1', '1'): 3556.565, ('1', '4', '5'): 5753.841}
    """
    from wilson_suite.wilson_main import abstractions as wm_abst

    states_list: list[wm_abst.VibState] = []
    states_dict: dict = data_dict.get('anharmonic_states')

    for state, energy in states_dict.items():
        states_list.append(wm_abst.VibState(harm_quanta_coeffs={state: 1.0}, energy=energy, state_label=','.join(state)))

    return states_list

def fillPropsData(data_dict):
    """
    
    """
    final_props = []
    for key in data_dict:
        if key in ['hess', 'dipgrad', 'diphess', 
                   'polgrad', 'polhess', 'B', 'coriolis']:
            final_props.append(MakeObjects.mk_empty_prop_by_trivname(key))
    return final_props



def get_from_pkl_features(pkl_file, lineshape_parameter):
    """
    pkld_file -- 'data_for_tests/FORM_conf1_B3LYP_aug_cc_pVTZ.pkl'

    returns list of features and dict of terms (hash to term)
    """
    from wilson_suite.wilson_utils.serialization import unpickle_smth_from
    unpickled = unpickle_smth_from(pkl_file)
    
    list_vibsstates = fillStatesData(unpickled)

    from wilson_suite.wilson_intensities.amplitudes.term_parts import VibStatesData
    from wilson_suite.wilson_intensities.amplitudes.vibene_differences import VibDiffCache

    include_list = tuple([0, 1, 2])
    vibstates_data = VibStatesData(allstates=tuple(list_vibsstates),
                                   harmonic_osc_states_labels=include_list)
    vibdiff_cache = VibDiffCache()

    from wilson_suite.wilson_intensities.amplitudes.evaluators import get_features_from_terms_for_eval
    from wilson_suite.fixtures import evv_experiment
    evv_exp = evv_experiment()
    terms = ws.derive.derive.get_fully_enhanced_terms(experiment=evv_exp)

    from wilson_suite.wilson_utils.builder_functions import make_SpectralAxisSet
    axes_choice: ws.main.spectrum_abstractions.SpectralAxisSet = make_SpectralAxisSet({'A': [1], 'B': [-1,2]}) # this makes A and B > 0
    sim = ws.main.workflow_abstractions.WilsonSimulation()
    sim.addTerms(terms)
    sim.setAxisChoiceAndTranslateTerms(axes_choice)
    from wilson_suite.wilson_utils.termdict_from_symb_term import derived_terms_flat
    terms_list = derived_terms_flat(sim.terms_in_axis_choice, tolistonly=True)
    
    hashmap = {t.h(): t for t in terms_list}

    features = get_features_from_terms_for_eval(derived_terms=terms_list,
                                                vibstates_data=vibstates_data,
                                                vibdiff_cache=vibdiff_cache, 
                                                lineshape_parameter=lineshape_parameter)
    return features, hashmap