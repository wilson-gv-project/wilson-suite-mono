from wilson_suite.wilson_main.abstractions import DataOriginInfo

def wilson_data_obtainer(requested_data_dict: dict[str,DataOriginInfo], 
                         get_geometry=False, get_displacements=False,
                         reindex_modes=False):
    """
    1. group by origin

    reindex_modes -- to switch nodes numbering for comparisons where it is not the same
    """
    if isinstance(requested_data_dict, dict):
        for k, v in requested_data_dict.items():
            if isinstance(v, DataOriginInfo):
                continue
            else:
                raise TypeError("requested_data_dict should be dict[str,DataOriginInfo]")
    else:
        raise TypeError("requested_data_dict should be dict[str,DataOriginInfo]")
    dict_with_data = {}

    origin_to_req_data: dict[DataOriginInfo, list] = {}

    if get_geometry:
        requested_data_dict.update({'atoms': requested_data_dict['nc_sqrt_eigval'],
                                    'equilibrium_geometry': requested_data_dict['nc_sqrt_eigval']})
    if get_displacements:
        requested_data_dict.update({'normal_modes': requested_data_dict['nc_sqrt_eigval']})

    for k, v in requested_data_dict.items():
        
        if v not in origin_to_req_data:
            origin_to_req_data[v] = [k]

        else:
            origin_to_req_data[v].append(k)

    
    for o in origin_to_req_data:

        if o.source_type in ['cfour', 'gaussian']:

            from CQCParse.parsing import parse_from_source
            from dataclasses import asdict

            """
            asdict(o) {'source_type': 'gaussian', 
                        'lvl_theory': 'B3LYP', 'basis_set': 'cc-pVQZ', 
                        'base_file_loc': '/home/vlev/monorepo/src/../data_for_tests/g16_formaldehyde_B3LYPcc_pVQZ.out'}
            """
            these_results_dict = parse_from_source(requested_data=origin_to_req_data[o], 
                                                   reindex_modes=reindex_modes,  **asdict(o))
            
            dict_with_data.update(these_results_dict)

        elif o.source_type in ['wilson']:

            raise NotImplementedError('In-house retrieval not yet implemented')
            # FIXME: MR to implement:
            #   - use base_file_loc as a parsing starting point. Support one or
            #     more of:
            #     a) Formatted or unformatted numpy arrays delimited by header info
            #     b) Other systematic (e.g. OpenRSP format) tensors
            #     c) Header info plus locator for values in formats a) or b)

        else:
            raise ValueError('Unsupported source type for data obtainer')

    return dict_with_data


def save_obtained_data(dict_with_data: dict, format: str, filename: str = 'obtained_data_dict', save_to_dir=''):
    """
    Saving data wrapper function
    
    :param dict_with_data: return dict from wilson_data_obtainer
    :param format: options are json or pkl
    :param filename: filename with or without extention; if no extention - will be added as `.format` value
    """
    if format not in ['json', 'pkl']:
        raise NotImplementedError("Cannot save this file format")
    
    if '.' not in filename:
        filename += '.' + format
    
    if format == 'json':
        save_datadict_json(dict_with_data, filename, save_to_dir)
    elif format == 'pkl':
        save_datadict_pkl(dict_with_data, filename, save_to_dir)


def save_datadict_json(dict_with_data: dict, filename: str, save_to_dir):
    raise NotImplementedError("Need to be able to handle non-serializable dicts, and that's not implemented")

def save_datadict_pkl(dict_with_data: dict, filename: str, save_to_dir):
    from wilson_suite.wilson_utils.serialization import pickle_this_to
    pickle_this_to(dict_with_data, filename, save_to=save_to_dir)
