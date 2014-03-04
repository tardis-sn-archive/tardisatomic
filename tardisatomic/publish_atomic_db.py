# functions to publish atomic datasets
from glob import glob
import os

import pandas as pd
from tardis import atomic
import urlparse

def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def species_include_list2string(species_include_list, max_char=20):
    """
    Format a list like ['Si I', 'Si II-V', ...] -> 'Si I, Si II-V\n....'
    """
    species_strings = []
    line_char_count = 0
    current_line_list = []
    for item in species_include_list:
        line_char_count += len(item)
        current_line_list.append(item)

        if line_char_count > max_char:
            species_strings.append(', '.join(current_line_list))
            current_line_list = []
            line_char_count = 0

    return '\n'.join(species_strings)
def create_data_sources_overview(data_sources_dict):
    data_source_string = ''
    data_source_template = '**{0}:**\n\n'
    for data_source in ['kurucz', 'chianti', 'tardis_artificial_missing_ion']:
        data_source_string += data_source_template.format(data_source) + \
                              species_include_list2string(data_sources_dict[data_source]) + '\n\n'

    return data_source_string

def create_atom_data_table_from_file(file_pattern, base_url=''):
    column_names = ['file name', 'uuid1', 'data sources', 'macroatom', 'zeta', 'synpp references', 'database version']

    atom_data_table = {item:[] for item in column_names}

    for fname in glob(file_pattern):
        atom_data = atomic.AtomData.from_hdf5(fname)
        base_fname = os.path.basename(fname)
        atom_data_table['file name'].append('`{0} <{1}>`_'.format(base_fname, urlparse.urljoin(base_url, base_fname)))
        atom_data_table['uuid1'].append(atom_data.uuid1)
        atom_data_table['data sources'].append(create_data_sources_overview(atom_data.data_sources))
        atom_data_table['database version'].append(atom_data.version)
        atom_data_table['macroatom'].append(atom_data.has_macro_atom)
        atom_data_table['zeta'].append(atom_data.has_zeta_data)
        atom_data_table['synpp references'].append(atom_data.has_synpp_refs)


    return pd.DataFrame(atom_data_table, columns=column_names)

