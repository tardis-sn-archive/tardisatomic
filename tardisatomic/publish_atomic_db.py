# functions to publish atomic datasets
from glob import glob
import os

import pandas as pd
from tardis import atomic
import urlparse

from tardisatomic import util

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

    species_strings.append(','.join(current_line_list))
    return '\n'.join(species_strings)


def create_data_sources_overview(data_sources_dict):
    data_source_string = ''
    data_source_template = '**{0}:**\n\n'
    for data_source in ['kurucz', 'chianti', 'tardis_artificial_missing_ion']:
        if data_source not in data_sources_dict:
            continue

        data_source_string += data_source_template.format(data_source) + \
                              species_include_list2string(data_sources_dict[data_source]) + '\n\n'

    return data_source_string

def publish_atom_data_sets(file_pattern,
                           base_url='http://moria.astro.utoronto.ca/~wkerzend/tardis_atomic_databases/',
                           default_atomic_database='kurucz_cd23_chianti_H_He',
                           table_fname='current_public_table.rst',
                           publish_directory='published'):


    column_names = ['file name', 'uuid1', 'data sources', 'macroatom', 'zeta', 'synpp references', 'database version']

    if os.path.exists(publish_directory):
        if raw_input('directory for publishing exists ({0}). Okay to overwrite [y/N]'.format(publish_directory)).strip().lower() == 'y':
            os.system('rm -r {0}'.format(publish_directory))
        else:
            print "Aborting!"
            return
    os.system('mkdir {0}'.format(publish_directory))

    atom_data_table = {item:[] for item in column_names}
    index = []
    for fname in glob(file_pattern):
        atom_data = atomic.AtomData.from_hdf5(fname)
        base_fname = os.path.basename(fname).replace('.h5', '.zip')
        os.system('cp {src} {dest}; zip {basename} {dest}'.format(src=fname, dest=os.path.basename(fname),
                                                                  basename=os.path.join(publish_directory, base_fname)))

        index.append(base_fname.replace('.zip', ''))
        atom_data_table['file name'].append('`{0} <{1}>`_'.format(base_fname, urlparse.urljoin(base_url, base_fname)))
        atom_data_table['uuid1'].append(atom_data.uuid1)
        atom_data_table['data sources'].append(create_data_sources_overview(atom_data.data_sources))
        atom_data_table['database version'].append(atom_data.version)
        atom_data_table['macroatom'].append(atom_data.has_macro_atom)
        atom_data_table['zeta'].append(atom_data.has_zeta_data)
        atom_data_table['synpp references'].append(atom_data.has_synpp_refs)


    dataset_table = pd.DataFrame(atom_data_table, columns=column_names, index=index)

    index.remove(default_atomic_database)
    index = [default_atomic_database] + index

    dataset_table = dataset_table.ix[index]

    open(table_fname, 'w').write(util.make_table([dataset_table.columns] + dataset_table.values.tolist()))




