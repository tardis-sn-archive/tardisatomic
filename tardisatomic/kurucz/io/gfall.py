import time
import re

import numpy as np
import pandas as pd

from astropy import units as u

def read_gfall_raw(fname):
    """
    Reading in a normal gfall.dat (please remove any empty lines)

    Parameters
    ----------

    fname: ~str
        path to gfall.dat

    Returns
    -------
        : pandas.DataFrame
            pandas Dataframe represenation of gfall
    """
    start_time = time.time()
    #FORMAT(F11.4,F7.3,F6.2,F12.3,F5.2,1X,A10,F12.3,F5.2,1X,A10,
    #3F6.2,A4,2I2,I3,F6.3,I3,F6.3,2I5,1X,A1,A1,1X,A1,A1,i1,A3,2I5,I6)

    kurucz_fortran_format =('F11.4,F7.3,F6.2,F12.3,F5.2,1X,A10,F12.3,F5.2,1X,'
                            'A10,F6.2,F6.2,F6.2,A4,I2,I2,I3,F6.3,I3,F6.3,I5,I5,'
                            '1X,I1,A1,1X,I1,A1,I1,A3,I5,I5,I6')

    number_match = re.compile(r'\d+(\.\d+)?')
    type_match = re.compile(r'[FIXA]')
    type_dict = {'F':np.float64, 'I':np.int64, 'X':'S1', 'A':'S10'}
    field_types = tuple([type_dict[item] for item in number_match.sub(
        '', kurucz_fortran_format).split(',')])

    field_widths = type_match.sub('', kurucz_fortran_format)
    field_widths = map(int, re.sub(r'\.\d+', '', field_widths).split(','))

    gfall = np.genfromtxt(fname, dtype=field_types, delimiter=field_widths,
                          skiprows=2)
    columns = ['wavelength', 'loggf', 'element_code', 'e_first', 'j_first',
               'blank1', 'label_first', 'e_second', 'j_second', 'blank2',
               'label_second', 'log_gamma_rad', 'log_gamma_stark',
               'log_gamma_vderwaals', 'ref', 'nlte_level_no_first',
               'nlte_level_no_second', 'isotope', 'log_f_hyperfine',
               'isotope2', 'log_iso_abundance', 'hyper_shift_first',
               'hyper_shift_second', 'blank3', 'hyperfine_f_first',
               'hyperfine_note_first', 'blank4', 'hyperfine_f_second',
               'hyperfine_note_second', 'line_strength_class', 'line_code',
               'lande_g_first', 'lande_g_second', 'isotopic_shift']

    gfall = pd.DataFrame(gfall)
    gfall.columns = columns
    print "took {0:.2f} seconds".format(time.time() - start_time)
    return gfall



def parse_gfall(gfall_df):
    """
    GFall pandas dataframe from read_gfall
    :param gfall_df:
    :return:
    """

    gfall_df = gfall_df.copy()

    double_columns = [item.replace('_first', '') for item in gfall_df.columns if
                     item.endswith('first')]

    # due to the fact that energy is stored in 1/cm
    order_lower_upper = (gfall_df.e_first.abs() <
                         gfall_df.e_second.abs())


    for column in double_columns:
        data = pd.concat([gfall_df['{0}_first'.format(
            column)][order_lower_upper], gfall_df['{0}_second'.format(
            column)][~order_lower_upper]])

        gfall_df['{0}_lower'.format(column)] = data

        data = pd.concat([gfall_df['{0}_first'.format(
            column)][~order_lower_upper], gfall_df['{0}_second'.format(
            column)][order_lower_upper]])

        gfall_df['{0}_upper'.format(column)] = data


        del gfall_df['{0}_first'.format(column)]
        del gfall_df['{0}_second'.format(column)]

    gfall_df.e_lower = (gfall_df.e_lower.values / u.cm).to(
        u.eV, u.spectral()) #convert to eV
    gfall_df.e_upper = (gfall_df.e_upper.values / u.cm).to(
        u.eV, u.spectral()) #convert to eV




    gfall_df.wavelength *= 10

    gfall_df.label_lower = gfall_df.label_lower.apply(
        lambda x: x.strip())
    gfall_df.label_upper = gfall_df.label_upper.apply(
        lambda x: x.strip())

    gfall_df['e_lower_predicted'] = gfall_df.e_lower < 0
    gfall_df.e_lower = gfall_df.e_lower.abs()
    gfall_df['e_upper_predicted'] = gfall_df.e_upper < 0
    gfall_df.e_upper = gfall_df.e_upper.abs()

    gfall_df['g_lower'] = 2*gfall_df.j_lower + 1
    gfall_df['g_upper'] = 2*gfall_df.j_upper + 1

    gfall_df['atomic_number'] = gfall_df.element_code.astype(int)
    gfall_df['ion_number'] = (
        (gfall_df.element_code.values -
         gfall_df.atomic_number.values) * 100).round().astype(int)

    del gfall_df['element_code']

    return gfall_df

def extract_levels(gfall_df, selected_columns=None):
    """
    Extract the levels from the gfall dataframe

    Parameters
    ----------

    gfall_df: ~pandas.DataFrame
    selected_columns: list
        list of which columns to select (optional - default=None which selects
        a default set of columns)

    Returns
    -------
        : ~pandas.DataFrame
            a level DataFrame
    """

    if selected_columns is None:
        selected_columns = ['atomic_number', 'ion_number', 'energy', 'g',
                            'label', 'theoretical']


    if 'e_lower' not in gfall_df.columns:
        raise ValueError('gfall dataframe needs to be parsed before this '
                         'function can be used')

    column_renames = {'e_{0}':'energy', 'g_{0}':'g', 'label_{0}':'label',
                      'e_{0}_predicted':'theoretical'}


    e_lower_levels = gfall_df.rename(
        columns=dict([(key.format('lower'), value)
                      for key, value in column_renames.items()]))



    e_upper_levels = gfall_df.rename(
        columns=dict([(key.format('upper'), value)
                      for key, value in column_renames.items()]))

    levels = pd.concat([e_lower_levels[selected_columns],
                        e_upper_levels[selected_columns]])


    levels = levels.drop_duplicates(['atomic_number', 'ion_number',
                                             'energy', 'g', 'label']).sort(
        ['atomic_number', 'ion_number', 'energy'])

    levels_clean = levels.drop_duplicates(['atomic_number', 'ion_number',
                                           'energy'])

    levels_clean_level_number = levels_clean.groupby(
        ['atomic_number', 'ion_number']).g.transform(
        lambda x: np.arange(len(x))).values
    levels_clean['level_number'] = levels_clean_level_number

    levels_clean = levels_clean.set_index(
        ['atomic_number', 'ion_number', 'energy'])

    aie_index = [tuple(item)
                 for item in levels[
            ['atomic_number', 'ion_number', 'energy']].values.tolist()]

    level_number = levels_clean.level_number.loc[aie_index]

    levels['level_number'] = level_number.values
    levels['level_id'] = np.arange(len(levels))

    return levels
    

def extract_lines(gfall_df, levels_df, selected_columns=None):

    if selected_columns is None:
        selected_columns = ['wavelength','loggf' , 'atomic_number', 'ion_number']



    if 'e_lower' not in gfall_df.columns:
        raise ValueError('gfall dataframe needs to be parsed before this '
                         'function can be used')

    levels_df_idx = levels_df.set_index(['atomic_number', 'ion_number',
                                         'energy', 'g', 'label'])

    lines = gfall_df[selected_columns].copy()

    level_lower_idx = gfall_df[['atomic_number', 'ion_number', 'e_lower',
                                'g_lower', 'label_lower']].values.tolist()
    level_lower_idx = [tuple(item) for item in level_lower_idx]

    level_upper_idx = gfall_df[['atomic_number', 'ion_number', 'e_upper',
                                'g_upper', 'label_upper']].values.tolist()
    level_upper_idx = [tuple(item) for item in level_upper_idx]

    lines['level_id_lower'] = levels_df_idx.level_id.loc[level_lower_idx].values
    lines['level_id_upper'] = levels_df_idx.level_id.loc[level_upper_idx].values

    lines['level_number_lower'] = levels_df_idx.level_number.loc[
        level_lower_idx].values
    lines['level_number_upper'] = levels_df_idx.level_number.loc[
        level_upper_idx].values

    return lines


def ingest_gfall(levels, lines, atomic_db):
    pass