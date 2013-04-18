# Importing certain species from the chiantidb into the Kurucz db

import h5py
import os
import numpy as np
from StringIO import StringIO
from astropy import table, constants, units
import chianti.core as ch
import pandas as pd
import pdb
from scipy import interpolate

kb_ev = constants.k_B.cgs.to('eV / K').value

basic_atom_data = h5py.File(os.path.join(os.path.dirname(__file__), 'data', 'atom_data_basic.h5'))['basic_atom_data']
symbol2z = dict(zip(basic_atom_data['symbol'], basic_atom_data['atomic_number']))

def read_chianti(symbol, ion_number, level_observed=True, temperatures = np.linspace(2000, 50000, 20)):
    ion_data = ch.ion('%s_%d' % (symbol.lower(), ion_number))
    levels_data = {}
    levels_data['level_number'] = ion_data.Elvlc['lvl']

    temperatures = np.array(temperatures)
    if level_observed:
        levels_data['energy'] = units.Unit('cm').to('eV', 1 / np.array(ion_data.Elvlc['ecm']), units.spectral())
    else:
        levels_data['energy'] = units.Unit('cm').to('eV', 1 / np.array(ion_data.Elvlc['ecmth']), units.spectral())
    levels_data['g'] = 2*np.array(ion_data.Elvlc['j']) + 1

    if levels_data['energy'][0] != 0.0:
        raise ValueError('Level 0 energy is not 0.0')

    levels_data = pd.DataFrame(levels_data)
    levels_data.set_index('level_number', inplace=True)

    last_bound_level = levels_data[levels_data['energy'] < ion_data.Ip].index[-1]

    levels_data = levels_data.ix[:last_bound_level]

    lines_data = {}

    lines_data['wavelength'] = ion_data.Wgfa['wvl']
    lines_data['level_number_lower'] = ion_data.Wgfa['lvl1']
    lines_data['level_number_upper'] = ion_data.Wgfa['lvl2']
    lines_data['A_ul'] = ion_data.Wgfa['avalue']
    nu = units.Unit('angstrom').to('Hz', lines_data['wavelength'], units.spectral())
    lines_data = pd.DataFrame(lines_data)

    g_lower = levels_data['g'].ix[lines_data['level_number_lower'].values].values
    g_upper = levels_data['g'].ix[lines_data['level_number_upper'].values].values

    A_coeff = (8 * np.pi**2 * constants.e.gauss.value**2 * nu**2)/ (constants.m_e.cgs.value * constants.c.cgs.value**3)
    lines_data['f_ul'] = lines_data['A_ul'] / A_coeff
    lines_data['f_lu'] = (lines_data['A_ul'] * g_upper) / (A_coeff * g_lower)
    lines_data['loggf'] = np.log10(lines_data['f_lu'] * g_lower)
    lines_data['wavelength'] = lines_data['wavelength'] / (1.0 + 2.735182E-4 + 131.4182 / lines_data['wavelength']**2
                                                           + 2.76249E8 / lines_data['wavelength']**4)

    lines_data = lines_data[lines_data['level_number_upper'] <= last_bound_level]
    collision_data_index = pd.MultiIndex.from_arrays((ion_data.Splups['lvl1'], ion_data.Splups['lvl2']))

    c_lvl1 = []
    c_lvl2 = []
    c_upper_lowers = []

    g_ratios = []
    delta_es = []


    for i, (lvl1, lvl2) in enumerate(zip(ion_data.Splups['lvl1'], ion_data.Splups['lvl2'])):
        if lvl2 > last_bound_level:
            continue

        c_lvl1.append(lvl1)
        c_lvl2.append(lvl2)
        c_upper_lower, g_ratio, delta_e = calculate_collisional_strength(ion_data.Splups, i, temperatures, levels_data)
        c_upper_lowers.append(c_upper_lower)
        g_ratios.append(g_ratio)
        delta_es.append(delta_e)

    c_upper_lowers = np.array(c_upper_lowers)
    g_ratios = np.array(g_ratios)
    delta_es = np.array(delta_es)


    collision_data = pd.DataFrame(c_upper_lowers, index=collision_data_index)
    collision_data['g_ratio'] = g_ratios

    #CAREFUL!!! delta_e has already been divided by k!!
    collision_data['delta_e'] = delta_es


    collision_data['level_number_lower'] = c_lvl1
    collision_data['level_number_upper'] = c_lvl2



    return levels_data, lines_data, collision_data

def calculate_collisional_strength(splups_data, splups_idx, temperature, level_data):
    """
        Function to calculation upsilon from Burgess & Tully 1992 (TType 1 - 4; Eq. 23 - 38)
    """

    c = splups_data['cups'][splups_idx]
    x_knots = np.linspace(0, 1, splups_data['nspl'][splups_idx])
    y_knots = splups_data['splups'][splups_idx]

    level_number_lower = splups_data['lvl1'][splups_idx]
    level_number_upper = splups_data['lvl2'][splups_idx]

    ttype = splups_data['ttype'][splups_idx]
    if ttype > 5: ttype -= 5

    kt = kb_ev * temperature

    delta_E = level_data.ix[level_number_upper]['energy'] - level_data.ix[level_number_lower]['energy']
    g_lower = level_data.ix[level_number_lower]['g']
    g_upper = level_data.ix[level_number_upper]['g']

    spline_tck = interpolate.splrep(x_knots, y_knots)

    if ttype == 1:
        x = 1 - np.log(c) / np.log(kt/delta_E + c)
        y_func = interpolate.splev(x, spline_tck)
        upsilon = y_func * np.log(kt/delta_E + np.exp(1))

    elif ttype == 2:
        x = (kt/delta_E) / (kt/delta_E + c)
        y_func = interpolate.splev(x, spline_tck)
        upsilon = y_func

    elif ttype == 3:
        x = (kt/delta_E) / (kt/delta_E + c)
        y_func = interpolate.splev(x, spline_tck)
        upsilon = y_func / (kt/delta_E + 1)

    elif ttype == 4:
        x = 1 - np.log(c) / np.log(kt/delta_E + c)
        y_func = interpolate.splev(x, spline_tck)
        upsilon = y_func * np.log(kt/delta_E + c)

    elif ttype == 5:
        raise ValueError('Not sure what to do with ttype=5')

    #### 1992A&A...254..436B Equation 20 & 22 #####

    c_upper_lower = 8.63e-6 * upsilon  / (g_upper * temperature**.5)
    g_ratio = float(g_upper) / float(g_lower)
    delta_Ek = delta_E / kb_ev


    return c_upper_lower, g_ratio, delta_Ek


def insert_to_db(symbol, ion_number, conn, temperatures=None):
    atomic_number = int(symbol2z[symbol])



    curs = conn.cursor()


    curs.execute('delete from levels where atom=? and ion=?', (atomic_number, ion_number - 1))
    curs.execute('delete from lines where atom=? and ion=?', (atomic_number, ion_number - 1))


    collision_data_cols = curs.execute('pragma table_info(collision_data)').fetchall()


    if temperatures is None:
        print "Trying to infer temperatures from column names"

    if collision_data_cols == []:
        raise IOError('The given database doesn\'t contain a collision_data table - please create it')

    else:
        temperatures_data = [int(item.strip('t')) for item in zip(*collision_data_cols)[1] if item.startswith('t')]
        print "Inferred temperatures are %s" % (temperatures_data,)

    levels_data, lines_data, collision_data = read_chianti(symbol, ion_number, temperatures=temperatures_data)

    for key, line in lines_data.iterrows():
        curs.execute('insert into lines(wl, atom, ion, level_id_upper, level_id_lower, f_lu, f_ul, loggf, source) '
                     'values(?, ?, ?, ?, ?, ?, ?, ?, "chianti")',
                     (line['wavelength'], atomic_number, ion_number-1,
                      line['level_number_upper']-1, line['level_number_lower']-1,
                      line['f_lu'], line['f_ul'], line['loggf']))


    for key, level in levels_data.iterrows():
        count_down = curs.execute('select count(id) from lines where atom=? and ion=? and level_id_upper=?',
                     (atomic_number, ion_number-1, int(key-1))).fetchone()[0]

        curs.execute('insert into levels(atom, ion, energy, g, level_id, metastable, source) values(?, ?, ?, ?, ?, ?, "chianti")',
                     (atomic_number, ion_number-1, level['energy'], level['g'], int(key-1), count_down == 0))




    insert_stmt = 'insert into collision_data(source, atom, ion, level_number_lower, level_number_upper, %s, g_ratio, delta_e) values(%s)'

    insert_stmt = insert_stmt % (', '.join(['t%06d' % item for item in temperatures_data]), ','.join('?' * (len(temperatures_data ) + 7)))

    for (level_number_lower, level_number_upper), collision_data in collision_data.iterrows():
        c_ul =  list(collision_data[:len(temperatures_data)].values)
        g_ratio = collision_data['g_ratio']
        delta_e = collision_data['delta_e']
        level_number_lower = int(collision_data['level_number_lower'] - 1)
        level_number_upper = int(collision_data['level_number_upper'] - 1)

        collision_line_data = ["chianti", atomic_number, ion_number-1, level_number_lower, level_number_upper] + c_ul + [g_ratio, delta_e]


        curs.execute(insert_stmt, collision_line_data)

    conn.commit()









def create_collision_data_table(conn, temperatures=np.arange(2000, 50000, 2000)):

    curs = conn.cursor()


    collision_data_table_stmt = """create table collision_data(id integer primary key,
                                                source text,
                                                atom integer,
                                                ion integer,
                                                level_number_upper integer,
                                                level_number_lower integer,
                                                g_ratio float,
                                                delta_e float,
                                                %s)
                                                """
    temperatures = temperatures.astype(np.int64)
    temperature_fields = ',\n'.join(['t%06d float' % temperature for temperature in temperatures])

    collision_data_table_stmt = collision_data_table_stmt % temperature_fields


    curs.execute('drop table if exists collision_data')
    curs.execute(collision_data_table_stmt)

    print "Created table collision_data"








