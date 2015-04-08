import logging
import sqlite3
import numpy as np
from astropy import constants, units
import h5py
import os

import pandas as pd


class BasicAtomicData(object):
    """
    Basic class for all data imports with TARDIS atomic.

    """

    __config = None

    def __init__(self, confg_dict):
        self.config = confg_dict

        try:
            self.load_data()
        except AttributeError:
            logging.critical('Error loading atomicdata. Object has no attribute load_data.')
            raise

        try:
            self.prepare()
        except:
            logging.warning(' Object has no attribute prepare.')

        try:
            self.filter()
        except:
            logging.warning(' Object has no attribute filter.')

        try:
            self.check()
        except:
            logging.warning(' Object has no attribute check.')

        try:
            self.save()
        except AttributeError:
            logging.critical('Error loading atomicdata. Object has no attribute load_data.')
            raise


    @property
    def config(self):
        return self.__config

    @config.setter
    def config(self, value):
        self.__config = value


    def __get_sql_conn(self):
        # get the SQL db
        try:
            sqlstr = self.config['SQL_ATOMIC_DB']
            sqlconn = sqlite3.connect(sqlstr)
            return sqlconn
        except KeyError:
            logging.critical('No atomic database given in the configuration')
            raise


class LoadLevelsAndLines(BasicAtomicData):
    def load_data(self):

        # get the SQL db

        self.sqlconn = self.__get_sql_conn()


        levels_sql_stmt = 'select atom, ion, level_id, energy, g, source from levels '
        try:
            where_stmt = ['atom <= %d' % self.config['MAX_ATOMIC_NUMBER']]
        except KeyError:
            where_stmt = ['atom <= 30']

        ion_filter_stmt = None
        atom_filter_stmt = None
        try:
            ion_filter_stmt = 'ion not in (%s)' % ', '.join(self.config['EXCLUDE_ION'])
        except KeyError:
            pass

        try:
            ion_filter_stmt = 'ion in (%s)' % ', '.join(self.config['INCLUDE_ION'])
        except KeyError:
            pass

        try:
            atom_filter_stmt = 'ion not in (%s)' % ', '.join(self.config['EXCLUDE_ATOM'])
        except KeyError:
            pass

        try:
            atom_filter_stmt = 'ion in (%s)' % ', '.join(self.config['INCLUDE_ATOM'])
        except KeyError:
            pass

        lines_data_sql_stmt = (
            'select id, wavelength, atom, ion, f_ul, f_lu, loggf, level_id_lower, level_id_upper, source from lines')
        lines_where_stmt = where_stmt.copy()
        if lines_where_stmt != []:
            lines_data_sql_stmt += ' where %s' % ' and '.join(lines_where_stmt)

        lines_data_sql_stmt += ' order by wavelength'

        lines_data_dtype = [('line_id', np.int), ('wavelength', np.float), ('atomic_number', np.int),
                            ('ion_number', np.int), ('f_ul', np.float),
                            ('f_lu', np.float), ('loggf', np.float), ('level_number_lower', np.int),
                            ('level_number_upper', np.int), ('source', '|S64')]

        lines_data = self.sqlconn.execute(lines_data_sql_stmt).fetchall()
        lines_data = np.array(lines_data, dtype=lines_data_dtype)

        lines_data = pd.DataFrame(lines_data)

        einstein_coeff = (4 * np.pi ** 2 * constants.e.gauss.value ** 2) / (
            constants.m_e.cgs.value * constants.c.cgs.value)

        lines_data['nu'] = units.Unit('angstrom').to('Hz', lines_data['wavelength'], units.spectral())
        lines_data['B_lu'] = einstein_coeff * lines_data['f_lu'] / (constants.h.cgs.value * lines_data['nu'])
        lines_data['B_ul'] = einstein_coeff * lines_data['f_ul'] / (constants.h.cgs.value * lines_data['nu'])
        lines_data['A_ul'] = 2 * einstein_coeff * lines_data['nu'] ** 2 / constants.c.cgs.value ** 2 * lines_data[
            'f_ul']

        lines_data_level_number_upper = lines_data.set_index(['atomic_number', 'ion_number', 'level_number_upper'])
        lines_data_level_number_upper = lines_data_level_number_upper.groupby(
            level=['atomic_number', 'ion_number', 'level_number_upper'])

        lines_data_level_number_lower = lines_data.set_index(['atomic_number', 'ion_number', 'level_number_lower'])
        lines_data_level_number_lower = lines_data_level_number_lower.groupby(
            level=['atomic_number', 'ion_number', 'level_number_lower'])

        if ion_filter_stmt is not None:
            where_stmt.append(ion_filter_stmt)
        if atom_filter_stmt is not None:
            where_stmt.append(atom_filter_stmt)

        if where_stmt != []:
            levels_sql_stmt += ' where %s' % ' and '.join(where_stmt)

        levels_sql_stmt += ' ORDER BY atom,ion, level_id'

        levels_dtype = [('atomic_number', np.int), ('ion_number', np.int), ('level_number', np.int),
                        ('energy', np.float), ('g', np.int), ('source', '|S64')]

        levels_data = self.sqlconn.execute(levels_sql_stmt).fetchall()
        levels_data = np.array(levels_data, dtype=levels_dtype)
        levels_data = pd.DataFrame(levels_data)
        levels_data.set_index(['atomic_number', 'ion_number', 'level_number'], inplace=True)

        try:
            lines_data_level_number_upper_meta = lines_data[
                lines_data.loggf > self.config['METASTABLE_LOGGF_THRESHOLD']]
        except KeyError:
            logging.warning('METASTABLE_LOGGF_THRESHOLD not set. Using default! loggf > -3 ')
            lines_data_level_number_upper_meta = lines_data[lines_data.loggf > -3]

        lines_data_level_number_upper_meta = lines_data_level_number_upper_meta.groupby(
            ['atomic_number', 'ion_number', 'level_number_upper'])

        count_down = lines_data_level_number_upper_meta['line_id'].count()
        levels_data['metastable'] = count_down == 0
        levels_data['metastable'] = pd.isnull(levels_data['metastable'])

        try:
            ionization_data = self.config['IONIZATION_DATA']
        except KeyError:
            logging.critical('No ionization data given in the configuration. To compute level ionization energy the '
                             'ionization data are required. Abort!')
            raise

        levels_atom_ion_index = zip(levels_data.index.get_level_values('atomic_number'),
                                    levels_data.index.get_level_values('ion_number').values + 1)

        levels_ionization_energy = ionization_data.ix[levels_atom_ion_index]
        levels_ionization_energy = levels_ionization_energy.values.reshape(levels_ionization_energy.shape[0])

        auto_ionizing_levels_mask = (levels_data.energy >= levels_ionization_energy).values

        auto_ionizing_levels = levels_data[auto_ionizing_levels_mask]

        auto_ionizing_line_ids = []

        # Culling lines_data with low gf values
        lines_data = lines_data[(lines_data.source != 'kurucz') | (lines_data.loggf > args.kurucz_loggf_threshold)]
        print "Cleaning auto-ionizing"
        for index in auto_ionizing_levels.index:
            index_in_lines = False
            try:
                line_group = lines_data_level_number_upper.get_group(index)
            except KeyError:
                pass
            else:
                auto_ionizing_line_ids += line_group.line_id.values.tolist()
                index_in_lines = True

            try:
                line_group = lines_data_level_number_lower.get_group(index)
            except KeyError:
                pass
            else:
                auto_ionizing_line_ids += line_group.line_id.values.tolist()
                index_in_lines = True

            if not index_in_lines:
                print "INDEX %s not in lines" % (index,)

        auto_ionizing_line_ids = np.unique(auto_ionizing_line_ids)
        lines_data = lines_data[~lines_data.reset_index().line_id.isin(auto_ionizing_line_ids).values]
        levels_data = levels_data[~auto_ionizing_levels_mask]

        print "cleaning levels which don't exist in lines"
        existing_levels_ids = []
        for idx, line in levels_data.reset_index().iterrows():
            index_in_lines = False
            index = line['atomic_number'], line['ion_number'], line['level_number']

            try:
                line_group = lines_data_level_number_upper.get_group(index)
            except KeyError:
                if index[0] == index[1]:
                    #print "Found fully ionized %s" % (index,)
                    index_in_lines = True
                elif line['source'].startswith('tardis_artificial'):
                    #print "Found artifical ion level %s" % (index, )
                    index_in_lines = True
                else:
                    pass
            else:
                index_in_lines = True

            try:
                line_group = lines_data_level_number_lower.get_group(index)
            except KeyError:
                if index[0] == index[1]:
                    #print "Found fully ionized %s" % (index,)
                    index_in_lines = True
                elif line['source'].startswith('tardis_artificial'):
                    #print "Found artifical ion level %s" % (index, )
                    index_in_lines = True
                else:
                    pass
                pass
            else:
                index_in_lines = True

            if not index_in_lines:
                pass
                #        print "INDEX %s not in lines" % (index,)
            else:
                existing_levels_ids.append(index)

        print "Found %d levels which don't exist in lines (total=%d levels)" % (
        len(levels_data) - len(existing_levels_ids), len(levels_data))
        levels_data = levels_data.ix[existing_levels_ids]
        levels_prev_index = pd.Series(index=levels_data.index.copy())
        levels_data.reset_index(inplace=True)
        levels_data.set_index(['atomic_number', 'ion_number'], inplace=True)

        print "Relabeling index numbers"

        for index in np.unique(levels_data.index.values):
            levels_prev_index.ix[index] = np.arange(levels_data.ix[index].level_number.count())
            levels_data.loc[index, 'level_number'] = np.arange(levels_data.ix[index].level_number.count())

        self._levels_data = levels_data
        self._levels_prev_index = levels_prev_index
        self._lines_data = lines_data




    def prepare(self):


        lines_data_level_number_upper = self._lines_data.set_index(
            ['atomic_number', 'ion_number', 'level_number_upper'])
        lines_data_level_number_upper_index = lines_data_level_number_upper.index.copy()
        lines_data_level_number_upper = lines_data_level_number_upper.groupby(
            level=['atomic_number', 'ion_number', 'level_number_upper'])

        lines_data_level_number_lower = self._lines_data.set_index(
            ['atomic_number', 'ion_number', 'level_number_lower'])
        lines_data_level_number_lower_index = lines_data_level_number_lower.index.copy()
        lines_data_level_number_lower = lines_data_level_number_lower.groupby(
            level=['atomic_number', 'ion_number', 'level_number_lower'])

        self._lines_data.level_number_lower = self._levels_prev_index.ix[lines_data_level_number_lower_index].values
        self._lines_data.level_number_upper = self._levels_prev_index.ix[lines_data_level_number_upper_index].values

        lines_data_level_number_upper = self._lines_data.set_index(
            ['atomic_number', 'ion_number', 'level_number_upper'])
        lines_data_level_number_upper_index = lines_data_level_number_upper.index.copy()
        lines_data_level_number_upper = lines_data_level_number_upper.groupby(
            level=['atomic_number', 'ion_number', 'level_number_upper'])

        lines_data_level_number_lower = self._lines_data.set_index(
            ['atomic_number', 'ion_number', 'level_number_lower'])
        lines_data_level_number_lower_index = lines_data_level_number_lower.index.copy()
        lines_data_level_number_lower = lines_data_level_number_lower.groupby(
            level=['atomic_number', 'ion_number', 'level_number_lower'])

        lines_data = self._lines_data.set_index('line_id')
        levels_data = self._levels_data.reset_index().set_index(['atomic_number', 'ion_number', 'level_number'])

        self.lines_data = lines_data
        self.levels_data = levels_data


    def save(self):
        with h5py.File(self.config['HDF5_FILE']) as hdf5_file:
            hdf5_file['lines_data'] = self.lines_data.reset_index().to_records(index=False)
            hdf5_file['levels_data'] = self.levels_data.reset_index().to_records(index=False)


class ZetaData(BasicAtomicData):
    def load_data(self):
        zeta_datafile = os.path.join(os.path.dirname(__file__), 'data', 'knox_long_recombination_zeta.dat')
        self.zeta_data = np.loadtxt(zeta_datafile, usecols=xrange(1, 23), dtype=np.float64)

    def save(self):
        with h5py.File(self.config['HDF5_FILE']) as hdf5_file:
            hdf5_file['zeta_data'] = self.zeta_data
            hdf5_file['zeta_data'].attrs['t_rad'] = np.arange(2000, 42000, 2000)
            hdf5_file['zeta_data'].attrs['source'] = 'Used with kind permission from Knox Long'


class MacroAtomData(BasicAtomicData):
    """
    Here comes the new macroatom

    """
    # ToDo: Wolfgang will add the new macroatom ;-)

    def load(self):
        pass


class CollisionData(BasicAtomicData):
    def load(self):
        self.sqlconn = self.__get_sql_conn()
        collision_data_exists = \
        self.sqlconn.execute('SELECT count(name) FROM sqlite_master WHERE name="collision_data"').fetchone()[0]
        if collision_data_exists == 0:
            logging.warning("WARNING: Collision data requested but not in the database")
            raise

        collision_columns = zip(*self.sqlconn.execute('PRAGMA table_info(collision_data)').fetchall())[1]
        temperature_columns = [item for item in collision_columns if item.startswith('t')]
        self._temperatures = [float(item[1:]) for item in temperature_columns]
        select_collision_stmt = "select atom, ion, level_number_upper, level_number_lower, g_ratio, delta_e, %s" \
                                " from collision_data"

        select_collision_stmt = select_collision_stmt % (', '.join(temperature_columns))

        collision_data_dtype = [('atomic_number', np.int), ('ion_number', np.int), ('level_number_upper', np.int),
                                ('level_number_lower', np.int), ('g_ratio', np.float), ('delta_e', np.float)]

        collision_data_dtype += [(str(item), np.float) for item in temperature_columns]

        collision_data = self.sqlconn.execute(select_collision_stmt).fetchall()
        self._collision_data = np.array(collision_data, dtype=collision_data_dtype)

    def save(self):
        with h5py.File(self.config['HDF5_FILE']) as hdf5_file:
            hdf5_file['collision_data'] = self._collision_data
            hdf5_file['collision_data'].attrs['temperatures'] = self._temperatures



