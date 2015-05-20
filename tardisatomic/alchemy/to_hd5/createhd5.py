from abc import ABCMeta
import platform
import hashlib
import os
from time import gmtime, strftime

import numpy as np

from pandas import HDFStore
from tardis.atomic import plugins as atomic_plugins
from tardisatomic.alchemy import Ion, Atom, Level, Transition, TransitionType, \
    TransitionValue, TransitionValueType
from sqlalchemy.orm import aliased
import pandas as pd
import gnupg


class BaseAtomicDatabase(object):
    __metaclass__ = ABCMeta

    _atomic_db_sql = None

    def _to_data_frame(self, sql_data):
        rec_data = [rec.__dict__ for rec in sql_data]
        return pd.DataFrame.from_records(rec_data)

    # ToDO:Finish that

    @property
    def atomic_db(self):
        return self._atomic_db_sql

    @atomic_db.setter
    def atomic_db(self, value):
        self._atomic_db_sql = value


class Atoms(atomic_plugins.Atoms, BaseAtomicDatabase):
    def __init__(self, **kwargs):
        self.exclude_species = kwargs.get('exclude_species', [])

    def load_sql(self):
        atomic_data = self.atomic_db.session.query(Atom).order_by(
            Atom.atomic_number).all()
        self.data = self._to_data_frame(atomic_data)
        try:
            self.data.drop('_sa_instance_state', 1, inplace=True)
        except:
            "No drop in atomic!"

        # To avoid the unicode Problem in hdf and python 2.7
        self.data['name'] = self.data['name'].apply(lambda x: str(x))
        self.data['symbol'] = self.data['symbol'].apply(lambda x: str(x))



class Ions(atomic_plugins.Ions, BaseAtomicDatabase):
    def __init__(self, **kwargs):
        self.exclude_species = kwargs.get('exclude_species', [])
        self.max_ionization_energy = kwargs.get('max_ionization_energy', np.inf)

    def load_sql(self):
        ion_data = self.atomic_db.session.query(Ion, Atom).join("atom").values(
            Atom.atomic_number, Ion.ion_number,
            Ion.ionization_energy)
        self.data = self._to_data_frame(ion_data)
        try:
            self.data.drop('_labels', 1, inplace=True)
        except:
            "No drop in Ions!"

        a = 1


class Levels(atomic_plugins.Levels, BaseAtomicDatabase):
    def __init__(self, **kwargs):
        self.exclude_species = kwargs.get('exclude_species', [])
        self.max_ionization_energy = kwargs.get('max_ionization_energy', np.inf)

    def load_sql(self):
        level_data = self.atomic_db.session.query(Level, Ion, Atom).join(
            "ion", "atom").values(Atom.atomic_number, Ion.ion_number,
                                  Level.level_number, Level.g, Level.energy)
        self.data = self._to_data_frame(level_data)
        try:
            self.data.drop('_labels', 1, inplace=True)
        except:
            "No drop in Levels!"
        a = 1


# ToDo: Change the super class tin TARDIS!!
class Lines(atomic_plugins.Lines, BaseAtomicDatabase):
    def __init__(self, **kwargs):
        self.exclude_species = kwargs.get('exclude_species', [])
        self.max_ionization_energy = kwargs.get('max_ionization_energy', np.inf)
        self.transition_types = kwargs.get('transition_types', ['lines'])

    def load_sql(self):
        target_level = aliased(Level)
        source_level = aliased(Level)
        trans_a = aliased(Transition)
        trans_b = aliased(Transition)
        trans_data = \
            self.atomic_db.session.query(TransitionValue,
                                         TransitionValue.value,
                                         Transition,
                                         Transition.id.label(
                                             "transition_id"),
                                         target_level.id.label(
                                             "target_level_id"),
                                         source_level.id.label(
                                             "source_level_id"),
                                         target_level.level_number.label(
                                             "target_level_number"),
                                         source_level.level_number.label(
                                             "source_level_number"),
                                         source_level.g.label(
                                             'g_lower'),
                                         target_level.g.label(
                                             'g_upper'),
                                         TransitionType.name,
                                         TransitionValueType.name
                                         ).join(Transition).join(
                target_level,
                Transition.target_level_id == target_level.id).join(
                source_level,
                Transition.target_level_id == source_level.id
            ).join(
                TransitionType).filter(TransitionType.name == 'Line')

        raw_line_transitions = self._to_data_frame(trans_data)
        raw_line_transitions.drop('Transition', axis=1, inplace=True)
        raw_line_transitions.drop('TransitionValue', axis=1, inplace=True)
        raw_line_transitions_wl = raw_line_transitions.ix[raw_line_transitions[
                                                              'name'] ==
                                                          'wavelength']

        raw_line_transitions_wl.rename(columns={'value': 'wavelength'},
                                       inplace=True)
        raw_line_transitions_loggf = raw_line_transitions.ix[
            raw_line_transitions[
                'name'] == 'loggf']
        raw_line_transitions_loggf.rename(columns={'value': 'loggf'},
                                          inplace=True)
        self.line_transitions = pd.merge(raw_line_transitions_loggf[[
            'transition_id', 'loggf']],
                                         raw_line_transitions_wl,
                                         on='transition_id')

        self.line_transitions['f_lu'] = pd.Series(np.zeros(
            len(self.line_transitions)))
        self.line_transitions['f_ul'] = pd.Series(np.zeros(
            len(self.line_transitions)))
        self.line_transitions['f_ul'] = np.power(10, self.line_transitions[
            'loggf'] \
                                                 .astype('float').div(
            self.line_transitions['g_upper'].astype(
                'float')))

        self.line_transitions['f_lu'] = np.power(10, self.line_transitions[
            'loggf'] \
                                                 .astype('float').div(
            self.line_transitions['g_lower'].astype(
                'float')))

        self.data = self.line_transitions
        try:
            self.data.drop('_labels', 1, inplace=True)
        except:
            "No drop in Ions!"

        self.data['name'] = self.data['name'].apply(lambda x: str(x))
        self.data['loggf'] = self.data['loggf'].apply(lambda x: float(x))
        self.data['wavelength'] = self.data['wavelength'].apply(lambda x:
                                                                float(x))
        a = 1

class CreateHDF(object):
    atomic_database_version = 'v2'

    def __init__(self, atom_db, hdf_file, data_types, **kwargs):
        """
        Creates the HDF file for TARDIS from the atomic database.

        :param atom_db: The AtomicDatabase object from
        tardisatomic.base containing all atomic data.

        :param exclude_species: Specifies the specie by atomic number which
        are exclude from the atomic data set. By default all species are
        included.

        :param max_ionization_level: Specifies the highest ionization level
        included in the atomic data set. All higher ionization levels are
        excluded. By default all ionization levels are included.

        :param transition_types: Specifies the transition types included in
        the atomic data set. Default: line
        """

        self.atomic_db = atom_db
        self.exclude_species = kwargs.get('exclude_species', None)
        self.max_ionization_level = kwargs.get('max_ionization_level', None)
        self.transition_types = kwargs.get('transition_types', None)

        self.hdf_buf = self._open_hdf_buf(hdf_file)

        self.data_type_dict = self._create_atomic_data_type_dict()

        self.data_type_dict = self._create_atomic_data_type_dict()

        self.atomic_plugins = self._create_atomic_objects(data_types,
                                                          self.data_type_dict)

        self._load_atomic_data(self.atomic_plugins, self.atomic_db)
        self._save_to_hdf(self.hdf_buf)

    def _open_hdf_buf(self, fname):
        if not os.path.isfile(fname):
            print "file does exist at this time. Create new HDF file."
            return pd.HDFStore(fname, 'w')
        else:
            print "Error: File hdf exist"
            raise

    @staticmethod
    def _create_atomic_data_type_dict():
        return {x.hdf_name: x for x in
                BaseAtomicDatabase.__subclasses__()}

    def _fix_unicode_column_names(self, df):
        column_name = [x.encode('ascii', 'ignore') for x in df.columns.values]
        df.columns = column_name

    def _save_to_hdf(self, hdf_buf):
        for plugin in self.atomic_plugins:
            try:
                self._fix_unicode_column_names(plugin.data)
                plugin.save_hdf(hdf_buf)
            except AttributeError, TypeError:
                print("Can't save {0} to HDF.".format(plugin.hdf_name))

    def _create_atomic_objects(self, atomic_data_types, data_type_dict,
                               **kwargs):
        data_types = []

        for key in atomic_data_types:
            data_types.append(data_type_dict[key](**kwargs))

        return data_types

    def _load_atomic_data(self, atomic_plugins, atomic_db):

        for plugin in atomic_plugins:
            print("Loading {0} from atom db.".format(plugin.hdf_name))
            plugin.atomic_db = atomic_db
            plugin.load_sql()

    def close_hdf(self, anonymous=False):
        if not anonymous:
            uname = platform.uname()
            user = os.getlogin()
        else:
            uname = 'anonymous'
            user = 'anonymous'

        ctime = strftime("%Y-%m-%d %H:%M:%S", gmtime())
        md5_hash = hashlib.md5()

        gpg = gnupg.GPG()
        gpg.encoding = 'utf-8'

        self.hdf_buf.put('metadata/system', pd.Series(uname))
        self.hdf_buf.put('metadata/username', pd.Series(user))
        self.hdf_buf.put('metadata/creation_time', pd.Series(ctime))

        for key  in self.hdf_buf.keys():
            md5_hash.update(str(self.hdf_buf[key].values))

        data_to_sign = "TARDISATOMIC\n" "creation_time: " + ctime + "\n" + \
                       "System:\n" + \
                       str(uname) + "\n" + "User: " \
                                      "" + \
                       user + "\n" + "md5: " + md5_hash.hexdigest() + "\n"
        signed_data = gpg.sign(data_to_sign)

        self.hdf_buf.put('metadata/gpg', pd.Series(signed_data.data))
        self.hdf_buf.put('metadata/md5', pd.Series(md5_hash.hexdigest()))
        self.hdf_buf.close()
