import os

import pytest
import tardisatomic
from tardisatomic.alchemy import Atom, Ion
from tardisatomic.base import AtomicDatabase
from sqlalchemy import inspect
from tardisatomic.nist.io import download_ionization, parse_ionization_data
from numpy import testing
from uncertainties import ufloat_fromstr
import re
from tardisatomic.alchemy.ingest.nist import NISTIonization
from sqlalchemy import and_

data_path = os.path.join(tardisatomic.__path__[0], 'tests', 'data')
test_db = os.path.join(data_path, 'example_tiny.db')
db_driver = 'sqlite:///' + test_db


@pytest.fixture
def atomic_database_empty():
    return AtomicDatabase('sqlite:///example_empty.db')


@pytest.fixture
def atomic_database():
    return AtomicDatabase(db_driver)


@pytest.fixture
def atomic_database_inspect(atomic_database):
    return inspect(atomic_database.engine)


@pytest.fixture
def raw_ions():
    return download_ionization('h-he')


def test_basic_db(atomic_database, atomic_database_inspect):
    tables = [u'atoms',
              u'data_sources',
              u'decay_types',
              u'decays',
              u'dtypes',
              u'ions',
              u'isotopes',
              u'levels',
              u'transition_types',
              u'transition_value_types',
              u'transition_values',
              u'transitions',
              u'units']

    assert atomic_database_inspect.get_table_names() == tables


def test_atoms_db(atomic_database):
    h1 = atomic_database.session.query(Atom).order_by(Atom.id).first()
    assert h1.id == 1
    assert h1.name == 'Hydrogen'
    assert atomic_database.session.query(Atom).count() == 118


def test_ions(raw_ions, atomic_database):
    ionization_string = \
    raw_ions.set_index(['atomic_number', 'ion_number']).ix[(1,
                                                            0)][
        'ionization_string'].strip()
    energy = float(re.findall("\d+.\d+", ionization_string)[0])

    testing.assert_almost_equal(energy, 13.5984, decimal=3)
    ionization_data = parse_ionization_data(raw_ions)
    assert list(ionization_data.columns._array_values()) == [u'atomic_number',
                                                             u'ion_number',
                                                             u'method',
                                                             u'ionization_energy',
                                                             u'ionization_energy_uncertainty']
    ion_engery = ionization_data.set_index(['atomic_number',
                                            'ion_number']).ix[(1,
                                                               0)][
        'ionization_energy']
    testing.assert_almost_equal(ion_engery)

    NISTIonization(atomic_database).ingest()
    atomic_databasesession.query(Ion).filter(and_(Ion.ion_number == 1, Atom.id \
                                                  == 1))


def test_levels():
    pass
