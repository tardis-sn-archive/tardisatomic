
import pytest

import numpy.testing as nptesting
import numpy as np

from tardisatomic.macro_atom_transition import (MacroAtomTransitions,
                                                PEmissionDown)


@pytest.fixture
def macro_atom_transitions(atomic_levels, atomic_lines, atomic_ionization_data):
    levels = atomic_levels.set_index(['atomic_number', 'ion_number',
                                      'level_number'])

    ionization_data = atomic_ionization_data.set_index(['atomic_number',
                                                        'ion_number'])
    return MacroAtomTransitions(levels, atomic_lines, ionization_data)

@pytest.fixture
def p_emission_down(atomic_levels, atomic_lines, atomic_ionization_data):
    levels = atomic_levels.set_index(['atomic_number', 'ion_number',
                                      'level_number'])

    ionization_data = atomic_ionization_data.set_index(['atomic_number',
                                                        'ion_number'])
    return PEmissionDown(levels, atomic_lines, ionization_data)

def test_create_macro_atom_structure1(macro_atom_transitions):
    assert len(macro_atom_transitions.macro_atom_data.columns) == 1
    assert macro_atom_transitions.macro_atom_data.columns[0] == 'transition_id'

def test_absolute_energy1(macro_atom_transitions):
    levels_abs_energy = macro_atom_transitions.levels['absolute_energy']
    assert np.all((levels_abs_energy.ix[2, 1].values -
                   macro_atom_transitions.ionization_data.ix[2, 1].values) >= 0.)

    assert np.any(levels_abs_energy.ix[2, 0] > 0)

def test_p_emission_down1(p_emission_down):
    p_emission_down.create_transitions_db()