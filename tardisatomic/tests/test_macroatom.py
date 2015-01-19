
import pytest

import numpy.testing as nptesting
import numpy as np

from tardisatomic.macro_atom_transition import (MacroAtomTransitions,
                                                PEmissionDown, PInternalDown,
                                                PInternalUp)


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


@pytest.fixture
def p_internal_down(atomic_levels, atomic_lines, atomic_ionization_data):
    levels = atomic_levels.set_index(['atomic_number', 'ion_number',
                                      'level_number'])

    ionization_data = atomic_ionization_data.set_index(['atomic_number',
                                                        'ion_number'])
    return PInternalDown(levels, atomic_lines, ionization_data)


@pytest.fixture
def p_internal_up(atomic_levels, atomic_lines, atomic_ionization_data):
    levels = atomic_levels.set_index(['atomic_number', 'ion_number',
                                      'level_number'])

    ionization_data = atomic_ionization_data.set_index(['atomic_number',
                                                        'ion_number'])
    return PInternalUp(levels, atomic_lines, ionization_data)




def test_create_macro_atom_structure1(macro_atom_transitions):
    assert len(macro_atom_transitions.macro_atom_data.columns) == 1
    assert macro_atom_transitions.macro_atom_data.columns[0] == 'transition_id'

def test_absolute_energy1(macro_atom_transitions):
    levels_abs_energy = macro_atom_transitions.levels['absolute_energy']
    assert np.all((levels_abs_energy.ix[2, 1].values -
                   macro_atom_transitions.ionization_data.ix[2, 1].values) >= 0.)

    assert np.any(levels_abs_energy.ix[2, 0] > 0)

def test_p_emission_down1(p_emission_down):
    p_emission_down_group = p_emission_down._compute_transition_group(1, 0, 9)
    assert not np.any(np.isnan(p_emission_down_group.values.astype(np.float64)))
    assert np.all(p_emission_down_group['source_level_number'][0] ==
                  p_emission_down_group['source_level_number'].values)
    assert np.all(p_emission_down_group['destination_level_number'][0] !=
                  p_emission_down_group['destination_level_number'].values[1:])

    assert np.all(p_emission_down_group['destination_level_number'].values[1:] !=
                  p_emission_down_group['source_level_number'][0])




def test_p_internal_down1(p_internal_down):
    p_internal_down_group = p_internal_down._compute_transition_group(1, 0, 9)
    assert not np.any(np.isnan(p_internal_down_group.values.astype(np.float64)))
    assert np.all(p_internal_down_group['source_level_number'][0] ==
                  p_internal_down_group['source_level_number'].values)
    assert np.all(p_internal_down_group['destination_level_number'][0] !=
                  p_internal_down_group['destination_level_number'].values[1:])

    assert np.all(p_internal_down_group['destination_level_number'].values[1:] !=
                  p_internal_down_group['source_level_number'][0])

    p_internal_down_group = p_internal_down._compute_transition_group(1, 0, 2)

    nptesting.assert_allclose(p_internal_down_group['p_coef'].values, 0.0)



def test_p_internal_up1(p_internal_up):
    p_internal_up_group = p_internal_up._compute_transition_group(1, 0, 9)
    assert not np.any(np.isnan(p_internal_up_group.values.astype(np.float64)))
    assert np.all(p_internal_up_group['source_level_number'][0] ==
                  p_internal_up_group['source_level_number'].values)
    assert np.all(p_internal_up_group['destination_level_number'][0] !=
                  p_internal_up_group['destination_level_number'].values[1:])

    assert np.all(p_internal_up_group['destination_level_number'].values[1:] !=
                  p_internal_up_group['source_level_number'][0])



