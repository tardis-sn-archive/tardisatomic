import pandas as pd
import numpy as np

class MacroAtomTransitions(object):
    """
    Base class for macro atom transitions

    Parameters:

    levels: ~pandas.DataFrame
        pandas dataframe with index on atomic_number, ion_number, level_number
    """
    _lines = None


    def __init__(self, levels, lines, ionization):
        macro_atom_transition_columns = ['atomic_number', 'source_ion_number',
                                         'source_level_number',
                                         'destination_ion_number',
                                         'destination_level_number',
                                         'transition_id']

        self.macro_atom_data = pd.DataFrame(
            columns=macro_atom_transition_columns)

        self.macro_atom_data.set_index(['atomic_number', 'source_ion_number',
                                        'source_level_number',
                                        'destination_ion_number',
                                        'destination_level_number'],
                                       inplace=True)
    ###
    # !!!! lines 2 indices for upper and lower !!! TODO
    ####
    def get_absolute_energy(self, atomic_number, ion_number, level_number):
        _t_levels = self.levels.reset_index()
        _t_ionization = self.ionization.reset_index()
        data_merged = _t_levels.merge(_t_ionization, on=['atomic_number', 'ion_number'], how='left')
        _energy_abs = data_merged['ionization_energy'] + data_merged['energy']
        _energy_abs = _energy_abs[np.isnan(_energy_abs)] = 0
        self.levels['energy_abs'] = _energy_abs

    @property
    def levels(self):
        return self._lines

    @levels.setter
    def levels(self, value):
        if value.index.names != ['atomic_number', 'ion_number', 'level_number']:
            raise ValueError('lines needs to have an index with '
                             '"atom, ion, level"')
        else:
            self._lines = value

    @property
    def lines(selfs):
        return selfs._lines

    @lines.setter
    def lines(self, value):
        if value.index.names != ['atomic_number', 'ion_number', 'line_id']:
            raise ValueError('lines needs to have an index with '
                             '"atom, ion, line"')
        else:
            self._lines = value


    @property
    def ionization(self):
        return self._ionization

    @ionization.setter
    def ionization(self, value):
        if value.index.names != ['atomic_number', 'ion_number']:
            raise ValueError('ionization needs to have an index with '
                             '"atom, ion"')
        else:
            self._ionization = value





class PInternalDown(MacroAtomTransitions):

    def __init__(self, levels, lines, bla):
        super(PInternalDown,self).__init__(levels, lines)

    transition_id = 0

    def __call__(self, levels, lines):
        return pd.DataFrame(atomic, source_ion, source_level, destionation_ion, destination_level, transitionid, coef)



def p_internal_down(levels, lines):
    pass
