import pandas as pd
import numpy as np

class MacroAtomTransitions(object):
    """
    Base class for macro atom transitions

    Parameters
    ----------

    levels: ~pandas.DataFrame
        pandas dataframe with index on atomic_number, ion_number, level_number

    lines: ~pandas.DataFrame
        pandas dataframe with all line_transitions

    ionization_data: ~pandas.DataFrame
        pandas dataframe with index on atomic_number, ion_number, level_number
    """

    def __init__(self, levels, lines, ionization_data):
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

        self.levels = levels
        self.lines = lines
        self.ionization_data = ionization_data


        self.levels['absolute_energy'] = self.get_absolute_energy()


    def get_absolute_energy(self):
        t_levels = self.levels.reset_index()
        t_ionization = self.ionization_data.reset_index()
        data_merged = t_levels.merge(t_ionization,
                                     on=['atomic_number', 'ion_number'],
                                     how='left')
        ionization_energy = data_merged['ionization_energy']
        ionization_energy[np.isnan(ionization_energy)] = 0.0
        energy_abs = ionization_energy + data_merged['energy']
        return energy_abs.values

    @property
    def levels(self):
        return self._levels

    @levels.setter
    def levels(self, value):
        if value.index.names != ['atomic_number', 'ion_number', 'level_number']:
            raise ValueError('lines needs to have an index with '
                             '"atom, ion, level"')
        else:
            self._levels = value

    @property
    def lines(self):
        return self._lines

    @lines.setter
    def lines(self, value):
        self._lines = value.reset_index()


    @property
    def lines_level_number_lower(self):
        """
        Lines with MultiIndex set to atom, ion, level_number_lower
        """
        return self.lines.set_index(['atomic_number', 'ion_number',
                                     'level_number_lower'])

    @property
    def lines_level_number_upper(self):
        """
        Lines with MultiIndex set to atom, ion, level_number_upper
        """

        return self.lines.set_index(['atomic_number', 'ion_number',
                                     'level_number_upper'])


    @property
    def ionization_data(self):
        return self._ionization

    @ionization_data.setter
    def ionization_data(self, value):
        if value.index.names != ['atomic_number', 'ion_number']:
            raise ValueError('ionization needs to have an index with '
                             '"atom, ion"')
        else:
            self._ionization = value




class PEmissionDown(MacroAtomTransitions):
    """
    Class to calculate the p emission down transitions

    """

    transition_id = 1

    def __init__(self, levels, lines, ionization_data):
        super(PInternalDown, self).__init__(levels, lines, ionization_data)



class PInternalDown(MacroAtomTransitions):
    """
    Class to calculate the p internal down transitions

    """

    transition_id = 1

    def __init__(self, levels, lines, ionization_data):
        super(PInternalDown, self).__init__(levels, lines, ionization_data)

