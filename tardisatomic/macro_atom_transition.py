import pandas as pd
import numpy as np

from astropy import constants

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
        self.macro_atom_transition_columns = ['atomic_number', 'source_ion_number',
                                         'source_level_number',
                                         'destination_ion_number',
                                         'destination_level_number',
                                         'transition_id']

        self.macro_atom_data = pd.DataFrame(
            columns=self.macro_atom_transition_columns)

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
        """
        Lines data with no MultiIndex
        """
        return self._lines

    @lines.setter
    def lines(self, value):
        self._lines = value.reset_index()


    @property
    def lines_level_number_lower(self):
        """
        Lines with MultiIndex set to atom, ion, level_number_lower
        """
        return self.lines.set_index(levels['atomic_number', 'ion_number',
                                     'level_number_lower'])


    @property
    def lines_level_number_lower_gby(self):
        """
        Lines grouped by MultiIndex
        """
        return self.lines.groupby(['atomic_number', 'ion_number',
                                          'level_number_lower'])



    @property
    def lines_level_number_upper(self):
        """
        Lines with MultiIndex set to atom, ion, level_number_upper
        """

        return self.lines.set_index(['atomic_number', 'ion_number',
                                     'level_number_upper'])


    @property
    def lines_level_number_upper_gby(self):
        """
        Lines grouped by MultiIndex
        """
        return self.lines.groupby(['atomic_number', 'ion_number',
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
        super(PEmissionDown, self).__init__(levels, lines, ionization_data)
        if 'metastable' not in self.levels:
            raise ValueError('metastable flag has not been set for levels yet')

    def create_transitions_db(self):
        """
        Calculating Database of all transitions
        """
        macro_atom_transitions = pd.DataFrame(
            columns=self.macro_atom_transition_columns)

        for index, data in self.levels.iterrows():
            temp_data = pd.DataFrame(columns=self.macro_atom_transition_columns)
            (energy, g, metastable, abs_energy) = data
            try:
                down_transitions = self.lines_level_number_upper_gby.get_group(
                    index)
            except KeyError:
                continue

            nu = down_transitions['nu'].values
            f_ul = down_transitions['f_ul'].values
            e_lower = self.levels['energy'].ix[index[:2]].ix[
                down_transitions['level_number_lower']].values
            p_coef = (nu**2 * f_ul / constants.c.cgs.value**2 *
                      (energy - e_lower))

            temp_data['p_coef'] = p_coef
            temp_data['transition_id'] = self.transition_id
            #temp_data['atomic_number'] =
            if len(down_transitions) > 3: 1/0


    def get_single_level_transitions(self, atomic_number, ion_number,
                                     level_number):

        """
        Get the transitions for a single level

        Parameters
        ----------

        atomic_number: int
            atomic number
        ion_number: int
            ion number (0 = not ionized, 1 = once ionized)
        level_number: int
            level number
        """

        index = atomic_number, ion_number, level_number
        (energy, g, metastable, abs_energy) = self.levels.ix[index]

        if metastable: return None

        try:
            down_transitions = self.lines_level_number_upper_gby.get_group(
                index)
        except KeyError:
            return None

        nu = down_transitions['nu'].values
        f_ul = down_transitions['f_ul'].values
        e_lower = self.levels['energy'].ix[index[:2]].ix[
            down_transitions['level_number_lower']].values

        p_coef = self._calculcate_transition_coefficients(nu, f_ul, e_lower,
                                                          energy)




    @staticmethod
    def _calculcate_transition_coefficients(nu, f_ul, e_lower, level_energy):
        """
        Calculate the transition coefficients


        """
        return nu**2 * f_ul / constants.c.cgs.value**2 * (level_energy - e_lower)

class PInternalDown(MacroAtomTransitions):
    """
    Class to calculate the p internal down transitions

    """

    transition_id = 1

    def __init__(self, levels, lines, ionization_data):
        super(PInternalDown, self).__init__(levels, lines, ionization_data)

