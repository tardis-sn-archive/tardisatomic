import pandas as pd
import numpy as np
import scipy.special as scisp
from astropy import constants, units

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

    def compute(self):
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

def p_internal_down(levels, lines):
    pass


class PCollisonalExcitation(MacroAtomTransitions):
    """
    Computes the van regemorter approximation on a temperature grid for the plasma array in TARDIS.


    """

    def __init__(self, levels, lines, ionization, ionization_cross_sections, T_grid):
        super(PcollisonalExcitation, self).__init__(levels, lines, ionization)
        self._cross_sections = ionization_cross_sections
        self._T_grid = T_grid


    def compute(self):
        self.macro_atom_data.reset_index(inplace=True).set_index(['atomic_number','source_ion_number','source_level_number'])
        self.macro_atom_data.index = self._levels.index

        self.macro_atom_data ['C_ul_conversion'] = None
        for T in self._T_grid:
            column_name = "t%06d" % T
            self.macro_atom_data[column_name] = None


        for row in self._lines.reset_index().iterrows():
            row_data = row[1]
            atom = int(row_data['atomic_number'])
            ion = int(row_data['ion_number'])
            nu = units.Unit('angstrom').to('Hz', row_data['wavelength'], units.spectral())
            f_lu = row_data['f_lu']
            level_number_upper = int(row_data['level_number_upper'])
            level_number_lower = int(row_data['level_number_lower'])
            g_upper = self._levels.ix[(atom, ion, level_number_upper)]
            g_lower = self._levels.ix[(atom, ion, level_number_lower)]

            c_lu = self._compute_van_regemorter(self._T_grid, f_lu, nu)
            C_ul_conversion = g_upper / float(g_lower)

            self.macro_atom_data.loc[(atom, ion, level_number_lower)]['destination_level_number'] = level_number_upper
            self.macro_atom_data.loc[(atom, ion, level_number_lower)]['destination_level_number'] = ion
            self.macro_atom_data.loc[(atom, ion, level_number_lower)]['C_ul_conversion'] = C_ul_conversion
            for T, value in zip(self._T_grid, c_lu):
                column_name = "t%06d" % T
                self.macro_atom_data.loc[(atom, ion, level_number_lower)][column_name] = value


    def _compute_van_regemorter(self, T, f_lu, nu_lu):
        g = 0.2  # This value is set to 2. We should select the value based on the main quantum number
        u = constants.h.cgs.value * nu_lu / constants.k_B.cgs.value / T
        I = 13.6  # eV
        c0 = 5.46510e-11
        integ = 0.276 * np.exp(u) * scisp.exp1(u)
        gamma = (g, integ )
        c = c0 * T ** (0.5) * 14.5 * (I / constants.h.cgs.value / nu_lu ) * f_lu * constants.h.cgs.value * nu_lu \
            / constants.k_B.cgs.value / T * np.exp(- u)
        return c


class PCollisonalIonization(MacroAtomTransitions):
    def __init__(self, levels, lines, ionization, ionization_cross_sections, T_grid):
        super(PcollisonalIonization, self).__init__(levels, lines, ionization)
        self._cross_sections = ionization_cross_sections
        self._T_grid = T_grid


    def compute(self):
        for row in self._levels.reset_index().iterrows():

        pass

    def _compute_seaton(self, T, ion, sigma_th, nu_th, ):
        gi = (lambda x: 0.3 if x >=2 else ((lambda x: 0.2 if x==1 else 0.1)(x)))(ion)
        seaton_const = 1.55e13 # in CGS
        return  T**(-0.5) * seaton_const * gi * sigma_th / constants.h.cgs.value / nu_th / constants.k_B.cgs / T

