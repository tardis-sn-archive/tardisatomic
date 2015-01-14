import pandas as pd
import numpy as np
import scipy.special as scisp
from astropy import constants, units


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

    # ##
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
        return self._levels

    @levels.setter
    def levels(self, value):
        if value.index.names != ['atomic_number', 'ion_number', 'level_number']:
            raise ValueError('lines needs to have an index with '
                             '"atom, ion, level"')
        else:
            self._levels = value

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
    def __init__(self, levels, lines, ionization):
        super(PInternalDown, self).__init__(levels, lines, ionization)

    transition_id = 0

    def __call__(self, levels, lines):
        return pd.DataFrame(atomic, source_ion, source_level, destionation_ion, destination_level, transitionid, coef)


def p_internal_down(levels, lines):
    pass


class PcollisonalExcitation(MacroAtomTransitions):
    """
    Computes the van regemorter approximation on a temperature grid for the plasma array in TARDIS.


    """

    def __init__(self, levels, lines, ionization, ionization_cross_sections, T_grid):
        super(PcollisonalExcitation, self).__init__(levels, lines, ionization)
        self._cross_sections = ionization_cross_sections
        self._T_grid = T_grid


    def compute(self):
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

            self.macro_atom_data.ix[(atom, ion, level_number_lower, ion, level_number_upper)][
                'C_ul_conversion'] = C_ul_conversion
            for T, value in zip(self._T_grid, c_lu):
                column_name = "t%06d" % T
                self.macro_atom_data.ix[(atom, ion, level_number_lower, ion, level_number_upper)][column_name] = value


    def _compute_van_regemorter(T, f_lu, nu_lu):
        g = 0.2  # This value is set to 2. We should select the value based on the main quantum number
        u = constants.h.cgs.value * nu_lu / constants.k_B.cgs.value / T
        I = 13.6  # eV
        c0 = 5.46510e-11
        integ = 0.276 * np.exp(u) * scisp.exp1(u)
        gamma = (g, integ )
        c = c0 * T ** (0.5) * 14.5 * (I / constants.h.cgs.value / nu_lu ) * f_lu * constants.h.cgs.value * nu_lu \
            / constants.k_B.cgs.value / T * np.exp(- u)
        return c