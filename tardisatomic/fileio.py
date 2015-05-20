import re
import sqlite3
import numpy as np
import time
import os
import pandas as pd
from astropy import units

tardisatomic_data = os.path.join(os.path.dirname(__file__), 'data')
default_chianti_ionization_data = os.path.join(tardisatomic_data, 'chianti_ionization.dat')
default_nist_ionization_data = os.path.join(tardisatomic_data, 'nist_ionization.dat')


def read_chianti_ionization_data(fname=default_chianti_ionization_data):
    """
    Reading the CHIANTI ionization data
    """
    ionization_data = np.recfromtxt(fname, names=['atomic_number', 'ion_number', 'ionization_energy'], comments='%')
    ionization_data['ionization_energy'] **= -1

    ionization_data['ionization_energy'] = units.Unit('cm').to('eV', ionization_data['ionization_energy'], units.spectral())


def read_nist_ionization_data(fname=default_nist_ionization_data, full_information=False):
    """
    reading the nist data

    Parameters
    ----------

    full_information: bool
    if set to False the ground level information will be culled
    """

    energy_pattern = re.compile('[\(\[]?(\d+\.?\d+)(\(\d+\))?[\)\]]?')

    label_pattern = re.compile('\d\w(.*)')
    label_j_pattern = re.compile('\*?<(\d)/?(\d)?>')
    def convert_label2j(value):
        value = value.strip()
        if value == '1':
            return 0.0
        j_string = label_pattern.match(value).groups()[0].strip()
        try:
            return float(j_string)
        except:
            pass

        numerator, denominator = label_j_pattern.match(j_string).groups()
        j_value = float(numerator) / float(denominator or 1.)
        return j_value


    def convert_energy(value):
        value = value.strip()
        if value == 'Hydrogen':
            return 0.0

        return float(energy_pattern.match(value).groups()[0])

    ionization_data = np.recfromtxt(fname, skip_header=3, delimiter='|', skip_footer=1,
                                    usecols=(0, 2, 6, 8), names=['atomic_number', 'ion_number', 'ground_level_j', 'energy'],
                                    converters={'ground_level_j': convert_label2j, 'energy': convert_energy},
                                    dtype=(int, int, float, float))
    ionization_data = pd.DataFrame(ionization_data)
    ionization_data['ground_level_g'] = 2*ionization_data['ground_level_j'] + 1
    if full_information:
        return ionization_data.to_records(index=False)
    else:
        del ionization_data['ground_level_g']
        del ionization_data['ground_level_j']
        ionization_data.columns = ['atomic_number', 'ion_number', 'ionization_energy']
        return ionization_data.to_records(index=False)

