import requests
from bs4 import BeautifulSoup
from uncertainties import ufloat_fromstr

from StringIO import StringIO

from tardisatomic.alchemy import Ion

import pandas as pd
import numpy as np

def download_ionization(spectra='h-uuu'):
    data = {'spectra':spectra, 'format':0, 'at_num_out':True,
            'ion_charge_out':True, 'spectra':'h-uuu', 'format':1,
            'e_out':0, 'unc_out':True, 'units':1}
    print "Downloading ionization data from NIST"
    r = requests.post('http://physics.nist.gov/cgi-bin/ASD/ie.pl', data=data)

    print "Parsing "
    bs_ion = BeautifulSoup(r.text, 'html5')

    anchors = bs_ion.findAll('a')
    for a in anchors:
        if a.string is None:
            continue
        a.previousSibling.replaceWith(a.previousSibling + a.string)


    [s.extract() for s in bs_ion('a')]


    return pd.read_table(StringIO(bs_ion('pre')[0].get_text()), delimiter='|',
                  comment='-', skiprows=3, names=['atomic_number',
                                                    'ion_number',
                                                    'ionization_string'],
                  index_col=False)



def parse_ionization_data(ionization_df):
    ionization_data = ionization_df.copy()

    ionization_data['ionization_string'] = (
        ionization_data['ionization_string'].apply(lambda x: x.strip()))

    def parse_ionization_method(ionization_string):
        if ionization_string.startswith('('):
            return 'theoretical'
        elif ionization_string.startswith('['):
            return 'interpolation'
        else:
            return 'measured'

    ionization_data['method'] = ionization_data.ionization_string.apply(
        parse_ionization_method)

    def parse_ionization_string(ionization_string):
        if ionization_string == '':
            return None
        return ufloat_fromstr(ionization_string.strip('(').replace('))', ')')
                              .strip('[]'))
    ionization_values = ionization_data.ionization_string.apply(
        parse_ionization_string)

    ionization_data['ionization_energy'] = [
        np.nan if item is None else item.nominal_value
        for item in ionization_values]
    ionization_data['ionization_energy_uncertainty'] = [
        np.nan if item is None else item.std_dev
        for item in ionization_values]

    del ionization_data['ionization_string']
    return ionization_data

def ingest_ionization(ionization_data, atomic_db):
    for i, row in ionization_data:
        Ion()