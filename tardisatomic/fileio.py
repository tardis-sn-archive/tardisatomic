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


def readGFALLRaw(fname):
    start_time = time.time()
    #FORMAT(F11.4,F7.3,F6.2,F12.3,F5.2,1X,A10,F12.3,F5.2,1X,A10,
    #3F6.2,A4,2I2,I3,F6.3,I3,F6.3,2I5,1X,A1,A1,1X,A1,A1,i1,A3,2I5,I6)
    kurucz_fortran_format =('F11.4,F7.3,F6.2,F12.3,F5.2,1X,A10,F12.3,'
                            'F5.2,1X,A10,F6.2,F6.2,F6.2,A4,I2,I2,I3,F6.3,I3,F6.3,'
                            'I5,I5,1X,A1,A1,1X,A1,A1,I1,A3,I5,I5,I6')
    number_match = re.compile(r'\d+(\.\d+)?')
    type_match = re.compile(r'[FIXA]')
    type_dict = {'F':np.float64, 'I':np.int64, 'X':'S1', 'A':'S10'}
    field_types = tuple([type_dict[item] for item in number_match.sub('', kurucz_fortran_format).split(',')])
    field_widths = type_match.sub('', kurucz_fortran_format)
    field_widths = map(int, re.sub(r'\.\d+', '', field_widths).split(','))

    gfall = np.genfromtxt(fname, dtype=field_types, delimiter=field_widths)
    print "took %.2f seconds" % (time.time() - start_time)
    return gfall

def GFAllRaw2DB(gfall_raw, conn):
    insert_stmt_template = """insert into gfall(
    wavelength ,
    loggf ,
    atomic_number ,
    ion_number ,
    e_upper , 
    e_lower , 
    j_upper , 
    j_lower , 
    label_upper , 
    label_lower ,
    log_gamma_rad , 
    log_gamma_stark , 
    log_gamma_vdw ,
    ref , 
    nlte_level_no_upper , 
    nlte_level_no_lower , 
    isotope_number ,
    log_f_hyperfine ,
    isotope_number_2 ,
    log_isotope_fraction ,
    hyper_shift_upper ,
    hyper_shift_lower ,
    hyperfine_f_upper ,
    hyperfine_f_upper_note ,
    hyperfine_f_lower ,
    hyperfine_f_lower_note ,
    line_strength_class ,
    line_code ,
    lande_g_upper ,
    lande_g_lower ,
    isotope_shift ,
    predicted) values(%s)""" % (','.join(32 * ['?']))

    for line in gfall_raw:
        if abs(line['f3']) > abs(line['f7']): first_upper = True
        else: first_upper = False

        if line['f7'] < 0:
            line['f7'] = abs(line['f7'])
            predicted = True
        else:
            predicted = False
        insert_data_id = [0, 1]

        if first_upper:
            insert_data_id += [3, 7, 4, 8, 6, 10]
        else:
            insert_data_id += [7, 3, 8, 4, 10, 6]
            #log gamma
        insert_data_id += [11, 12, 13, 14]

        #nlte level idx
        if first_upper:
            insert_data_id += [15, 16]
        else:
            insert_data_id += [16, 15]
            #iso, hyperfine
        insert_data_id += [17, 18, 19, 20]

        #hyperfine shifts,etc
        if first_upper:
            insert_data_id += [21, 22, 23, 24, 25, 26]
        else:
            insert_data_id += [22, 21, 25, 26, 23, 24]

        #codes
        insert_data_id += [27, 28]

        # lande g factor
        if first_upper:
            insert_data_id += [29, 30]
        else:
            insert_data_id += [30, 29]

        insert_data_id += [31]

        insert_data = [line[idx] for idx in insert_data_id]
        insert_data += [predicted]
        insert_data = [int(item) if type(item) == np.int64 else item for item in insert_data]
        atom = int(line[2])
        ion = int(round(100 * (line[2]-atom)))
        insert_data.insert(2, ion)
        insert_data.insert(2, atom)
        conn.execute(insert_stmt_template, tuple(insert_data))


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
        return ionization_data.to_records(index=False)

