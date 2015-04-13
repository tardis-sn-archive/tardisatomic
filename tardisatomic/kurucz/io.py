import time
import re

import numpy as np

def read_gfall_raw(fname):
    start_time = time.time()
    #FORMAT(F11.4,F7.3,F6.2,F12.3,F5.2,1X,A10,F12.3,F5.2,1X,A10,
    #3F6.2,A4,2I2,I3,F6.3,I3,F6.3,2I5,1X,A1,A1,1X,A1,A1,i1,A3,2I5,I6)

    kurucz_fortran_format =('F11.4,F7.3,F6.2,F12.3,F5.2,1X,A10,F12.3,'
                            'F5.2,1X,A10,F6.2,F6.2,F6.2,A4,I2,I2,I3,F6.3,I3,F6.3,'
                            'I5,I5,1X,A1,A1,1X,A1,A1,I1,A3,I5,I5,I6')

    number_match = re.compile(r'\d+(\.\d+)?')
    type_match = re.compile(r'[FIXA]')
    type_dict = {'F':np.float64, 'I':np.int64, 'X':'S1', 'A':'S10'}
    field_types = tuple([type_dict[item] for item in number_match.sub(
        '', kurucz_fortran_format).split(',')])

    field_widths = type_match.sub('', kurucz_fortran_format)
    field_widths = map(int, re.sub(r'\.\d+', '', field_widths).split(','))

    gfall = np.genfromtxt(fname, dtype=field_types, delimiter=field_widths,
                          skiprows=2)

    print "took %.2f seconds" % (time.time() - start_time)
    return gfall

def gfall_raw_2_db(gfall_raw, conn):
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