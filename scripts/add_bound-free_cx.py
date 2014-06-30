__author__ = 'michi'
import h5py as h5
import numpy as np

fname_in = "kurucz_cd23_chianti_H_He.h5"
fname_out = "kurucz_cd23_chianti_H_He_bf.h5"
finput = 'bf_cx.dat'

fin = h5.File(fname_in,'r')
fout = h5.File(fname_out,'w')


for k in fin.keys():
    fout[k] = fin[k].value

levels_data = fin['levels_data'].value

cxinput = np.loadtxt(finput, comments="#", dtype=[('atom', '|S6'), ('ion', 'S6'), ('level', 'S6'), ('cx', np.float64)])

ion_cx_data = np.zeros(levels_data.shape, dtype=[('atomic_number', '<i8'), ('ion_number', '<i8'),
                                                 ('level_number', '<i8'),('cross_section', '<f8')])

def is_eq_or_astrix(v,c):
    if v =='*':
        return True
    elif int(v) == int(c):
        return True
    else:
        return False

for i, level in enumerate(levels_data):
    ion_cx_data[i][0] = level[0]
    ion_cx_data[i][1] = level[1]
    ion_cx_data[i][2] = level[2]
    for ii, cx in enumerate(cxinput):
        if is_eq_or_astrix(cx[0],level[0]):
            if is_eq_or_astrix(cx[1],level[1]):
                if is_eq_or_astrix(cx[2],level[2]):
                   ion_cx_data [i][3] = cx[3]



fout['ion_cx_data'] = ion_cx_data

cx_sp_data = np.zeros((1,), dtype=[('atomic_number', '<i8'), ('ion_number', '<i8'), ('level_number', '<i8'),
                                   ('wavelength', '<f8'), ('cross_section', '<f8')])
fout['ion_cx_sp_data'] = cx_sp_data
fin.close()
fout.close()
