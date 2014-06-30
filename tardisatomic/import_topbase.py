import numpy as np
from astropy import table, constants, units
import pandas as pd
import pdb
import h5py
import os
import pickle
import  re

basic_atom_data = h5py.File(os.path.join(os.path.dirname(__file__), 'data', 'atom_data_basic.h5'))['basic_atom_data']
symbol2z = dict(zip(basic_atom_data['symbol'], basic_atom_data['atomic_number']))

class TopBaseData:

    def __init__(self):
        self.level_data = None
        self.f_data = None
        self.phot_crosssections_data = None
        self.index_dict = {}

    def load_data(self, fvalue_fname, level_fname, photocx_fname):
        self._read_topbase_fvalues(fvalue_fname)
        self._read_topbase_level(level_fname)
        self._read_topbase_crosssections(photocx_fname)

    def _compute_level_index_key(self, slp, lv):
        return slp *1000 +lv

    def _creat_new_level_index(self,atom, ion, slp, lv):
        key = self._compute_level_index_key(slp, lv)
        if atom not in self.index_dict:
            self.index_dict[atom] = {}
            self.index_dict[atom][ion] = {}
            self.index_dict[atom][ion]['max'] = 0
            self.index_dict[atom][ion]['level_index'] = {}
        else:
            if ion not in self.index_dict[atom]:
                self.index_dict[atom][ion] = {}
                self.index_dict[atom][ion]['max'] = 0
                self.index_dict[atom][ion]['level_index'] = {}


        level_dic = self.index_dict[atom][ion]
        if key in level_dic['level_index']:
            return level_dic['level_index'][key]
        else:
            index = level_dic['max'] +1
            level_dic['max']= index
            level_dic['level_index'][key] = index
            self.index_dict[atom][ion]['level_index'][key] = index
            return index


#    @np.vectorize(excluded=['self'])
    def _get_level_index(self, atom, ion, slp, lv):
        key = self._compute_level_index_key(slp, lv)
        return self.index_dict[atom][ion]['level_index'][key]

    _get_level_index_vec = np.vectorize(_get_level_index, excluded=['self'])

    def read_level(self,fname):
       self.level_data = self._read_topbase_level(fname)

    def read_fvalues(self,fname):
        self.f_data = self._read_topbase_fvalues(fname)

    def read_photoionization_cross_section(self,fname):
        self.phot_crosssections_data = self._read_topbase_crosssections(fname)

    def __ion_atom_filter(self,atomic_number,ion_number, data):
        if isinstance(data, pd.DataFrame):
            amask = data.__array__()[:,0] < 0
            imask = amask.copy()
        elif isinstance(data, tuple):
            _one = self. __ion_atom_filter(atomic_number,ion_number, data[0])
            _two = self. __ion_atom_filter(atomic_number,ion_number, data[1])
            return (_one, _two)
        else:
            amask  = data[:]['i'] < 0
            imask  = data[:]['i'] < 0
        if isinstance(atomic_number,list):
            for n in atomic_number:
                _tmp = data[:]['atomic_number'] == n
                amask = np.logical_or(amask, _tmp)
        else:
            amask =np.logical_or(amask, data[:]['atomic_number'] == atomic_number)

        if isinstance(ion_number, list):
            for d in ion_number:
                _tmp = data[:]['ion_number'] = d
                imask = np.logical_or(imask, _tmp)
        else:
            imask = np.logical_or(imask, data[:]['ion_number'] == ion_number)

        return data[np.logical_and(amask,imask)]

    def get_levels(self, atomic_number, ion_number):
        return self.__ion_atom_filter(atomic_number, ion_number,self.level_data)

    def get_fvalues(self, atomic_number, ion_number):
        return self.__ion_atom_filter(atomic_number, ion_number, self.f_data)

    def get_photo_crosssection(self, atomic_number, ion_number):
        return self.__ion_atom_filter(atomic_number, ion_number, self.phot_crosssections_data)

    def _read_topbase_level(self,fname):
        e_datalist = []

        e_pattern_value = re.compile(r"\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+")
        pattern_unused = re.compile(r"^=|^[A-Z,a-z]")

        try:
            tmp_datat = pickle.load(open('tmp_level.txt', 'rb'))
            tmp_datat_index = pickle.load(open('tmp_level_index.txt', 'rb'))
            self.index_dict = tmp_datat_index
            print('Using tmp_level file')
        except IOError:
            tmp_datat = None

        if tmp_datat is None:
            with open(fname, 'r') as efile:
                for  data in efile:
                    if pattern_unused.match(data):
                        print("unused line")
                    elif e_pattern_value.match(data):
                        I = data[:7].strip()
                        NZ = data[8:10].strip()
                        NE = data[10:13]
                        iSLP = data[13:18].strip()
                        iLV = data[18:22].strip()
                        iCONF = data[23:39].strip()
                        Eryd = units.Unit('Ry').to('eV',float(data[40:51].strip()))
                        TEryd = units.Unit('Ry').to('eV', float(data[52:64].strip()))
                        gi = data[64:69].strip()
                        _tmp = [I,NZ, NE,iSLP,iLV,iCONF,Eryd,TEryd,gi]
                        e_datalist.append(_tmp)

                length = len(e_datalist)
                e_dataarray = np.recarray((length,),dtype=[('i',int),('level_number', int),('atomic_number', int),('ion_number', int),('iSLP', int),('iLV', int),
                                            ('iCONF', str),('E', np.float64),('energy', np.float64),('g', int)])

                for i,data in enumerate(e_datalist):
                    try:

                        e_dataarray[i]['i'] = int(data[0])
                        e_dataarray[i]['atomic_number'] = int(data[1])
                        e_dataarray[i]['ion_number'] = int(data[2])
                        e_dataarray[i]['iSLP'] = int(data[3])
                        e_dataarray[i]['iLV'] = int(data[4])
                        e_dataarray[i]['level_number'] = self._creat_new_level_index(int(data[1]),int(data[2]),int(data[3]),int(data[4]))#int('%d%d'%(e_dataarray[i]['iSLP'],e_dataarray[i]['iLV']))
                        e_dataarray[i]['iCONF'] = str(data[5])
                        e_dataarray[i]['E'] = np.float64(data[6]) #E
                        e_dataarray[i]['energy'] = np.float64(data[7]) #TE
                        e_dataarray[i]['g'] = int(np.float64(data[8]))
                    except ValueError:
                        print(data)
                        print(i)


            levels_data = pd.DataFrame(e_dataarray)
            levels_data.set_index('level_number', inplace=True)
            pickle.dump( levels_data, open('tmp_level.txt', 'wb'))
            pickle.dump( self.index_dict, open('tmp_level_index.txt', 'wb'))
        else:
            levels_data = tmp_datat

        return levels_data


    def _read_topbase_crosssections(self, fname):
        c_datalist = []
        c_datathlist = []

        c_pattern_unused = re.compile(r"^=|^[A-Z,a-z]")
        c_pattern_newline = re.compile(r"\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([+-]\d+\.\d+[Ee][+-]?\d+)\s+(\d+)")
        c_pattern_value = re.compile(r"\s*([+-]?\d+\.\d+[Ee][+-]?\d+)\s+([+-]?\d+\.\d+[Ee][+-]?\d+)")
        try:
            tmp_datat = pickle.load(open('tmp_photocx.txt', 'rb'))
            print('Using tmp_photocx file')
        except IOError:
            tmp_datat = None

        if tmp_datat is None:

            with open(fname,'r') as cfile:
                for data in cfile:
                    if c_pattern_unused.match(data):
                        print("unused line")
                    elif c_pattern_newline.match(data):
                        splitt_data = c_pattern_newline.match(data).groups()
                        I = splitt_data[0]
                        nz = splitt_data[1]
                        ne = splitt_data[2]
                        ISLP = splitt_data[3]
                        ILV = splitt_data[4]
                        Eion = units.Unit('Ry').to('eV',np.float64(splitt_data[5]))
                        number = splitt_data[6]
                        i = 0
                        print(splitt_data)
                    elif c_pattern_value.match(data):
                        i+=1
                        if i > number:
                            raise
                        splitt_data = c_pattern_value.match(data).groups()
                        Ecx = units.Unit('Ry').to('eV',np.float64(splitt_data[0]))
                        cx =  units.Unit('Mbarn').to('cm**2',float(splitt_data[1]))
                        _tmp = [I,nz,ne,ISLP,ILV,Eion,number,Ecx,cx]
                        c_datalist.append(_tmp)

                        if i ==1:
                            c_datathlist.append(_tmp)


            length = len(c_datalist)
            c_dataarray = np.recarray((length,),dtype=[('i', int),('atomic_number', int), ('ion_number',int), ('level_number', int), ('iSLP', int), ('iLV', int),
                                                       ('Eion', np.float64), ('number', int), ('Ecx',np.float64),('nu',np.float64) ,('cx', np.float64)])
            for i,data in enumerate(c_datalist):
                c_dataarray[i]['i'] = int(data[0])
                c_dataarray[i]['atomic_number'] = int(data[1])
                c_dataarray[i]['ion_number'] = int(data[2])
                c_dataarray[i]['iSLP'] = int(data[3])
                c_dataarray[i]['iLV'] = int(data[4])
                c_dataarray[i]['level_number'] = self._get_level_index(int(data[1]),int(data[2]),int(data[3]),int(data[4]))#int('%d%d'%(c_dataarray[i]['iSLP'],c_dataarray[i]['iLV']))
                c_dataarray[i]['Eion'] = np.float64(data[5])
                c_dataarray[i]['number'] = int(data[6])
                c_dataarray[i]['Ecx'] = units.Unit('Ry').to('eV',np.float64(data[7]))
                c_dataarray[i]['nu'] =  units.Unit('angstrom').to('Hz',(constants.h.cgs.value / c_dataarray[i]['Ecx']), units.spectral())
                c_dataarray[i]['cx'] = np.float64(data[8])

            length = len(c_datathlist)
            c_datatharray = np.recarray((length,),dtype=[('i', int),('atomic_number', int), ('ion_number',int),('level_number', int)
                                        , ('iSLP', int), ('iLV', int), ('Eion', np.float64), ('number', int), ('Ecx',np.float64),('nu', np.float64), ('cx', np.float64)])
            for i,data in enumerate(c_datathlist):
                c_datatharray[i]['i'] = int(data[0])
                c_datatharray[i]['atomic_number'] = int(data[1])
                c_datatharray[i]['ion_number'] = int(data[2])
                c_datatharray[i]['iSLP'] = int(data[3])
                c_datatharray[i]['iLV'] = int(data[4])
                c_dataarray[i]['level_number'] = self._get_level_index(int(data[1]),int(data[2]),int(data[3]),int(data[4]))#int('%d%d'%(c_datatharray[i]['iSLP'],c_datatharray[i]['iLV']))
                c_datatharray[i]['Eion'] = np.float64(data[5])
                c_datatharray[i]['number'] = int(data[6])
                c_datatharray[i]['Ecx'] = np.float64(data[7])
                c_datatharray[i]['Ecx'] = units.Unit('Ry').to('eV',np.float64(data[7]))
                c_datatharray[i]['nu'] = units.Unit('angstrom').to('Hz',(constants.h.cgs.value / c_datatharray[i]['Ecx']), units.spectral())
                c_datatharray[i]['cx'] = np.float64(data[8])


            pickle.dump( [c_dataarray, c_datatharray], open('tmp_photocx.txt', 'wb'))
        else:
            c_dataarray  = tmp_datat[0]
            c_datatharray = tmp_datat[1]

        return c_dataarray, c_datatharray

    def _read_topbase_fvalues(self, fname):

        try:
            tmp_datat = pickle.load(open('tmp_fvalue.txt', 'rb'))
            print('Using tmp_fvalue file')
        except IOError:
            tmp_datat = None

        if tmp_datat is None:
            f_datalist = []
            pattern_unused = re.compile(r"^=|^[A-Z,a-z]")
            f_pattern_value = re.compile(r"\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+")

            with open(fname,'r') as ffile:
                for data in ffile:
                    if pattern_unused.match(data):
                        print("unused line")
                    elif f_pattern_value.match(data):
                        I = data[:7].strip()
                        NZ = data[8:10].strip()
                        NE = data [11:13].strip()
                        iSLP = data[14:18].strip()
                        jSLP = data[19:23].strip()
                        iLV = data[24:27].strip()
                        jLV = data[28:31].strip()
                        gF = data[32:41].strip()
                        gA = data[42:51].strip()
                        WL = data[52:61].strip()
                        gi = data[62:66].strip()
                        gj = data[67:].strip()

                        _tmp = [I, NZ, NE, iSLP, jSLP, iLV, jLV, gF, gA, WL, gi, gj]
                        if np.abs(float(gF)) <= 0:
                            print('----')
                            print(_tmp)
                            print('f value ignored (abs(gf) <= 0 )')
                            print('----')
                        else:
                            f_datalist.append(_tmp)

                length = len(f_datalist)
                f_dataarray = np.recarray((length,),dtype=[('i', int),('atomic_number', int),('ion_number', int),('level_number_lower', int),
                                                           ('level_number_upper', int),('iLV', int),('jLV', int),('gF', np.float64),
                                          ('gA', np.float64),('wavelength', np.float64),('gi', int),('gj', int),('gu', int),('gl', int) ,('nu', np.float64), ('F_lu',np.float64), ('F_ul',np.float64),
                    ('loggf',np.float64), ('jSLP', int), ('iSLP', int),('is_ul', bool), ('test', float) ])




                @np.vectorize
                def compute_level_index(slp, lv):
                    return int('%d%d'%(slp, lv))

                for i,data in enumerate(f_datalist):
                     try:
                        # nu = units.Unit('angstrom').to('Hz', np.float64(data[9]), units.spectral())
                         #A_coeff = (8 * np.pi**2 * constants.e.gauss.value**2 * nu**2)/ (constants.m_e.cgs.value * constants.c.cgs.value**3)
                         #gF_lu = - np.float64(data[7])
                         #F_lu = gF_lu /  int(np.float64(data[11])) #f_lu/ g_l
                         #F_ul = np.float64(data[7]) / int(np.float64(data[10])) #f_ul/ g_u
                         #loggf = np.log10(F_ul*int(np.float64(data[10])))
                         f_dataarray[i]['i'] = int(data[0])
                         f_dataarray[i]['atomic_number'] = int(data[1])
                         f_dataarray[i]['ion_number'] = int(data[2])
                         f_dataarray[i]['iSLP'] = int(data[3])
                         f_dataarray[i]['jSLP'] = int(data[4])
                         f_dataarray[i]['iLV'] = int(data[5])
                         f_dataarray[i]['jLV'] = int(data[6])
                         f_dataarray[i]['gF'] = np.abs(np.float64(data[7]))
                         f_dataarray[i]['gA'] = np.float64(data[8])
                         f_dataarray[i]['wavelength'] = np.float64(data[9])
                         f_dataarray[i]['gi'] = int(np.float64(data[10]))
                         f_dataarray[i]['gj'] = int(np.float64(data[11]))
                         f_dataarray[i]['is_ul']  = np.float64(data[7]) >0
                         if f_dataarray[i]['is_ul']:
                             f_dataarray[i]['gu'] = f_dataarray[i]['gi']
                             f_dataarray[i]['gl'] = f_dataarray[i]['gj']
                         else:
                             f_dataarray[i]['gu'] = f_dataarray[i]['gj']
                             f_dataarray[i]['gl'] = f_dataarray[i]['gi']
                         #f_dataarray[i]['level_number_lower'] = int('%d%d'%(int(data[3]),int(data[5])))
                         #f_dataarray[i]['level_number_upper'] = int('%d%d'%(int(data[4]),int(data[6])))
                     except ValueError:
                         print(data)
                         print(i)

                f_dataarray[:]['nu'] = units.Unit('angstrom').to('Hz',f_dataarray[:]['wavelength'], units.spectral())
                f_dataarray[:]['F_lu'] = f_dataarray[:]['gF'] / map(float, f_dataarray[:]['gl'])
                f_dataarray[:]['F_ul'] = f_dataarray[:]['gF']  / map(float, f_dataarray[:]['gu'])
                f_dataarray[:]['loggf'] = np.log10(f_dataarray[:]['F_ul'] * map(float, f_dataarray[:]['gu']))
                f_dataarray[:]['test'] = (f_dataarray[:]['F_ul'] * map(float, f_dataarray[:]['gu'])) - f_dataarray[:]['gF']
                f_dataarray[:]['level_number_lower'] = self._get_level_index_vec(self ,f_dataarray[:]['atomic_number'],f_dataarray[:]['ion_number'],f_dataarray[:]['iSLP'], f_dataarray[:]['iLV'])#compute_level_index(f_dataarray[:]['iSLP'], f_dataarray[:]['iLV'])
                f_dataarray[:]['level_number_upper'] = self._get_level_index_vec(self ,f_dataarray[:]['atomic_number'],f_dataarray[:]['ion_number'],f_dataarray[:]['jSLP'], f_dataarray[:]['jLV'])#compute_level_index(f_dataarray[:]['jSLP'], f_dataarray[:]['iLV'])

                pickle.dump( f_dataarray, open('tmp_fvalue.txt', 'wb'))
        else:
            f_dataarray = tmp_datat

        return f_dataarray





def insert_one_species_to_db(symbol, ion_number, conn, topbasedatastore, temperatures):


    atomic_number = int(symbol2z[symbol])

    print('deleting atom=%d ion=%d'%(atomic_number, ion_number))
    curs = conn.cursor()
    curs.execute('delete from levels where atom=? and ion=?', (atomic_number, ion_number - 1 ))
    curs.execute('delete from lines where atom=? and ion=?', (atomic_number, ion_number  - 1))

    #collision_data_cols = curs.execute('pragma table_info(collision_data)').fetchall()

    #if temperatures is None:
    #    print "Trying to infer temperatures from column names"

    #if collision_data_cols == []:
    #    raise IOError('The given database doesn\'t contain a collision_data table - please create it')

    #else:
    #    temperatures_data = [int(item.strip('t')) for item in zip(*collision_data_cols)[1] if item.startswith('t')]
    #    print "Inferred temperatures are %s" % (temperatures_data,)


    datastore = topbasedatastore
    levels_data = datastore.get_levels(atomic_number,ion_number)
    lines_data = datastore.get_fvalues(atomic_number,ion_number)
    cx_data, cx_datath = datastore.get_photo_crosssection(atomic_number,ion_number)
    #collision_data = datastore. # set to zero

    for line in lines_data:
        print(line)
        curs.execute('insert into lines(wl, atom, ion, level_id_upper, level_id_lower, f_lu, f_ul, loggf, g_upper, g_lower) '
                     'values(?, ?, ?, ?, ?, ?, ?, ?,?,?)',
                     (line['wavelength'], line['atomic_number'], line['ion_number'] -1 ,
                      line['level_number_upper'], line['level_number_lower'],
                      line['F_lu'], line['F_ul'], line['loggf'],line['gu'],line['gl'] ))


    for key, level in levels_data.iterrows():
        count_down = curs.execute('select count(id) from lines where atom=? and ion=? and level_id_upper=?',
                     (atomic_number, ion_number - 1, int(key-1))).fetchone()[0]

        curs.execute('insert into levels(atom, ion, energy, g, level_id, metastable) values(?, ?, ?, ?, ?, ?)',
                     (atomic_number, ion_number-1, level['energy'], level['g'], int(key-1), count_down == 0))


    photocx_data_tabel_stmt = r"create table cx_data(id integer primary key," \
                              r"atom integer," \
                              r"ion integer," \
                              r"level_number integer," \
                              r"cross_section float)"

    photocxsp_data_tabel_stmt = r"create table cx_sp_data(id integer primary key," \
                              r"atom integer," \
                              r"ion integer," \
                              r"level_number integer," \
                              r"nu float," \
                              r"cross_section float)"

    curs.execute(photocx_data_tabel_stmt)
    curs.execute(photocxsp_data_tabel_stmt)


    for data in cx_data:
        curs.execute('insert into cx_data(atom, ion, level_number, cross_section)values(?,?,?,?)',
                     (data['atomic_number'], data['ion_number'] - 1, data['level_number'], data['cx']))

    for data in cx_datath:
        curs.execute('insert into cx_sp_data(atom, ion, level_number, nu ,cross_section) values(?,?,?,?,?)',
                     (data['atomic_number'], data['ion_number'] - 1, data['level_number'],data['nu'] , data['cx']))

    conn.commit()



def insert_to_db(species,conn, fvalue_fname, level_fname, photocx_fname, temperatures=None):
    topbasedatastore = TopBaseData()
    topbasedatastore.read_level(level_fname)
    topbasedatastore.read_fvalues(fvalue_fname)
    topbasedatastore.read_photoionization_cross_section(photocx_fname)

    for species in species:
        print("Inserting %s %s to database" % (str(species[0]), str(species[1])))
        insert_one_species_to_db(species[0],species[1],conn ,topbasedatastore, temperatures)
