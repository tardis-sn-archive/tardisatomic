# getting the linelist, atoms data and ionization data into one file.
import sqlite3
import os
import re
import time
import numpy as np
def read_raw_wiki_ionize(fname):
    raw_data = file(fname).readlines()
    elements = [int(re.search('^\d+',item).group()) for item in raw_data[::4]]
    ion_data = []
    for elem, raw_ion_data in zip(elements, raw_data[3::4]):
        for ion, chi in enumerate(raw_ion_data.split()[1:]):
            ion_data += [[elem, ion + 1, float(chi)]]
    return ion_data

os.system('rm kurucz_atoms.db3')
try:
    conn.close()
except:
    pass
    
conn = sqlite3.connect('kurucz_atoms.db3')


##### IONIZATION DATA READ
print "Reading ionization data"
ionization_stmt = """create table ionization(id integer primary key,
                        atom integer,
                        ion integer, 
                        ionize_ev float)"""

conn.execute(ionization_stmt)

for atom, ion, energy in read_raw_wiki_ionize('wikipedia_ionization_energies.dat'):
    conn.execute('insert into ionization(atom, ion, ionize_ev) values(?, ?, ?)', (atom, ion, energy))


##### ATOMIC DATA READ
print "Reading general atomic data"
atomic_wt_stmt = """create table atoms(atom integer primary key,
                symbol string,
                name string,
                weight float)"""


conn.execute(atomic_wt_stmt)
for z, symbol, name, weight in np.genfromtxt('atomic_wts_formatted.dat', dtype=(np.int64, 'S2', 'S20', np.float64)):
    conn.execute('insert into atoms(atom, symbol, name, weight) values(?, ?, ?, ?)', (int(z), symbol, name, weight))
    


##### LINE LIST  READ

print "reading line list"
hc = 4.135667516e-15 * 2.99792458e10
llist_create_stmt = """create table main.lines(wl FLOAT,
                                loggf float,
                                atom integer,
                                ion integer,
                                e_upper float,
                                g_upper integer,
                                label_upper text,
                                level_id_upper integer default -1,
                                e_lower float,
                                g_lower integer,
                                label_lower text,
                                level_id_lower integer default -1)"""
conn.execute(llist_create_stmt)                                
conn.execute("attach 'gfall.db3' as kurucz_lines")
llist_select_stmt = """select 10*wl, loggf, elem as atom, ion, 
                    e_upper * ? , cast(2*j_upper + 1 as integer) as g_upper, label_upper, 
                    e_lower *  ?, cast(2*j_lower + 1 as integer) as g_lower, label_lower
                from kurucz_lines.gfall
                where loggf > -3
                and elem in (6, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30)"""
             
llist_insert_stmt = """insert into main.lines(wl, loggf, atom, ion, e_upper, g_upper, label_upper, e_lower, g_lower, label_lower) """
conn.execute(llist_insert_stmt + llist_select_stmt, (hc, hc))


print "%d rows remaining" % (conn.execute('select count(atom) from lines').fetchone()[0])

print "Creating level and level indices"

level_select_stmt= """select * from (select distinct atom, ion, e_upper as energy, g_upper as g, label_upper as label from lines
                        union
                      select distinct atom, ion, e_lower as energy, g_lower as g, label_lower as label from lines) order by atom, ion, energy"""
                      
level_create_stmt = "create table levels(atom integer, ion integer, energy float, g integer, label text, level_id integer)"
conn.execute('drop table if exists levels')
conn.execute(level_create_stmt)

curs = conn.execute(level_select_stmt)

old_atom = None
old_ion = None
levels_data = []
for atom, ion, energy, g , label in curs:
    if atom == old_atom and ion == old_ion:    
        i +=1
    else:
        old_atom = atom
        old_ion = ion
        i=0
    conn.execute('insert into levels(atom, ion, energy, g, label, level_id) values(?, ?, ?, ?, ?, ?)', (atom, ion, energy, g , label, i)) 
conn.execute('create index level_unique_idx on levels(atom, ion, energy, g, label)')
conn.commit()

#LINKING LINES TO LEVELS
print "linking levels -- this may take several minutes"
starttime = time.time()
update_link_stmt = """
        update
            lines
        set
            level_id_upper = 
            (select level_id from levels
                where lines.atom=levels.atom and
                    lines.ion=levels.ion and
                    lines.e_upper=levels.energy and
                    lines.g_upper=levels.g and
                    lines.label_upper=levels.label),
            level_id_lower = 
            (select level_id from levels
                where lines.atom=levels.atom and
                    lines.ion=levels.ion and
                    lines.e_lower=levels.energy and
                    lines.g_lower=levels.g and
                    lines.label_lower=levels.label)
        """
conn.execute(update_link_stmt)

conn.commit()
conn.close()
