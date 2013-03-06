import os
import sqlite3
import math
import numpy as np
import macro_atom_transition
from tardisatomic import import_ionDB

gfall_db = os.path.join(os.path.dirname(__file__), 'data', 'gfall.db3')
zeta_datafile = os.path.join(os.path.dirname(__file__), 'data', 'knox_long_recombination_zeta.dat')

try:
    import sqlparse
    sqlparse_available = True
except ImportError:
    sqlparse_available = False


#constants

hc = 4.135667516e-15 * 2.99792458e10 # 2 ev


linelist_select_stmt = """
SELECT
    10*wl,
    loggf,
    elem AS atom,
    ion,
    e_upper * %(hc).20f,
    cast(2*j_upper + 1 AS integer) AS g_upper,
    label_upper,
    e_lower * %(hc).20f,
    cast(2*j_lower + 1 AS integer) AS g_lower,
    label_lower
FROM
    kurucz_lines.gfall
"""

linelist_insert_stmt = """
insert into
    main.lines(wl,
        loggf,
        atom,
        ion,
        e_upper,
        g_upper,
        label_upper,
        e_lower,
        g_lower,
        label_lower)"""
        
        
linelist_create_stmt = """
CREATE TABLE
    main.lines(
    id integer primary key,
    wl float,
    loggf float,
    atom integer,
    ion integer,
    e_upper float,
    g_upper integer,
    label_upper text,
    level_id_upper integer default -1,
    global_level_id_upper integer default -1,
    f_ul float,
    e_lower float,
    g_lower integer,
    label_lower text,
    level_id_lower integer default -1,
    global_level_id_lower integer default -1,
    f_lu float)
    """

update_oscillator_stmt = """
UPDATE
    lines
SET
    f_ul = pow(10, loggf) / g_upper,
    f_lu = pow(10, loggf) / g_lower
"""
    
    
def new_linelist_from_gfall(new_dbname, gfall_fname=None, select_atom=None):
    print "Reading lines from Kurucz gfall"
    conn = sqlite3.connect(new_dbname)
    conn.create_function('pow', 2, math.pow)
    if gfall_fname is None:
        gfall_fname = gfall_db
    #attaching gfall database
    conn.execute("attach '%s' as kurucz_lines" % gfall_fname)
    curs = conn.cursor()
    curs.execute('drop table if exists lines')
    curs.execute(linelist_create_stmt)
    if select_atom is None:
        elem_select_stmt = ""
    else:
        elem_select_stmt = " and elem in (%s)" % (','.join(map(str, select_atom)),)
    insert_fromgfall_stmt = linelist_insert_stmt + linelist_select_stmt % {'hc':hc, 'where_stmt':elem_select_stmt}
    
    if sqlparse_available:
        print sqlparse.format(insert_fromgfall_stmt, reindent=True)
    else:
        print insert_fromgfall_stmt
    
    curs.execute(insert_fromgfall_stmt)
    
    conn.commit()
    print "%d lines in database" % (conn.execute('select count(atom) from lines').fetchone()[0])
    print "updating oscillator strengths"
    
    conn.execute(update_oscillator_stmt)
    conn.commit()
    return conn


level_create_stmt = """
CREATE TABLE
    levels(
        id integer primary key,
        atom integer,
        ion integer,
        energy float,
        g integer,
        label text,
        level_id integer,
        metastable bool default NULL)"""


level_select_stmt= """
SELECT
    *
FROM
    (SELECT DISTINCT
        atom,
        ion,
        e_upper as energy,
        g_upper as g,
        label_upper as label
    FROM
        lines
    UNION
    
    SELECT DISTINCT
        atom,
        ion,
        e_lower as energy,
        g_lower as g,
        label_lower as label 
    FROM
        lines)
    ORDER BY atom, ion, energy
    """

def add_fully_ionized_levels(conn):
    #Clean first
    print "Adding fully ionized levels for H and He"
    clean_fully_ionized_stmt = "DELETE FROM levels WHERE atom == ion "
    conn.execute(clean_fully_ionized_stmt)

    add_full_ionized_levels_stmt ="""
INSERT INTO
    levels (atom,ion,energy,g,metastable,level_id)
    SELECT DISTINCT
        atom,atom,0,1,1,0
    FROM
        levels
    WHERE NOT EXISTS
        (SELECT 1 FROM levels WHERE atom == ion)
    AND
        atom in (1, 2)
    """
    conn.execute(add_full_ionized_levels_stmt)
    conn.commit()


def create_levels(conn):
    print "Creating level and level indices"
    curs = conn.cursor()
    curs.execute('drop table if exists levels')
    curs.execute(level_create_stmt)
   #TODO: The level id is not nice, but necessary for the make_kurucz_hdf5.
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
    conn.execute('create index level_global_idx on levels(id)')
    print('Creating fully ionized levels')
    return conn

def link_levels(conn):
    print "linking levels -- this may take several minutes"
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
    

    update_global_link_stmt = """
        update
            lines
        set
            global_level_id_upper = 
            (select id from levels
                where lines.atom=levels.atom and
                    lines.ion=levels.ion and
                    lines.level_id_upper=levels.level_id),
                    
            global_level_id_lower = 
            (select id from levels
                where lines.atom=levels.atom and
                    lines.ion=levels.ion and
                    lines.level_id_lower=levels.level_id)
        """



    conn.execute(update_link_stmt)
    conn.execute('create index line_level_upper_idx on lines(atom, ion, level_id_upper)')
    conn.execute('create index line_level_lower_idx on lines(atom, ion, level_id_lower)')
    conn.execute('create index level_idx on levels(atom, ion, level_id)')
    conn.execute('create index global_level_idx on levels(id)')
    print "Linking lines and levels with global id"
    conn.execute(update_global_link_stmt)
    conn.execute('')
    conn.commit()
    return conn




create_macro_atom_stmt="""
    create table macro_atom(
        id integer_primary_key references levels(id),
        count_down default NULL,
        count_up default NULL,
        reference_down int_ndarray default NULL,
        reference_up int_ndarray default NULL,
        line_id_down int_ndarray default NULL,
        line_id_up int_ndarray default NULL,
        p_internal_down float_ndarray,
        p_emission_down float_ndarray,
        p_internal_up float_ndarray
    )"""


down_transition_stmt = """
SELECT
    global_level_id_upper as id,
    count(global_level_id_lower) as number_down,
    group_concat_intarray(global_level_id_lower) as reference_down,
    group_concat_intarray(id) as line_id_down,
    calculate_p_internal_down(wl, g_lower, g_upper, f_lu, e_lower) as p_int_down,
    calculate_p_emission_down(wl, g_lower, g_upper, f_lu, e_upper - e_lower) as p_em_down
FROM
    lines 
GROUP BY
    global_level_id_upper
"""

up_transition_stmt = """
SELECT
    global_level_id_lower as id,
    count(global_level_id_upper) as number_up,
    group_concat_intarray(global_level_id_upper) as reference_up,
    group_concat_intarray(id) as line_id_up,
    calculate_p_internal_up(wl, f_lu, e_lower) as p_int_up
FROM
    lines
GROUP BY
    global_level_id_lower

"""

update_up_stmt="""
update macro_atom
SET count_up = 
    (select 
        number_up 
     from
        up_transitions
     where
        up_transitions.id=macro_atom.id),
        
    reference_up = 
    (select 
        reference_up 
     from
        up_transitions
     where
        up_transitions.id=macro_atom.id),

    line_id_up = 
    (select 
        line_id_up 
     from
        up_transitions
     where
        up_transitions.id=macro_atom.id),

    p_internal_up = 
    (select 
        p_int_up 
     from
        up_transitions
     where
        up_transitions.id=macro_atom.id)
    """
    
update_down_stmt="""
update macro_atom
SET count_down = 
    (select 
        number_down 
     from
        down_transitions
     where
        down_transitions.id=macro_atom.id),
        
    reference_down = 
    (select 
        reference_down 
     from
        down_transitions
     where
        down_transitions.id=macro_atom.id),
    line_id_down = 
    (select 
        line_id_down
     from
        down_transitions
     where
        down_transitions.id=macro_atom.id),

    p_internal_down = 
    (select 
        p_int_down
     from
        down_transitions
     where
        down_transitions.id=macro_atom.id),

    p_emission_down = 
    (select 
        p_em_down 
     from
        down_transitions
     where
        down_transitions.id=macro_atom.id)
    """




def create_temporary_transition_table(conn):
    curs = conn.cursor()
    conn.create_aggregate('calculate_p_internal_down', 5, macro_atom_transition.calculate_p_internal_down)
    conn.create_aggregate('calculate_p_emission_down', 5, macro_atom_transition.calculate_p_emission_down)
    conn.create_aggregate('calculate_p_internal_up', 3, macro_atom_transition.calculate_p_internal_up)
    conn.create_aggregate('group_concat_intarray', 1, macro_atom_transition.group_concat_intarray)
    print "Creating temporary tables with down and up transition probabilities"
    curs.execute('create temporary table down_transitions as %s' % down_transition_stmt)
    curs.execute('create temporary table up_transitions as %s' % up_transition_stmt)
    curs.execute('create index down_idx on down_transitions(id)')
    curs.execute('create index up_idx on up_transitions(id)')
    conn.commit()
    return conn

    

update_macro_atom_clean = """
UPDATE
    macro_atom
SET
    count_%(trans_type)s = 0,
    reference_%(trans_type)s = ?,
    line_id_%(trans_type)s = ?,
    p_internal_%(trans_type)s = ?
    %(p_emission_stmt)s
WHERE
    reference_%(trans_type)s isnull"""

def create_macro_atom(conn):
    empty_array = -1
    curs = conn.cursor()
    print "Creating and populating the macro atom table"
    curs.execute('drop table if exists macro_atom')    
    curs.execute(create_macro_atom_stmt)
    curs.execute('insert into macro_atom(id) select id from (select id from down_transitions union select id from up_transitions)')
    curs.execute(update_down_stmt)
    curs.execute(update_up_stmt)
    print "Cleaning up macro atom table and setting counts to 0 where transitions don't exist"
    print update_macro_atom_clean % dict(trans_type='down', p_emission_stmt=', p_emission_down = ?')
    curs.execute(update_macro_atom_clean % dict(trans_type='down', p_emission_stmt=', p_emission_down = ?'),
                    (empty_array, empty_array, empty_array, empty_array))
    
    curs.execute(update_macro_atom_clean % dict(trans_type='up', p_emission_stmt=''),
                    (empty_array, empty_array, empty_array))
    conn.execute('create index macro_atom_global_idx on macro_atom(id)')
    conn.commit()
    return conn
    
def flag_metastable(conn):
    flag_metastable_stmt = """
    UPDATE
        levels
    SET
        metastable=(SELECT CASE
                    WHEN
                        count_down = 0
                    THEN 1
                    WHEN
                        count_down > 0
                    THEN 0
                    END
              FROM
                macro_atom
              WHERE
                macro_atom.id=levels.id)
    """
    print "Flagging metastable levels"
    conn.execute(flag_metastable_stmt)
    conn.commit()
    return conn


create_zeta_table_stmt = """
CREATE TABLE zeta(
    id integer primary key,
    atom integer,
    ion integer,
    zeta float_ndarray)"""

def read_zeta(conn):
    print "reading zeta values and inserting into db from %s" % zeta_datafile
    zeta_data = np.loadtxt(zeta_datafile, usecols=xrange(1,23), dtype=np.float64)
    conn.execute(create_zeta_table_stmt)
    for line in zeta_data:
        atom = int(line[0])
        ion = int(line[1])
        z_data = line[2:]
        conn.execute('insert into zeta(atom, ion, zeta) values(?, ?, ?)',
                     (atom, ion, sqlite3.Binary(z_data.tostring())))
    conn.commit()
    return conn

def ion_xs(conn):
    print("Creating the ionization cross section table.")
    ion_xs_create_table = """
    CREATE TABLE ion_cx(id INTEGER PRIMARY KEY, level_id INTEGER, atom INTEGER, ion INTEGER, cx_threshold FLOAT, FOREIGN KEY(level_id, atom, ion) REFERENCES levels(level_id, atom, ion))
    """
    ion_xs_creat_supporter = """
    CREATE TABLE ion_cx_supporter(id INTEGER PRIMARY KEY, ion_cx_id INTEGER , nu FLOAT , xs FLOAT, FOREIGN KEY(ion_cx_id) REFERENCES ion_cx(id))
    """


    curs = conn.cursor()
    #curs.execute('PRAGMA foreign_keys = ON')
    curs.execute(ion_xs_create_table)
    curs.execute(ion_xs_creat_supporter)


    atomdata = np.array(curs.execute('SELECT atom, ion, level_id FROM levels').fetchall())

    cx = import_ionDB.analytic_cross_section(atomdata[:,0],atomdata[:,0] - atomdata[:,1],0.8)
    final = np.concatenate((atomdata,cx[:,None]),axis=1)
    for i in final:
        curs.execute('INSERT OR IGNORE INTO ion_cx (level_id, atom, ion, cx_threshold) VALUES (?,?,?,?)',(int(i[2]),int(i[0]),int(i[1]),i[3]))
    conn.commit()
    return conn






    