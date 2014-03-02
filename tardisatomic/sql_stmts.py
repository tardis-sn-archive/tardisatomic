__author__ = 'wkerzend'

linelist_select_stmt = """
SELECT
    convert_air2vacuum(10*wavelength),
    loggf,
    atomic_number,
    ion_number,
    e_upper * %(hc).20f,
    cast(2*j_upper + 1 AS integer) AS g_upper,
    label_upper,
    e_lower * %(hc).20f,
    cast(2*j_lower + 1 AS integer) AS g_lower,
    label_lower,
    "kurucz"
FROM
    kurucz_lines.gfall
"""


linelist_insert_stmt = """
insert into
    main.lines(wavelength,
        loggf,
        atom,
        ion,
        e_upper,
        g_upper,
        label_upper,
        e_lower,
        g_lower,
        label_lower,
        source)"""


linelist_create_stmt = """
CREATE TABLE
    main.lines(
    id integer primary key,
    wavelength float,
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
    f_lu float,
    source text)
    """

update_oscillator_stmt = """
UPDATE
    lines
SET
    f_ul = pow(10, loggf) / g_upper,
    f_lu = pow(10, loggf) / g_lower
"""