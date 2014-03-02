*********************
Running TARDIS Atomic
*********************

TARDIS Atomic works with two different file formats. The main operations are done on an sqlite database. In the last step
this sqlite database is then converted to the TARDIS hdf5 file.

Generate initial Kurucz Database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Download the gfall.dat from `<http://kurucz.harvard.edu/LINELISTS/GFALL/>`_. The next step is to convert `gfall.dat` to
its sqlite form `gfall.db3`::

    python tardisatomic/scripts/read_gfall2db gfall.dat gfall.db3

The next step is to convert the `gfall.db3` to a real atomic database (with levels and lines and links between them). This is
done by `make_kurucz_db`::

    python tardisatomic/scripts/make_kurucz_db kurucz.db3 gfall.db3
