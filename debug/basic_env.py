import h5py
import os
import sys
import pandas as pd


from astropy import units, constants, table

import cPickle as pickle
from tardisatomic import sql_stmts
import numpy as np
import sqlite3
from collections import OrderedDict

import hashlib
import uuid
import itertools

from tardis import util

f = h5py.File("atom.hd5", "r")
levels_data = f['levels_data'].value
lines_data = f['lines_data'].value
ionization_data = f['ionization_data'].value

level_index = ['atomic_number','ion_number', 'level_number']
lines_index =['atomic_number','ion_number','line_id']
ionization_index = ['atomic_number','ion_number']


levels_pd = pd.DataFrame(levels_data)
levels_pd.set_index(level_index, inplace=True)
lines_pd = pd.DataFrame(lines_data)
lines_pd.set_index(lines_index,inplace=True)
ionization_pd = pd.DataFrame(ionization_data)
ionization_pd.set_index(ionization_index, inplace=True)




