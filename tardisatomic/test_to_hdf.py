from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
#execfile('tardisatomic/alchemy/ingest/kurucz.py')
#from tardisatomic.alchemy.ingest.kurucz import *

from tardisatomic.alchemy.to_hd5.createhd5 import CreateHd5
from tardisatomic.base import AtomicDatabase
import os
import shutil

atomic_db = AtomicDatabase('sqlite:///example3.db')

hdf = CreateHd5(atomic_db, 'test.hdf')
hdf._load_transitions()
