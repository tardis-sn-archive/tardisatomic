#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse

import numpy as np

from tardisatomic.alchemy.to_hd5.createhd5 import CreateHDF

parser = argparse.ArgumentParser()

parser.add_argument('hdf5', help='HDF5 file that is created')
parser.add_argument('atom_db', help='atom_db')

parser.add_argument('--loggf_threshold', default=-3.,
                    help='log(gf) threshold for the kurucz dataset')
parser.add_argument('--max_atomic_number', default=30, type=int, help=(
    'Setting the maximum atomic number to be stored in the HDF5'
    'file. For now 30 is the limit as I don\'t have reliable ionization data '
    'for '
    'anything above.'))

parser.add_argument('--exclude_atoms', default=False, action="store_true",
                    help='Exclude basic atomic data')

parser.add_argument('--max_ionization_level', default=np.inf, type=int)

parser.add_argument('--exclude_ions', default=False, action="store_true",
                    help='Exclude ion data')
parser.add_argument('--exclude_levels', default=False, action="store_true",
                    help='Exclude level data')
parser.add_argument('--exclude_lines', default=False, action="store_true",
                    help='Exclude line data')
parser.add_argument('--anonymous', default=False, action="store_true",
                    help='Excludes all meta information from the hdf file ')
atom_group = parser.add_mutually_exclusive_group()

args = parser.parse_args()


def included_data_list(args):
    _tmp = []
    if not args.exclude_atoms is True:
        _tmp.append('atoms')
    if not args.exclude_ions is True:
        _tmp.append('ions')
    if not args.exclude_levels is True:
        _tmp.append('levels')
    if not args.exclude_ions is True:
        _tmp.append('lines')

    return _tmp


def open_hdf_file(args):
    atom_db = args.atom_db
    hdf_file = args.hdf_file
    included_data = included_data_list(args)
    exclude_species = args.exclude_atoms
    max_ionization_level = args.max_ionization_level

    return CreateHDF(atom_db, hdf_file, included_data,
                     exclude_species=exclude_species,
                     max_ionization_level=max_ionization_level)


hdf = open_hdf_file(args)
hdf.close_hdf(anonymous=args.anonymous)
