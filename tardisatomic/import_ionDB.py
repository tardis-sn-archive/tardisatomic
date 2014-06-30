__author__ = 'michi'

import tardisatomic
from astropy import constants
import math

import argparse
import numpy as np
import h5py
import os
import sys
import sqlite3
import re


class Cache(object):

    def __init__(self):
        self.cache = {}
        self.func = None

    def cachedFunc(self, *args):
        if args not in self.cache:
            self.cache[args] = self.func(*args)
        return self.cache[args]

    def __call__(self, func):
        self.func = func
        return self.cachedFunc



@np.vectorize
@Cache()
def analytic_cross_section(n,z,g):
    """
   Computes the photoionization cross sections based on  Karzas and Latter 1961.

    Parameters
    ----------

    n : `int`
       Main quantum number

    z : `int`
       the charge in


    g : 'float'
       Kramers-Gaunt factor dimension-less



    """


    fineStructureConstant = (constants.e.gauss.value)**2/(constants.hbar.cgs.value * constants.c.cgs.value)
    bohrRadius = constants.hbar.cgs.value / (constants.m_e.cgs.value * constants.c.cgs.value * fineStructureConstant)
    sigam_edge = ((64 * math.pi * n * g)/(3 * math.sqrt(3)*np.power(z,2))) * fineStructureConstant * np.power(bohrRadius,2)

    return sigam_edge
