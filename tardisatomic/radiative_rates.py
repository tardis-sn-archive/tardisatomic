import mpmath as mp
import numpy as np
from astropy import constants
from scipy import special

def gamma_incomplete_0(z):
    return -special.expi(-z) + 0.5*(- np.log(1/z)+np.log(-z)) - np.log(z)


def calculate_spontaneous_recomb(nu_edge, T):
    return 8 * np.pi * nu_edge ** 3 / constants.c.cgs.value ** 2 * float(
        mp.gammainc(0, (nu_edge * constants.h.cgs.value) / (constants.k_B.cgs.value * T)))


def calculate_spontaneous_recomb_modified(nu_edge, T):
    return 8 * np.pi * nu_edge ** 2 / constants.c.cgs.value ** 2 / constants.h.cgs.value * \
    np.exp(-nu_edge * constants.h.cgs.value / constants.k_B.cgs.value / T)


def calculate_photoionization_coef(nu_edge, T):
    """
    calculate_photoionization_coef gives two values for the computation of the photoionization coefficient.
     The photoionization coefficient
    """
    f1 = 4 * np.pi / constants.h.cgs.value * nu_edge / 3
    helper1 = -nu_edge * constants.h.cgs.value / constants.k_B.cgs.value / T
    f2 = np.exp(helper1) * constants.k_B.cgs.value * T * (
    nu_edge ** 2 * constants.h.cgs.value ** 2 - nu_edge * constants.h.cgs.value * constants.k_B.cgs.value * T + 2 *
    constants.k_B.cgs.value ** 2 * T ** 2 ) + nu_edge ** 3 * constants.h.cgs.value ** 3 *\
        special.expi(helper1)
    return f1, f2
