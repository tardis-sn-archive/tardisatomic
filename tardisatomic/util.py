from astropy import units as u


def convert_vacuum_to_air(wavelength_vacuum):
    sigma2 = (1e4 * u.angstrom / wavelength_vacuum)**2.
    fact = 1. +  5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/( 57.362 - sigma2)

    return wavelength_vacuum / fact

def convert_air_to_vacuum(wavelength_air):
    sigma2 = (1e4/wavelength_air)**2.
    fact = 1. +  5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/( 57.362 - sigma2)
    return wavelength_air * fact