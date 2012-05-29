from tardis import constants
import numpy as np
import sqlite3
f2A_coefficient = (2 * np.pi * constants.e**2) / (constants.me * constants.c)
f2B_coefficient = (np.pi * constants.e**2) / (constants.me * constants.c)

class group_concat_intarray(object):
    def __init__(self):
        self.list = []
    
    def step(self, value):
        self.list.append(value)
        
    def finalize(self):
        return sqlite3.Binary(
            np.array(self.list, dtype=np.int64).tostring()
        )

class calculate_p_internal_down(object):
    def __init__(self):
        self.p_internals_down = []
    def step(self, wl, g_lower, g_upper, f_lu, energy_lower):
        p_internal_down = (wl * 1e-8)**-2  * \
            f2A_coefficient * \
            (float(g_lower) / g_upper) * \
            f_lu * (energy_lower / constants.erg2ev)
        self.p_internals_down.append(p_internal_down)
    
    def finalize(self):
        return sqlite3.Binary(
            np.array(self.p_internals_down, dtype=np.float64).tostring()
        )
        #for testing purposes only
        #return ','.join(map(str, self.p_internals_down))

class calculate_p_internal_up(object):
    def __init__(self):
        self.p_internals_up = []
        
    def step(self, wl, f_lu, energy_lower):
        p_internal_up = (wl * 1e-8) * f2B_coefficient * f_lu * energy_lower
        self.p_internals_up.append(p_internal_up)
    
    def finalize(self):
        return sqlite3.Binary(
            np.array(self.p_internals_up, dtype=np.float64).tostring()
        )
        #for testing purposes only
        #return ','.join(map(str, self.p_internals_up))


    
class calculate_p_emission_down(object):
    def __init__(self):
        self.p_emissions_down = []
    def step(self, wl, g_lower, g_upper, f_lu, energy_delta):
        p_emission_down = (wl * 1e-8)**-2  * \
            f2A_coefficient * \
            (float(g_lower) / g_upper) * \
            f_lu * (energy_delta / constants.erg2ev)
        self.p_emissions_down.append(p_emission_down)
    
    def finalize(self):
        return sqlite3.Binary(
            np.array(self.p_emissions_down, dtype=np.float64).tostring()
        )
        #for testing purposes only        
        #return ','.join(map(str, self.p_emissions_down))
