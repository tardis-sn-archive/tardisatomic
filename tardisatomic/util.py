from astropy import units as u

def convert_vacuum_to_air(wavelength_vacuum):
    sigma2 = (1e4 / wavelength_vacuum)**2.
    fact = 1. +  5.792105e-2 / (238.0185 - sigma2) + 1.67917e-3 / ( 57.362 - sigma2)

    return wavelength_vacuum / fact

def convert_air_to_vacuum(wavelength_air):
    sigma2 = (1e4/wavelength_air)**2.
    fact = 1. +  5.792105e-2 / (238.0185 - sigma2) + 1.67917e-3 / ( 57.362 - sigma2)

    return wavelength_air * fact

"""
>> print make_table([['Name', 'Favorite Food', 'Favorite Subject'],
                     ['Joe', 'Hamburgers', 'Cars'],
                     ['Jill', 'Salads', 'American Idol'],
                     ['Sally', 'Tofu', 'Math']])

+------------------+------------------+------------------+
| Name             | Favorite Food    | Favorite Subject |
+==================+==================+==================+
| Joe              | Hamburgers       | Cars             |
+------------------+------------------+------------------+
| Jill             | Salads           | American Idol    |
+------------------+------------------+------------------+
| Sally            | Tofu             | Math             |
+------------------+------------------+------------------+
by cieplak@stackoverflow
"""

def make_table(grid):
    cell_width = 2 + max(reduce(lambda x,y: x+y, [[len(item) for item in row] for row in grid], []))
    num_cols = len(grid[0])
    rst = table_div(num_cols, cell_width, 0)
    header_flag = 1
    for row in grid:
        rst = rst + '| ' + '| '.join([normalize_cell(x, cell_width-1) for x in row]) + '|\n'
        rst = rst + table_div(num_cols, cell_width, header_flag)
        header_flag = 0
    return rst

def table_div(num_cols, col_width, header_flag):
    if header_flag == 1:
        return num_cols*('+' + (col_width)*'=') + '+\n'
    else:
        return num_cols*('+' + (col_width)*'-') + '+\n'

def normalize_cell(string, length):
    return string + ((length - len(string)) * ' ')