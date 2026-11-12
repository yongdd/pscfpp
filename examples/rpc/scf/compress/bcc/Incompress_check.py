#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: kexinchen
"""

from pscfpp.field import Field

# Read the r-grid concentration file. The number of header lines is not
# fixed (a field without a space group, such as the one written after
# COMPRESS, has no group_name entry), so parse the header rather than
# skip a fixed number of lines.
field = Field('out/c_rf')

#concentration matrix
c = field.data

# Check that a row was read for every grid point
mesh = field.header['mesh'] if 'mesh' in field.header else field.header['ngrid']
nGrid = 1
for n in mesh:
    nGrid *= n
if len(c) != nGrid:
    raise Exception('Read ' + str(len(c)) + ' rows for a mesh of '
                    + str(nGrid) + ' grid points')

#Cound the summation of concentration at each grid differ from 1
diff = 0
for i in range(len(c)):
    sum_c = 0
    for j in range(len(c[i])):
        #sum the the concentation at grid i
        sum_c += c[i][j]
    diff += abs(sum_c -1)

print("The sum of inaccuracy of concentration on each grid is ")
print(diff)
