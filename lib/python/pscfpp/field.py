"""! Module for parsing field files. """

##
# Container for data in a PSCF field file.
#
#  A Field object can parse a PSCF "field file" and store its contents.
#  A field file contains both a header section and a data section. 
#  This class can read and parse all of the field file formats used by
#  PSCF programs including the basis, rgrid, and kgrid formats used by
#  the pscf_rpc and pscf_rpg programs to represent periodic fields, and
#  the simpler file format used by pscf_1d for one-dimensional fields.
#  For detail on these file formats see the section of the web manual
#  on \ref user_field_page "field files". Users can access and modify
#  stored values of the parameters in the header and values in the data 
#  section after parsing by using specific statements (commands), and 
#  can write the entire object to a file in proper format.
#
#  Every PSCF field file format contains a header section and a data 
#  section. The header contains globale information about the system
#  and spatial discretization. The data section is a table of field
#  components in which each row contains components associated with
#  a single spatial degree of freedom, and each column contains 
#  components associated with a single monomer type. The format 
#  of header is different for different types of field file, which 
#  makes it possible for the Field class to infer information about
#  the type of a field file while parsing the header,
#  
#  A Field object has three attributes named 'header', 'data', and 
#  'type'. The 'header' attribute is a dictionary that contains the
#  contents of the header section, using the variable names or labels
#  as dictionary keys. The 'data' attribute is a list of lists in 
#  which each element of the outer list corresponds to one row of
#  the data section in the field file. When accessing elements of 
#  the data attribute using two subscripts, the first subscript (the 
#  row index) is thus index for a spatial degree of freedom (i.e., 
#  a grid point or basis function) and the second (the column index) 
#  is an index for a monomer type. 
#
#  The 'type' attribute is a string with possible values 'fd1d', 
#  'basis', 'rgrid' or 'kgrid'. 
#
#  The listings of waves and stars written by the WRITE_WAVES and
#  WRITE_STARS commands of pscf_rpc and pscf_rpg use the header of a
#  periodic field file (without N_monomer) followed by a table with one
#  row per wave or star. They can also be read, and are given the types
#  'waves' and 'stars'. The space group files written by WRITE_GROUP
#  (the format of the files in data/groups) are not field files and are
#  rejected with an exception.
#
#  **Construction:**
#
#  The Field constructor can parse any PSCF field file, and infer its 
#  type. To read and parse a field file with name 'field', one could 
#  enter:
#  \code
#     from pscfpp.field import *
#     f = Field('field')
#  \endcode
#
#  **Printing and writing to file:**
#
#  The __str__ function can be used to convert a Field object to a
#  multi-line string with a consistent with that of the file from which
#  it was created. The write function can be used to write that string
#  to a specified file. To write a field stored in a Field object named
#  f to a file named 'out', one would enter:
#  \code
#     f.write('out')
#  \endcode
#
#  **Accessing elements:**
#
#  The header attribute is a Python dictionary containing
#  the contents of the header section of the field file. One can
#  access elements using names of header variables as dictionary
#  keys, as in the expressions
#  \code
#     f.header['format']    
#     f.header['N_monomer']
#  \endcode
#  The available keys are different for different file formats.
#  Values for some quantities that are generally represented by more 
#  than one value are stored as lists, such as the list 'cell_param' 
#  of cell parameter values to describe the unit cell in  a periodic 
#  system.
#
#  The data attribute is a list of lists, in which elements can be 
#  accessed using square brackets with integer indices.  The expressions
#  \code
#    f.data[10]
#    f.data[10][1]
#  \endcode
#  access an entire row of the data section or a single component,
#  respectively.
# 
#  Elements of the header and data attributes can also be modified,
#  simply by putting corresponding expressions on the left side of 
#  an equality (assignment) operator. For example, the expressions
#  \code
#     f.header['dim'] = 3  
#     f.header['cell_param'] = [2.45]   
#     f.header['group_name'] = 'P_1'
#  \endcode
#  modify variables defined in the header of a field file for a 
#  periodic system.
#
#  **Modifying data table structure:**
#
#  The methods addColumn, deleteColumn, and  reorder can change the 
#  structure of the data section by adding, deleting or reordering 
#  columns associated with different monomer types. 
#
class Field:

   ##
   # Constructor.
   #
   # \param filename  a filename string
   #
   def __init__(self, filename):
      self.data = []
      self.header = {}
      self.type = ''
      with open(filename) as f:
         self.read(f)

   ##
   # Read and parse the passed-in file
   # 
   # This function reads the passed-in open-file object line
   # by line and update the read items into instance
   # variables. The reading stops when all lines in file are
   # read.
   #
   # \param openFile  a file object, open for reading
   #
   def read(self, openFile):
      line = openFile.readline()
      l = line.split()
      if l[0] == 'format':
         self.header['format'] = [getValue(l[1]), getValue(l[2])]
      elif l[0] == 'nx':
         self.type ='fd1d'
         self.header['nx'] = getValue(l[1])
      elif l[0] == 'dim':
         name = str(getattr(openFile, 'name', 'this file'))
         raise Exception('Not a field file: ' + name
                         + ' is a space group file (as written by'
                         + ' WRITE_GROUP), which the Field class'
                         + ' does not read.')
      else:
         raise Exception('Not valid field file.')

      if self.type != 'fd1d':
         # Label and value lines alternate. The header ends with the
         # mesh dimensions (r-grid and k-grid files) or the number of
         # basis functions (basis files); group_name is optional, so
         # the number of header lines is not fixed.
         while True:
            line = openFile.readline()
            if line == '':
               raise Exception('Incomplete field file header.')
            l = line.split()
            if len(l) == 0:
               continue
            name = l[0]
            line = openFile.readline()
            l = line.split()
            if name == 'mesh' or name == 'ngrid' or name == 'cell_param':
               d = []
               for i in range(0, len(l)):
                  d.append(getValue(l[i]))
               self.header[name] = d
            else:
               self.header[name] = getValue(l[0])
            if name == 'N_basis' or name == 'N_star':
               # Basis field file, or listing of stars (WRITE_STARS);
               # the two are distinguished by the first data row
               self.type = 'basis'
               break
            if name == 'N_wave':
               # Listing of waves (WRITE_WAVES)
               self.type = 'waves'
               break
            if name == 'mesh' or name == 'ngrid':
               break
         # First data line
         line = openFile.readline()
         l = line.split()
      else:
         # Each header line contains a label and a value. The number of
         # header lines is not fixed: w and c field files have nx and
         # nm, files written by WRITE_BLOCK_C have nx, n_block and
         # n_solvent, and propagator slice files have only nx. The data
         # section starts with the first line that begins with a number.
         while True:
            line = openFile.readline()
            l = line.split()
            if line == '':
               break
            if len(l) == 0:
               continue
            if not isinstance(getValue(l[0]), str):
               break
            if len(l) != 2:
               raise Exception('Invalid header line in field file: '
                               + line.strip())
            self.header[l[0]] = getValue(l[1])

      while line != '':
         if l:
            self.data.append([getValue(v) for v in l])
         line = openFile.readline()
         l = line.split()
         

      if self.type == 'basis':
         # Integer-valued fields are valid basis files too. Recognise
         # star listings by their column count, preserving a basis file
         # when its N_monomer makes both interpretations possible.
         starColumns = self.header['dim'] + 6
         if 'N_star' in self.header:
            starColumns += 1
         basisColumns = self.header.get('N_monomer', -1) + self.header['dim'] + 1
         if (self.data and len(self.data[0]) == starColumns
             and len(self.data[0]) != basisColumns
             and all(isinstance(x, int) for x in self.data[0])):
            self.type = 'stars'

      if self.type == '':
         # A mesh header is shared by real and Fourier grid formats.
         # Numeric spelling is not a format marker: real coefficients
         # may be integer literals, and Fourier indices may carry '+'.
         nMonomer = self.header.get('N_monomer', 0)
         if not self.data or nMonomer <= 0:
            raise Exception('Missing mesh field data or N_monomer.')
         nColumn = len(self.data[0])
         if nColumn == nMonomer:
            self.type = 'rgrid'
         elif nColumn == 1 + 2*nMonomer:
            self.type = 'kgrid'
         else:
            raise Exception('Invalid number of columns in mesh field data.')
         for row in self.data:
            if len(row) != nColumn:
               raise Exception('Inconsistent number of mesh field columns.')
            if self.type == 'kgrid' and not isinstance(row[0], int):
               raise Exception('A k-grid row must begin with an integer index.')

   ##
   # Return string representation of this Field.
   #
   # This function return the string representation of the Field object 
   # in the appropriate field file format.
   #
   def __str__(self):
      out = ''
      if self.type == 'fd1d':
         # Width of the label column: 7, or more for longer labels
         width = 7
         for x in self.header:
            width = max(width, len(x) + 3)
      for x, y in self.header.items():
         if self.type == 'fd1d':
            o = f'{x:<{width}}' + str(y) + '\n'
         else:
            if x == 'format':
               o = x + '   ' + str(y[0]) + '   ' + str(y[1]) + '\n'
            else:
               o = x + '\n'
               if x == 'cell_param':
                  val = f'{y[0]:.10e}'
                  o += f'{val:>20}'
                  if len(y) == 1:
                     o += '\n'
                  for i in range(1, len(y)):
                     val = f'{y[i]:.10e}'
                     o += f'{val:>18}'
                     if i == len(y)-1:
                        o += '\n'
               elif x == 'mesh' or x == 'ngrid':
                  o += f'{y[0]:>20}'
                  if len(y) == 1:
                     o += '\n'
                  for i in range(1, len(y)):
                     o += f'{y[i]:>8}'
                     if i == len(y)-1:
                        o += '\n'
               else:
                  o += f'{str(y):>20s}' + '\n'

         out += o

      for i in range(0, len(self.data)):
         row = ''
         for j in range(0, len(self.data[i])):
            if self.type == 'fd1d':
               # Rows start with a grid point index, except in files
               # with a single field and no index column
               if j == 0 and isinstance(self.data[i][j], int):
                  row += f'{self.data[i][j]:>6}'
               else:
                  val = f'{self.data[i][j]:.11e}'
                  row += f'{val:>20}'

            if self.type == 'basis':
               if j < self.header['N_monomer']:
                  val = f'{self.data[i][j]:.10e}'
                  row += f'{val:>20}'
               else:
                  # Each integer column is preceded by a space, so that
                  # columns stay separated when a value fills its width
                  if j == self.header['N_monomer']:
                     row += '   '
                  row += ' ' + f'{self.data[i][j]:>4}'

            if self.type == 'rgrid':
               val = f'{self.data[i][j]:.15e}'
               row += ' ' + f'{val:>22}'

            if self.type == 'waves':
               # Columns: i, waveId, D wave indices, starId, Re(c), Im(c)
               nc = len(self.data[i])
               if j < 2:
                  row += ' ' + f'{self.data[i][j]:>7}'
               elif j < nc - 3:
                  row += ' ' + f'{self.data[i][j]:>4}'
               elif j == nc - 3:
                  row += ' ' + f'{self.data[i][j]:>5}'
               else:
                  row += '  ' + f'{self.data[i][j]:>15.7e}'

            if self.type == 'stars':
               # Columns: basisId, starId, size, beginId, endId,
               # starInvert, cancel (N_star files only), D wave indices
               widths = [5, 5, 4, 7, 7, 3]
               if 'N_star' in self.header:
                  widths.append(3)
               if j < len(widths):
                  w = widths[j]
               else:
                  w = 5
               row += ' ' + f'{self.data[i][j]:>{w}}'

            if self.type == 'kgrid':
               if j == 0:
                  row += f'{self.data[i][j]:>6}'
               else:
                  val = f'{self.data[i][j]:.12e}'
                  if j%2 != 0:
                     row += ' ' + f'{val:>21}'
                  else:
                     row += ' ' + f'{val:>19}'

         out += row + '\n'

      return out

   ##
   # Write out field to a file.
   #
   # This function writes out the field file string to the 
   # specified file with the name of the passed-in parameter, 
   # filename.
   #
   # \param filename  a filename string.
   #
   def write(self, filename):
      with open(filename, 'w') as f:
         f.write(self.__str__())

   ##
   # Add a new column to the data list.
   #
   # This function adds a new column to the data list of this Field
   # object. The parameter 'index' is the desired position of the new
   # column. The parameter 'element' can be either a single value that
   # is applied to every element of the entire column, or a list of
   # of values that represents the entire column of new values. If a
   # list of values is supplied, the length of the list must be equal
   # to the number of rows in the data list. 
   #
   # Example:
   # \code
   #    f.addColumn(2, 1.2)    
   #    f.addColumn(2, [list of values for the whole column])
   # \endcode
   #
   # \param index  integer column index 
   # \param element  value (float) or values (list) for elements of column
   #
   def addColumn(self, index, element):
      if isinstance(element, list):
         if not (len(element) == len(self.data)):
             raise Exception('Incorrect dimension for parameter element')
         for i in range(0, len(self.data)):
            self.data[i].insert(index, element[i])
      else:
         for i in range(0, len(self.data)):
            self.data[i].insert(index, element)

   ##
   # Delete an existing column from the data list.
   #
   # This function deletes an existing column from the data list of 
   # this Field object. Input parameter "index" is the index of the
   # column that should removed.
   #
   # \param index  column index, numbered from 0
   #
   def deleteColumn(self, index):
      for i in range(0, len(self.data)):
         self.data[i].pop(index)

   ##
   # Reorder the data list.
   #
   # This function reorders the data list of the Field object. By
   # passing in one parameter, order, an integer list that represents
   # the new order of the columns in the data list, the data list of
   # the Field object can be updated with the desired order. The 
   # length of the list 'order' must be equal to the total number 
   # of columns in the data section, including any that do not 
   # contain field component values. 
   #
   # Note: This function treats all columns in the file format 
   # equivalently, whether they contain field component values or
   # other information. Specifically, when treating a field file 
   # in basis file format for a system with C monomer types in a 
   # system of spatial dimension D, the field file format contains 
   # C + D + 1 columns, in which only the first C contain field 
   # components and the remaining D+1 contain Miller indices for
   # a characteristic wave of each star and the number of waves in
   # the star. To re-order columns that contain field components
   # one must enter a permutation of the integers [0,...,C+D] in
   # which only the first C columns (columsn 0,...,C-1) are 
   # re-ordered and the remaining D+1 integers must appear in 
   # their original order.
   #
   # Example: To re-order the field component columns for a 3D periodic
   # system containing three monomer, one might use something like
   # \code
   #     f.reorder([2, 0, 1, 3, 4, 5, 6])
   # \endcode
   # Here, columns 0-2 represent field components, which have been
   # re-ordered, while columns 3-6 represent information about basis
   # functions (Miller indices and the number of waves per star)
   # which are left in their original order.
   #
   # \param order  a list of position integers.
   #
   def reorder(self, order):
      nc = len(order)
      newData = []
      for i in range(0, len(self.data)):
         if not (nc == len(self.data[i])):
            raise Exception('Incorrect length for permutation order')
         d = []
         for j in range(0,len(order)):
            d.append(self.data[i][order[j]])
         newData.append(d)
      self.data = newData

##
# Distinguish the correct type of the value from a string.
#
# This function is a helper function to read a open-file
# object. It distinguishes the correct type of the value
# from a string to either integer, float or string.
#
# Return value:
#
# The exact value with the correct type.
#
# \param v  a string represents a value.
#
def getValue(v):
   if v.isdigit() or (v[:1] in ('-', '+') and v[1:].isdigit()):
      val = int(v)  
   else:
      try:
         val = float(v)
      except ValueError:
         val = str(v)
   
   return val

