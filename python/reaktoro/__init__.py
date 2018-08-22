from reaktoro.PyReaktoro import *

import numpy

#==============================================================================
# Extending Reaktoro classes with handy methods 
#==============================================================================

#------------------------------------------------------------------------------
# Extending class ChemicalOutput
#------------------------------------------------------------------------------
# Define a method to convert the file data into an array.
ChemicalOutput.array = lambda self: numpy.loadtxt(self.filename(), skiprows=1)


# Define a method to convert the data in the file into a dictionary, with the 
# headings being the keys, and the column data as the values.
def _ChemicalOutput_dict(self):
    data = numpy.loadtxt(self.filename(), skiprows=1)
    columns = [ data[:, i] for i in xrange(data.shape[1]) ]
    return { heading:column for heading,column in zip(self.headings(), columns) } 

ChemicalOutput.dict = _ChemicalOutput_dict
