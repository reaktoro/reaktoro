from PyReaktoro import *

import numpy

#==============================================================================
# Extending Reaktoro classes with handy methods 
#==============================================================================

#------------------------------------------------------------------------------
# Extending class ChemicalOutput
#------------------------------------------------------------------------------
# Define a method to convert the file data into an array.
ChemicalOutput.array = lambda self: numpy.loadtxt(self.filename(), skiprows=1)
