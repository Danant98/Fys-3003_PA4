"""
@author: Daniel Elisabethsønn Antonsen, UiT Institute of Physics and Technology

Main file for assignment
"""
# Importing libraries and modules
import os
import numpy as np
# Opening files
def openData(filename):
    """
    Function to load in data from file
    """
    return np.loadtxt(os.path.join("resources", filename), comments="%", dtype=np.float32)
# Loading data
msisFile = openData("MSIS.dat")
iriFile = openData("IRI.dat")
# Loading data for height and densities
height = msisFile[:, 0:1]
nO2 = None
nO = None
nN2 = None
ne = iriFile[:, 1:2]




