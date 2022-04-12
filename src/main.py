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
height = msisFile[:, 0:1] # Height, km
ne = iriFile[:, 1:2] # Electron density, m^(-3)
nO2 = msisFile[:, 3:4] # Number density of molecular oxygen, cm^(-3)
nO = msisFile[:, 1:2] # Number density of atomic oxygen, cm^(-3)
nN2 = msisFile[:, 2:3] # Number density of molecular nitrogen, cm^(-3)
temperature = msisFile[:, 5:6] # Temperature, K




if __name__ == '__main__':
    print(temperature)
