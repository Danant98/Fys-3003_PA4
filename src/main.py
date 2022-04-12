"""
@author: Daniel Elisabethsønn Antonsen, UiT Institute of Physics and Technology

Main file for assignment
"""
# Importing libraries and modules
import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
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
ne = iriFile[:, 1:2] # Electron number density, m^(-3)
nO2 = msisFile[:, 3:4] * (1/1E6) # Number density of molecular oxygen, m^(-3)
nO = msisFile[:, 1:2] * (1/1E6) # Number density of atomic oxygen, m^(-3)
nN2 = msisFile[:, 2:3] * (1/1E6) # Number density of molecular nitrogen, m^(-3)
temperature = msisFile[:, 5:6] # Temperature, K

# Defining the time using numpy array
t = np.arange(0, 3600, 1) # Time, s
 # Radiative recombination-rate, (m^3)/s
radiativeRecRate = 3.7E-18 * (250 / temperature)**(0.7)
# Defining reaction rates, (m^3)/s
alpha1 = 2.1E-13 * (temperature / 300)**(-0.85)
alpha2 = 1.9E-13 * (temperature / 300)**(-0.5)
alpha3 = 1.8E-13 * (temperature / 300)**(-0.39)






if __name__ == '__main__':
    #plt.show()
    pass
