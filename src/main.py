"""
@author: Daniel Elisabethsønn Antonsen, UiT Institute of Physics and Technology

Main file for assignment
"""
# Importing libraries and modules
import os
import numpy as np
import scipy as sp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
# Defining functions to open files containing the data
def openData(filename):
    """
    Function to load data from files
    """
    return (np.loadtxt(os.path.join("resources", filename), comments="%", dtype=np.float32))
# Loading data
msisFile = openData("MSIS.dat")
iriFile = openData("IRI.dat")
# Loading data for height and densities
height = msisFile[100:, 0:1] # Height, km
ne = iriFile[:, 1:2] # Electron number density, m^(-3)
nO2 = msisFile[100:, 3:4] * (1/1E6) # Number density of molecular oxygen, m^(-3)
nO = msisFile[100:, 1:2] * (1/1E6) # Number density of atomic oxygen, m^(-3)
nN2 = msisFile[100:, 2:3] * (1/1E6) # Number density of molecular nitrogen, m^(-3)
temperature = msisFile[100:, 5:6] # Neutral temperature, K
ionTemp = iriFile[:, 2:3] # Ion temperature, K
electronTemp = iriFile[:, 3:4] # Electron temperature, K
# Loading the amount of ions
ionO = iriFile[:, 4:5]
ionH = iriFile[:, 5:6]
ionHe = iriFile[:, 6:7]
ionO2 = iriFile[:, 7:8]
ionNO = iriFile[:, 8:9]
ionN = iriFile[:, 10:11]
# Defining the time as a numpy array with spacing 1s
t = np.linspace(0, 3600, 1) # Time, s
# Defining the constant ionization-rate
ionRate = 1E8
# Integrating 
def coupledODEforE(initialValues, t):
    """
    Creating the coupled ODE for E-region     
    """



    return None


"""
Solving the continuity equations for two heights; approx 110km and 230km. Both integrated over 3600s.
Is to be used as a stable background for initial conditions for the different densities.
""" 


# Running only of this file is runned
if __name__ == '__main__':
    plt.show()
