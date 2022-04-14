"""
@author: Daniel Elisabethsønn Antonsen, UiT Institute of Physics and Technology

Main file for assignment
"""
# Importing libraries and modules
import os
import numpy as np
import scipy as sp
import scipy.integrate as integral 
import matplotlib.pyplot as plt
# Defining functions to open files containing the data
def openData(filename):
    """
    Function to load in data from file
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
# Defining the time using numpy array
t = np.arange(0, 3600, 1) # Time, s
# Defining the constant ionization-rate
ionRate = 1E8
# Integrating 
def solveDiff(h, time):
    """
    Function to solve the differential equation for a given height 'h' and a given time-interval 'time'

    Input 'h' with reference height at 100 km and 'time'
    """
    height = height[h]
     # Radiative recombination-rate, (m^3)/s
    radiativeRecRate = 3.7E-18 * (250 / electronTemp[height])**(0.7)
    # Defining dissociate recombination reaction rates, (m^3)/s
    alpha1 = 2.1E-13 * (electronTemp[height] / 300)**(-0.85)
    alpha2 = 1.9E-13 * (electronTemp[height] / 300)**(-0.5)
    alpha3 = 1.8E-13 * (electronTemp[height] / 300)**(-0.39)
    # Defining the average alpha
    avgAlpha = (alpha1 * (ionNO / ne) + alpha2 * (ionO2 / ne) + alpha3 * (nN2 / ne))
    # Integrating
    integralNee = integral.cumulative_trapezoid(ionRate - avgAlpha * (ne)(ne), time)
    return (integralNee)





if __name__ == '__main__':
    #plt.show()
    pass
