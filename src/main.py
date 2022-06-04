"""
@author: Daniel Elisabethsønn Antonsen, UiT Institute of Physics and Technology

Main file for programming assignment 4; Ion Chemistry
"""
# Importing libraries and modules
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
# Opening files
def openData(filename):
    """
    Function to load data from files
    """
    return np.loadtxt(os.path.join("resources", filename), comments="%")

# Loading data
msisFile = openData("MSIS.dat")
iriFile = openData("IRI.dat")

# Loading data for height and densities
height = msisFile[100:, 0:1] # Height, km
O2 = msisFile[100:, 3:4] * 1E6 # Number density molecular oxygen, m^-3
O = msisFile[100:, 1:2] * 1E6 # Number density atomic oxygen, m^-3
N2 = msisFile[100:, 2:3] * 1E6 # Number density molecular nitrogen, m^-3
ne = iriFile[:, 1:2] # Number density electrons, m^-3
Oplus = iriFile[:, 4:5] * ne/100 # Number density O ions, m^-3
O2plus = iriFile[:, 7:8] * ne/100 # Number density O2 ions, m^-3
NOplus = iriFile[:, 8:9] * ne/100 # Number density NO ions, m^-3
N2plus = 0 # Number density N2 ions, m^-3
NO = 0 # Number density NO, m^-3
ti = iriFile[:, 2:3] # Ion temperature, K
te = iriFile[:, 3:4] # Electron temperature, K
tn = msisFile[100:, 5:6] # Neutral temperature, K
tr = (tn + ti) / 2
# Defining time in seconds
time = np.arange(0, 4200, 1)

# 
def odes(x, t, h, changeT=False):
    """
    Function defining the ODEs
    """
    Te = te[int(h)]
    Tr = tr[int(h)]
    nO = O[int(h)]
    nO2 = O2[int(h)]
    nN2 = N2[int(h)]

    # Defining reactionrate (m^3/s)
    alpha1 = 2.1E-13 * (Te / 300)**(-0.85)
    alpha2 = 1.9E-13 * (Te / 300)**(-0.5)
    alpha3 = 1.8E-13 * (Te / 300)**(-0.39)
    alphar = 3.7E-18 * (250 / Te)**(0.7)

    k1 = 2E-18
    k2 = 2E-17 * (Tr / 300)**(-0.4)
    k3 = 4.4E-16
    k4 = 5E-22
    k5 = 1.4E-16 * (Tr / 300)**(-0.44)
    k6 = 5E-17 * (Tr / 300)**(-0.8)
    # Defining function to represent the ionization-rate
    def ion_rate(t):
        """
        Function to represent the ionization-rate over some time interval. (/m^3/s)
        R eturns the ionization-rate for the given interval
        """
        if t < 3600:
            return 1E8
        elif t >= 3600 and t <= 3700:
            return 1E10
        else:
            return 0

    qe = ion_rate(t)

    # Defining ionization-rate (/m^3/s)
    qN2plus = qe * (0.92 * nN2) / (0.92*nN2 + nO2 + 0.56*nO)
    qOplus = qe * (0.56 * nO) / (0.92*nN2 + nO2 + 0.56*nO)
    qO2plus = qe * (nO2) / (0.92*nN2 + nO2 + 0.56*nO)

    # Assigning each ODE to a vector element
    ne = x[0]
    nOplus = x[1]
    nO2plus = x[2]
    nN2plus = x[3]
    nNO = x[4]
    nNOplus = x[5]

    # Defining ODEs 
    dne_dt = qe - ne * (alpha1 * nNOplus + alpha2 * nO2plus + alpha3 * nN2plus + alphar * nOplus)
    dOplus_dt = qOplus - nOplus * (k1 * nN2 + k2 * nO2 + alphar * ne)
    dO2plus_dt = qO2plus + (k2 * nOplus * nO2 + k6 * nN2plus * nO2 - k3 * nNO * nO2plus - k4 * nN2 * nO2plus - alpha2 * nO2plus * ne)
    dN2plus_dt = qN2plus - nN2plus * (alpha3 * ne + k5 * nO + k6 * nO2)
    dNO_dt = (nO2plus * nN2 * k4) - (k3 * nNO * nO2plus)
    dNOplus_dt = (k1 * nOplus * nN2) + (k3 * nNO * nO2plus) + (k4 * nO2plus * nN2) + (k5 * nN2plus * nO) - (nNOplus * alpha1 * ne)

    return np.array([dne_dt, dOplus_dt, dO2plus_dt, dN2plus_dt, dNO_dt, dNOplus_dt]).squeeze()


# Functions to return initial values
def initialvalues(index):
    """
    Initial values at a given height
    """
    return np.array([ne[index], Oplus[index], O2plus[index], N2plus, NO, NOplus[index]], dtype=np.float64)

def solveODEs(h, t):
    """
    Solving coupled ODEs using odeint function from scipy
    Returning solutions in form of an array
    """
    iv = initialvalues(h)
    solveODE = odeint(odes, iv, t, args=(h,))
    return solveODE

# Running function for three different heights 
H110km = solveODEs(10, time[:3600])
H170km = solveODEs(70, time[:3600])
H230km = solveODEs(130, time[:3600])

# Plotting functions from 0 to 3600s 
fig, ax = plt.subplots(1, 3, sharey=True)
# Plotting for height 110km
ax[0].plot(time[:3600], H110km[:, :])
ax[0].set_ylabel(r"Density [$m^{-3}$]")
ax[0].set_xlabel(r"Time [s]")
ax[0].set_title(r"Height 110km")
ax[0].legend([r"$e^-$", r"$O^{+}$", r"$O_{2}^{+}$", r"$N_{2}^{+}$", r"$NO$", r"$NO^{+}$"])
# Plotting for height 170km
ax[1].plot(time[:3600], H170km[:, :])
ax[1].set_xlabel(r"Time [s]")
ax[1].set_title(r"Height 170km")
ax[1].legend([r"$e^-$", r"$O^{+}$", r"$O_{2}^{+}$", r"$N_{2}^{+}$", r"$NO$", r"$NO^{+}$"])
# Plotting for height 230km
ax[2].plot(time[:3600], H230km[:, :])
ax[2].set_xlabel(r"Time [s]")
ax[2].set_title(r"Height 230km")
ax[2].legend([r"$e^-$", r"$O^{+}$", r"$O_{2}^{+}$", r"$N_{2}^{+}$", r"$NO$", r"$NO^{+}$"])
# Figure layout
fig.suptitle(r"Densities for constant ionization-rate of $1\cdot 10^{8} (/m^{3}/s)$")
fig.tight_layout()

# 
ionrateH110km = solveODEs(10, time[3600:])
ionrateH170km = solveODEs(70, time[3600:])
ionrateH230km = solveODEs(130, time[3600:])

# Plotting for variable ionization-rate from 3600s to 4200s
fig1, ax1 = plt.subplots(1, 3, sharey=True)
# Plotting for height 110km
ax1[0].plot(time[3600:], ionrateH110km[:, :])
ax1[0].set_xlabel(r"Time [s]")
ax1[0].set_ylabel(r"Density [$m^{-3}$]")
ax1[0].set_title(r"Height 110km")
ax1[0].legend([r"$e^-$", r"$O^{+}$", r"$O_{2}^{+}$", r"$N_{2}^{+}$", r"$NO$", r"$NO^{+}$"])
# Plotting for height 170km 
ax1[1].plot(time[3600:], ionrateH170km[:, :])
ax1[1].set_xlabel(r"Time [s]")
ax1[1].set_title(r"Height 170km")
ax1[1].legend([r"$e^-$", r"$O^{+}$", r"$O_{2}^{+}$", r"$N_{2}^{+}$", r"$NO$", r"$NO^{+}$"])
# Plotting for height 230km
ax1[2].plot(time[3600:], ionrateH230km[:, :])
ax1[2].set_xlabel(r"Time [s]")
ax1[2].set_title(r"Height 230km")
ax1[2].legend([r"$e^-$", r"$O^{+}$", r"$O_{2}^{+}$", r"$N_{2}^{+}$", r"$NO$", r"$NO^{+}$"])
# Figure layout
fig1.suptitle("Densities for changing ionization-rate")
fig1.tight_layout()

if __name__ == "__main__":
    plt.show()

