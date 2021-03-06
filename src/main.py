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

# Defining time in seconds
time = np.arange(0, 4200, 1)
# Decay time
decay_time = np.arange(0, 500, 1)

# Defining a function containg our ODEs
def odes(x, t, h, htemp=False, sinusoidal=False):
    """
    Function defining the ODEs
    """
    # Making if statement to check if we would change the temperature at a given altitude
    if htemp:
        if h < 50:
            Te = te[int(h)] + 1000
            Tr = (tn[int(h)] + ti[int(h)] + 1000) / 2
        elif h > 50:
            Te = te[int(h)] + 2000
            Tr = (tn[int(h)] + ti[int(h)] + 2000) / 2
    else:
        Te = te[int(h)]
        Tr = (tn[int(h)] + ti[int(h)]) / 2
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

    # Defining function to represent the ionization-rate for electrons
    def ion_rate(t):
        """
        Function to represent the ionization-rate over some time interval. (/m^3/s)
        Returns the ionization-rate for the given interval
        """
        if t < 3600:
            return 1E8
        elif t >= 3600 and t <= 3700:
            return 1E10
        else:
            return 0
    
    def sinusoidal_ionrate(t):
        """
        Function to represent the ionization-rate for a sinusoidal ionization-rate. (/m^3/s)
        Returns the ionization-rate for the given interval
        """
        qehat = 2E10
        if t <= 3600:
            return 1E8
        elif t > 3600 and t <= 3700:
            return qehat*(np.sin(2*np.pi*t/20))**2
        else:
            return 0

    # Changing ionozation-rate 
    if sinusoidal:
        qe = sinusoidal_ionrate(t)
    else:
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

def solveODEs(h, t, htemp=False, sinusoidal=False):
    """
    Solving coupled ODEs using odeint function from scipy
    Returning solutions in form of an array
    """
    iv = initialvalues(h)
    solveODE = odeint(odes, iv, t, args=(h, htemp, sinusoidal))
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

# Changing ionozation-rate
ionrateH110km = solveODEs(10, time)
ionrateH170km = solveODEs(70, time)
ionrateH230km = solveODEs(130, time)

# Plotting for variable ionization-rate from 3600s to 4200s
fig1, ax1 = plt.subplots(1, 3, sharey=True)
# Plotting for height 110km
ax1[0].plot(time[3600:], ionrateH110km[3600:, :])
ax1[0].set_xlabel(r"Time [s]")
ax1[0].set_ylabel(r"Density [$m^{-3}$]")
ax1[0].set_title(r"Height 110km")
ax1[0].legend([r"$e^-$", r"$O^{+}$", r"$O_{2}^{+}$", r"$N_{2}^{+}$", r"$NO$", r"$NO^{+}$"])
# Plotting for height 170km 
ax1[1].plot(time[3600:], ionrateH170km[3600:, :])
ax1[1].set_xlabel(r"Time [s]")
ax1[1].set_title(r"Height 170km")
ax1[1].legend([r"$e^-$", r"$O^{+}$", r"$O_{2}^{+}$", r"$N_{2}^{+}$", r"$NO$", r"$NO^{+}$"])
# Plotting for height 230km
ax1[2].plot(time[3600:], ionrateH230km[3600:, :])
ax1[2].set_xlabel(r"Time [s]")
ax1[2].set_title(r"Height 230km")
ax1[2].legend([r"$e^-$", r"$O^{+}$", r"$O_{2}^{+}$", r"$N_{2}^{+}$", r"$NO$", r"$NO^{+}$"])
# Figure layout
fig1.suptitle("Densities for changing ionization-rate")
fig1.tight_layout()

# Higher electron and ion temperature
hightTemp110km = solveODEs(10, time, True)
hightTemp170km = solveODEs(70, time, True)
hightTemp230km = solveODEs(130, time, True)

# Plotting the densites for higher electron and ion tempperature
fig2, ax2 = plt.subplots(1, 3, sharey=True)
# Plotting for height 110km
ax2[0].plot(time[3600:], hightTemp110km[3600:, :])
ax2[0].set_xlabel(r"Time [s]")
ax2[0].set_ylabel(r"Density [$m^{-3}$]")
ax2[0].set_title(r"Height 110km")
ax2[0].legend([r"$e^-$", r"$O^{+}$", r"$O_{2}^{+}$", r"$N_{2}^{+}$", r"$NO$", r"$NO^{+}$"])
# Plotting for height 170km
ax2[1].plot(time[3600:], hightTemp170km[3600:, :])
ax2[1].set_xlabel(r"Time [s]")
ax2[1].set_title(r"Height 170km")
ax2[1].legend([r"$e^-$", r"$O^{+}$", r"$O_{2}^{+}$", r"$N_{2}^{+}$", r"$NO$", r"$NO^{+}$"])
# Plotting for height 230km
ax2[2].plot(time[3600:], hightTemp230km[3600:, :])
ax2[2].set_xlabel(r"Time [s]")
ax2[2].set_title(r"Height 230km")
ax2[2].legend([r"$e^-$", r"$O^{+}$", r"$O_{2}^{+}$", r"$N_{2}^{+}$", r"$NO$", r"$NO^{+}$"])
# Figure layout
fig2.suptitle(r"Densities for changing ionization-rate, higher e- and i-temperature")
fig2.tight_layout()

# Defining functions for alpha and beta decay
def alpha_decay(x, h, t):
    """
    Function to calculate the alpha-decay of electron density
    """
    # Finding initial values for a given height
    Te = te[int(h)]
    NE = x[3700, 0]
    noplus = x[3700, 5]
    o2plus = x[3700, 2]
    n2plus = x[3700, 3]
    # Defining varables for alphabar
    alpha1 = 2.1E-13 * (Te / 300)**(-0.85)
    alpha2 = 1.9E-13 * (Te / 300)**(-0.5)
    alpha3 = 1.8E-13 * (Te / 300)**(-0.39)
    alphabar = alpha1*(noplus/NE) + alpha2*(o2plus/NE) + alpha3*(n2plus/NE)
    # Expression for ne as a function of time
    E = NE/(1 + alphabar*NE*t)
    return E

def beta_decay(x, h, t):
    """
    Function to calculate the beta decay for electron density
    """
    Tr = (tn[int(h)] + ti[int(h)]) / 2
    o2 = O2[int(h)]
    k2 = 2E-17 * (Tr / 300)**(-0.4)
    NE = x[3700, 0]
    # Defining beta
    beta = k2*o2
    # Function for electron density 
    func = NE*np.exp(-beta*t)
    return func

# Running function for alpha and beta decay
alphaH110km = alpha_decay(ionrateH110km, 10, decay_time)
alphaH230km = alpha_decay(ionrateH230km, 130, decay_time)
betaH110km = beta_decay(ionrateH110km, 10, decay_time)
betaH230km = beta_decay(ionrateH230km, 130, decay_time)

# Plotting the alpha and beta decay with the electron density decay (Linear)
fig3, ax3 = plt.subplots(1, 2, sharey=True)
# Plotting for height 110km
ax3[0].plot(time[3700:], alphaH110km)
ax3[0].plot(time[3700:], betaH110km)
ax3[0].plot(time[3700:], ionrateH110km[3700:, 0], linestyle='--', color='red')
ax3[0].set_xlabel(r"Time [s]")
ax3[0].set_ylabel(r"Electron density $n_e$ [$m^{-3}$]")
ax3[0].set_title(r"Height 110km")
ax3[0].legend(["alpha decay", "beta decay", "regular decay"])
# Plotting for height 230km
ax3[1].plot(time[3700:], alphaH230km)
ax3[1].plot(time[3700:], betaH230km)
ax3[1].plot(time[3700:], ionrateH230km[3700:, 0], linestyle='--', color='red')
ax3[1].set_xlabel(r"Time [s]")
ax3[1].set_title(r"Height 230km")
ax3[1].legend(["alpha decay", "beta decay", "regular decay"])
# Format figure
fig3.suptitle(r"Comparison of alpha, beta and found decay for electron density (Linear)")
fig.tight_layout()

# Plotting the alpha and beta decay with the electron density decay (Semilog)
figlog, axlog = plt.subplots(1, 2, sharey=True)
# Plotting for height 110km
axlog[0].semilogy(time[3700:], alphaH110km)
axlog[0].semilogy(time[3700:], betaH110km)
axlog[0].semilogy(time[3700:], ionrateH110km[3700:, 0], linestyle='--', color='r')
axlog[0].set_xlabel(r"Time [s]")
axlog[0].set_ylabel(r"Electron density $n_e$ [$m^{-3}$]")
axlog[0].set_title(r"Height 110km")
axlog[0].legend(["alpha decay", "beta decay", "regular decay"])
# Plotting for height 230km
axlog[1].semilogy(time[3700:], alphaH230km)
axlog[1].semilogy(time[3700:], betaH230km)
axlog[1].semilogy(time[3700:], ionrateH230km[3700:, 0], linestyle='--', color='red')
axlog[1].set_xlabel(r"Time [s]")
axlog[1].set_title(r"Height 230km")
axlog[1].legend(["alpha decay", "beta decay", "regular decay"])
plt.ylim(1E10, 1E12)
# Figure format
figlog.suptitle(r"Comparison of alpha, beta and found decay for electron density (semilog)")
figlog.tight_layout()

# Sinusoidal ionization-rate
sinusoidalH110km = solveODEs(10, time, sinusoidal=True)
sinusoidalH120km = solveODEs(20, time, sinusoidal=True)
sinusoidalH170km = solveODEs(70, time, sinusoidal=True)
sinusoidalH230km = solveODEs(130, time, sinusoidal=True)

# Plotting with sinusoidal ionization-rate
figs, axs = plt.subplots(1, 2, sharey=True)
# Plotting for height 110km
axs[0].semilogy(time[3600:], sinusoidalH110km[3600:])
axs[0].set_xlabel(r"Time [s]")
axs[0].set_ylabel(r"Density [$m^{-3}$]")
axs[0].set_title(r"Height 110km")
axs[0].legend([r"$e^-$", r"$O^{+}$", r"$O_{2}^{+}$", r"$N_{2}^{+}$", r"$NO$", r"$NO^{+}$"])
# Plotting for height 120km
axs[1].semilogy(time[3600:], sinusoidalH120km[3600:])
axs[1].set_xlabel(r"Time [s]")
axs[1].set_title(r"Height 120km")
axs[1].legend([r"$e^-$", r"$O^{+}$", r"$O_{2}^{+}$", r"$N_{2}^{+}$", r"$NO$", r"$NO^{+}$"])
plt.ylim(bottom=1E4)
# Format figure
figs.suptitle(r"Sinusoidal ionization-rate (110km and 120km)")
figs.tight_layout()
#
figs1, axs1 = plt.subplots(1, 2, sharey=True)
# Plotting for height 170km
axs1[0].semilogy(time[3600:], sinusoidalH170km[3600:])
axs1[0].set_xlabel(r"Time [s]")
axs1[0].set_ylabel(r"Density [$m^{-3}$]")
axs1[0].set_title(r"Height 170km")
axs1[0].legend([r"$e^-$", r"$O^{+}$", r"$O_{2}^{+}$", r"$N_{2}^{+}$", r"$NO$", r"$NO^{+}$"])
# Plotting for height 230km
axs1[1].semilogy(time[3600:], sinusoidalH230km[3600:])
axs1[1].set_xlabel(r"Time [s]")
axs1[1].set_title(r"Height 230km")
axs1[1].legend([r"$e^-$", r"$O^{+}$", r"$O_{2}^{+}$", r"$N_{2}^{+}$", r"$NO$", r"$NO^{+}$"])
plt.ylim(bottom=1E6)
# Format figure
figs1.suptitle(r"Sinusoidal ionization-rate (170km and 230km)")
figs1.tight_layout()

if __name__ == "__main__":
    plt.show()

