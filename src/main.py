"""
@author: Daniel Elisabethsønn Antonsen, UiT Institute of Physics and Technology

Main module assignment
"""
# Importing libraries and modules
import os
import numpy as np
# Opening files
def openData(filename):
    return np.loadtxt(os.path.join("resources", filename))

listOfData = ["MSIS.dat", ]


