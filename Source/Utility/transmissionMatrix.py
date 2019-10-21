#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 13:55:55 2019

@author: universe
"""

import numpy as np

def transmissionMatrixElement(pNoAbsorption1, \
                              pNoAbsorption2, \
                              pNoTransmission1, \
                              pNoTransmission2):
    "Calculates the matrix element of the transmission matrix for the detector"
    return pNoAbsorption1*(1-pNoTransmission1)*\
(1+pNoAbsorption2*pNoTransmission2)/(1-pNoAbsorption1*pNoTransmission1*\
pNoAbsorption2*pNoTransmission2)

def transmissionMatrix(r1, r2, r1up, r1down, r2up, r2down):
    "Calculate the transmission Matrix"
    return np.array([[transmissionMatrixElement(r1, r2, r1up, r2down),0],\
                      [0,transmissionMatrixElement(r1, r2, r1down, r2up)]])

def flipPolarization(polVector):
    "Multiply the polarization vector with the first Pauli-Matrix"
    return np.array(polVector[1], polVector[0])