#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 13:55:55 2019

@author: universe
"""

import numpy as np

def transmissionMatrix(tup, tdown):
    "Calculate the transmission Matrix"
    return np.array([[tup,0],\
                     [0,tdown]])

def lossMatrix(lossUp, lossDown):
    "Calculate the loss Matrix"
    return np.array([[1-lossUp,0],\
                    [0,1-lossDown]])

def flipPolarization(polVector):
    "Multiply the polarization vector with the first Pauli-Matrix"
    return np.array([polVector[1], polVector[0]])