#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 10:41:55 2019

@author: lenny
"""

import Source.Utility.DetectorSystematics_Utility as DSut
import Source.Utility.Math_Utility as mut
import Source.Utility.detectorMatrices as DM
import numpy as np

def systematics(initialPolVector, lossUp, lossDown, tupup, tupdown, tdowndown, \
                tdownup, N_total, detEfficiencyUp, detEfficiencyDown, \
                BackgroundUp, BackgroundDown):
    "return yields, polarization, purities and fractions of correctly \
    identified spins"
    #initialise transmission- and loss-matrix
    Tup = DM.transmissionMatrix(tupup, tupdown)
    Tdown = DM.transmissionMatrix(tdownup, tdowndown)
    Loss = DM.lossMatrix(lossUp,lossDown)
    
    #apply matrices to polarisation vector 
    temp = np.dot(Loss, initialPolVector)
    polUp = np.dot(Tup, temp)
    polDown = np.dot(Tdown, temp)
    
    #Extract observables
    NUp = DSut.totalDetectorCount(polUp, N_total, BackgroundUp, detEfficiencyUp)
    NDown = DSut.totalDetectorCount(polDown, N_total, BackgroundDown, detEfficiencyDown)
    NtotDetector = NUp + NDown
    polarizationVector = np.array([NUp/NtotDetector, NDown/NtotDetector])
    Polarization = mut.getPolarization(polarizationVector[0], polarizationVector[1])
    return NUp, NDown, Polarization, Tup, Tdown, Loss, polarizationVector
    