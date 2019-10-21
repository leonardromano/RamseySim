#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 15:11:48 2019

@author: lenny
"""

import numpy as np
import Source.Utility.Math_Utility as mut
from Source.Utility.SystematicsConstants import *

def totalProbability(polVector):
    "Calculate the total probability of a Partcle to be in a certain state"
    return polVector[0] + polVector[1]

def totalCount(polVector, N_total):
    "Calculate the total Count of particles to be in a certain state weighted\
    with the probability to be in said state"
    return totalProbability(polVector)*N_total

def totalDetectorCount(polVector, N_total, Background, detEfficiency):
    "Calculate the total count in the detector"
    return detEfficiency*(totalCount(polVector, N_total)+Background)

def purity(initialPolVector, polState, polVector, N_total, \
           Background, detEfficiency):
    "Calculate the purity or Fraction of the detector signal"
    if polState == "up":
        return initialPolVector[0]*N_total/\
    totalDetectorCount(polVector, N_total, Background, detEfficiency)
    elif polState == "down":
        return initialPolVector[1]*N_total/\
    totalDetectorCount(polVector, N_total, Background, detEfficiency)
    
def getLeftAndRightPolarizationVectors(initialPolVector, detectorSwitch, \
                                       detectorSystemState):
    "Returns the polarization vector on each side of the detector \
    for a certain state of the detector system"
    detectorSwitch.setLeftRight(initialPolVector)
    return detectorSwitch.getItemWithKey(detectorSystemState)

def getNumberOfReflections(N_total, reflecProb):
    "returns the number of reflections"
    n = 0
    while(0.5 < N_total*reflecProb**n):
        n+=1
    print("Maximum number of reflections: " + n)
    return n

def getPolError(polAfter, initialPolVector, initial_pol_error, \
                N_total_detector, N_ges, redTransMatrixUp, \
                redTransMatrixDown, LossMatrix, detEfficiencyUp, \
                detEfficiencyDown):
    "returns the polarization error of the detector"
    return sigma_gauge*N_ges/N_total_detector*\
np.sqrt((initial_pol_error/2*dpAfterdpInitial(np.dot(redTransMatrixUp, \
                                                      LossMatrix), \
    np.dot(redTransMatrixDown, LossMatrix), polAfter[0]))**2 + \
    (lossUp*initialPolVector[0])**2*\
    ((tupup_err*detEfficiencyUp*polAfter[1])**2+\
     (tdownup_err*detEfficiencyDown*polAfter[0])**2) +\
     (lossDown*initialPolVector[1])**2*\
     ((tupdown_err*detEfficiencyUp*polAfter[1])**2+\
      (tdowndown_err*detEfficiencyDown*polAfter[0])**2)+\
      (lossUp_err*\
       dpAfterdloss(initialPolVector[0],polAfter, redTransMatrixUp[0,0], \
                      redTransMatrixDown[1,1]))**2 +\
       (lossDown_err*\
       dpAfterdloss(initialPolVector[1], polAfter, redTransMatrixUp[1,1], \
                      redTransMatrixDown[0,0]))**2)
    
def dpAfterdpInitial(TupDelta, TdownDelta, pup):
    "returns the normalized derivative of the polarization with respect to the initial polarization"
    return mut.alternatingTrace(TupDelta)-mut.alternatingTrace(TupDelta+TdownDelta)*pup

def dpAfterdloss(p0, pAfter, redTransMatrixElementUp, redTransMatrixElementDown):
    "returns the normalized derivative of the polarization with respect to the loss Up"
    return p0*(redTransMatrixElementUp*pAfter[1]-redTransMatrixElementDown*pAfter[0])