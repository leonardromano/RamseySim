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
                detEfficiencyDown, lossUp, lossDown, tupup_err, tupdown_err, \
                tdowndown_err, tdownup_err, lossUp_err, lossDown_err, sigma_gauge):
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

def PrimaryTransmissionCoefficient(x, y):
    """"Returns the primary transmission Coefficient T1 for intrinsic
    transmission probabilities x and y, referring to the preferred (upup/downdown) and 
    unpreffered (updown/downup) transmission probability"""
    return (1 + (y/0.05) * ((1+x)**T1_par2 /T1_par0/x**T1_par1 - 1))**-1

def PrimaryTransmissionError(x,y, sigma_x = 0, sigma_y = 0):
    "return the statistical error of PrimaryTransmissionCoefficient(x, y)"
    if x==0:
        return 0
    else:
        return PrimaryTransmissionCoefficient(x, y)**2 *np.sqrt\
            (((y/0.05)*(1+x)**T1_par2 /T1_par0/x**T1_par1)**2 * \
             ((np.log(1+x)*T1_err2)**2 + (np.log(x)*T1_err1)**2 + \
              (T1_err0/T1_par0)**2 + ((T1_par2/(1+x) - T1_par1/x)*sigma_x)**2) + \
                 (((1+x)**T1_par2 /T1_par0/x**T1_par1 - 1)*sigma_y/0.05)**2)

def SecondaryTransmissionCoefficient(x, y):
    """"Returns the secondary transmission Coefficient T2 for intrinsic
    transmission probabilities x and y, referring to the preferred (upup/downdown) and 
    unpreffered (updown/downup) transmission probability"""
    return T2_par0*y*x**(-T2_par2) *(1-y)**T2_par1

def SecondaryTransmissionError(x,y, sigma_x = 0, sigma_y = 0):
    "return the statistical error of SecondaryTransmissionCoefficient(x, y)"
    if (y==0 or x==0 or y==1):
        return 0
    else:
        return SecondaryTransmissionCoefficient(x, y)*\
    np.sqrt((T2_err0/T2_par0)**2 + (T2_err1*np.log(1-y))**2+ \
            (T2_err2*np.log(x))**2 + (T2_par2*sigma_x/x)**2 + \
            (sigma_y*(1/y + T2_par1/(1-y)))**2)

def LossCoefficient(x,y):
    "Returns the Loss Coefficient Lambda for intrinsic transmission probabilities x and y"
    return (L_par0 - L_par1*y)*(1+x)**(-L_par2)

def LossError(x,y, sigma_x = 0, sigma_y = 0):
    "Return the statistical error of LossCoefficient(x,y)"
    if (L_par0 - L_par1*y==0):
        return 0
    else:
        return LossCoefficient(x, y)*\
            np.sqrt((np.log(1+x)*L_err2)**2 + (L_par2*sigma_x/x)**2 + \
                    (L_par0 - L_par1*y)**(-2)*((L_err0)**2 + (y*L_err1)**2 +\
                                               (L_par1*sigma_y)**2))