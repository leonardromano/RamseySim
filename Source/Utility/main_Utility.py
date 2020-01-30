#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 12:29:51 2019

@author: universe
"""

from Source.Utility.Constants import *
import numpy as np
import Source.Utility.Ramsey as Ramsey
from Source.Utility.Math_Utility import weightedAverage
import Source.Simulation as simulate
import Source.DetectorSystematics_main as systematics

def initOmega(dipole_moment):
    omega_edm = 2*E/hbar * dipole_moment
    omega_larmor_test_data1 = omega_0 + omega_edm
    omega_larmor_test_data2 = omega_0 - omega_edm
    return omega_edm, omega_larmor_test_data1, omega_larmor_test_data2

def initialiseOutputObjects():
    "[omega_larmor, omega_larmor_error, spin_down_sys_count_error, spin_up_sys_count_error]"
    return [np.array([]), np.array([]), np.array([]), np.array([]),]

def registerOutput(OutputObject, Output):
    for i in range(4):
        OutputObject[i] = np.append(OutputObject[i], Output[i])
        
def coincidalPointFraction(xdata, ydata, yerror, func, args):
    dist = np.abs(ydata - func(xdata, *args))
    return sum([1 for i in range(len(dist)) if dist[i] < yerror[i]])/random_numb * 100

def presentResults(domain, parameters, polarization, polError, OutputObj1, \
                   OutputObj2, omega_edm, dipole_moment, afterDetector = False):
    #Calculate average larmor frequency obtained from all runs
    omega_mean1, omega_mean_err1 = weightedAverage(OutputObj1[0], OutputObj1[1]**(-2))
    omega_mean2, omega_mean_err2 = weightedAverage(OutputObj2[0], OutputObj2[1]**(-2))   
    EDM = (omega_mean1-omega_mean2)*hbar/(4*E)
    comment = ''
    #Check if Output has to be adjusted. If more options are added a switch will be more readable
    if afterDetector == True:
        comment = ' after detector systematics'

    print('')
    print('Estimated omega from edm' + comment + ': ' + str(abs(omega_mean1-omega_mean2)/2) \
          + ', True value: ' + str(omega_edm))
    print('')
    print('Estimated EDM: ' + str(EDM) \
          + ', True value: ' + str(dipole_moment))
    print('')
    print('Fraction [percent] of data points coincidal with the fit function ' + \
          '(within error bar): ' + str(coincidalPointFraction(domain, \
           polarization, polError, Ramsey.analytical, parameters)))
    print('')
    print('Larmor frequency 1:     ' + str(omega_mean1) + ' +- ' + str(omega_mean_err1))
    print('Larmor frequency 2:     ' + str(omega_mean2) + ' +- ' + str(omega_mean_err2))
    print('')
    print('EDM:                  ' + str(EDM) \
          + ' +- ' + str(np.sqrt(omega_mean_err1**2 + omega_mean_err2**2) * hbar / (2*E)))
    return EDM


def mainProcedure(OutputObjects, omegaL, dT, dNu, dNd, addNoise, \
                  omega_errorbars, DoPlot, PlotResiduals, fomega, fp, \
                  DoDetectorSystematics= False, flipPolarization = False, \
                  CustomTransmissionProbabilities = False, \
                  Polarizers = None, sigma_gauge = None, key = "0"):
    #Do Ramsey "Simulation" to produce raw data and store it
    param, std, result, omega, pol, omega_fit, pol_error, \
    omega_error, polarizationVector = simulate.and_fit\
    (omegaL, [omega_0, T + dT, T1, 0, 1], random_numb, \
     dNu, dNd, addNoise, omega_errorbars, DoPlot, PlotResiduals, 2e-13, \
     fomega, fp, N_ges)
    registerOutput(OutputObjects[0], [param[0], std[0], dNd, dNu])
    parameters = [param]
    polarisation = [pol]
    polarisation_errors = [pol_error]
    print('Reduced chi-squared: {}'.format(result.redchi))
    if(DoDetectorSystematics == True):
        # Perturb the data with the detector charasteristics and store perturbed data
        param_det, std_det, sys_result, pol_det, pol_det_error, polVector = \
        systematics.main(omega, omega_error, polarizationVector, pol_error, \
                         omega_fit, flipPolarization, dT, param, DoPlot, \
                         PlotResiduals, CustomTransmissionProbabilities, \
                         Polarizers, sigma_gauge, key)
        if(std_det[0] == None):
            "Avoids crashing on unrealistic uncertainties"
            std_det[0] = std[0]
            print("Unrealistic uncertainty occured!")
        registerOutput(OutputObjects[1], [param_det[0], std_det[0], dNd, dNu])
        parameters.append(param_det)
        polarisation.append(pol_det)
        polarisation_errors.append(pol_det_error)
        print('Reduced chi-squared: {}'.format(sys_result.redchi))
    return omega, parameters, polarisation, polarisation_errors