#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 11:43:22 2019

@author: universe
"""
import Source.Utility.DetectorSystematics as run
import Source.Utility.Math_Utility as mut
import Source.Utility.DetectorSystematics_Utility as DSut
import Source.Utility.Fitting as fit
import Source.Utility.Ramsey as Ramsey
import Source.Utility.Plotting as plot
import Source.Utility.QAplots as QAplots
import numpy as np
from Source.Utility.setDetectorParameters import setDetectorParameters
from Source.Utility.detectorMatrices import flipPolarization as flip
from Source.Utility.SystematicsConstants import *
from Source.Utility.Constants import *

def main(omega_data, omega_error, initialPolVector, initial_pol_error, \
         omega_range, flipPolarization, t, start_params, DoPlot = False, \
         PlotResiduals = False, CustomTransmissionProbabilities =False, \
         Polarizers= [1, 0, 1, 0], sigma_gauge = 2.0, key = "0"):
    "Fits the polarization after beam-propagation through the detector and returns parameters"
    
    #initialise all parameters
    detEfficiency = mut.linearFunction(np.array([detEfficiencyUp0, \
                                                 detEfficiencyDown0]), \
    np.array([delta_detEfficiencyUp, delta_detEfficiencyDown]), t)
    background = mut.linearFunction(np.array([BackgroundUp0, BackgroundDown0]), \
    np.array([delta_BackgroundUp, delta_BackgroundDown]), t)
    transmission, transmission_err, loss, loss_err = \
    setDetectorParameters(CustomTransmissionProbabilities, *Polarizers, key)
    
    #spinflip before the polarization enters the detector system
    if flipPolarization:
        initialPolVector = flip(initialPolVector)
    
    NUp, NDown, polarization, Tup, Tdown, Loss, polarizationVector = \
    run.systematics(initialPolVector, *loss, *transmission, N_ges, \
                    *detEfficiency, *background)
    
    #calculate statistical error of raw polarization
    pol_error = DSut.getPolError(polarizationVector, initialPolVector, \
                                 initial_pol_error, NUp+NDown, N_ges, \
                                 detEfficiency[0]*Tup, \
                                 detEfficiency[1]*Tdown, Loss, \
                                 *detEfficiency, *loss, \
                                 *transmission_err, *loss_err, sigma_gauge)
    
    # Fit ramsey pattern to calculated test data
    fit_params, param_std_dev, result  = fit.ramsey_fringes\
    (omega_data, polarization, omega_error, pol_error, p0 = start_params)
    
    # Calculate fitted ramsey pattern
    polarization_fit      = Ramsey.analytical(omega_range, *fit_params)
    
    if PlotResiduals == True:
        QAplots.plotResiduals(omega_data, polarization, pol_error, \
                              Ramsey.analytical(omega_data, *fit_params))
    
    if DoPlot == True:
        plot.fit_and_data(omega_data, polarization, pol_error, omega_error, \
             omega_range, polarization_fit)
    
    return fit_params, param_std_dev, result, polarization, pol_error, \
polarizationVector
