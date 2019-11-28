#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 11:43:22 2019

@author: universe
"""
import Source.Utility.DetectorSystematics as run
import Source.Utility.Math_Utility as mut
import Source.Utility.DetectorSystematics_Utility as DSut
import Source.Fitting as fit
import Source.Utility.Ramsey as Ramsey
import Source.Utility.Plotting as plot
import Source.Utility.QAplots as QAplots
from Source.Utility.setDetectorParameters import setDetectorParameters
from Source.Utility.detectorMatrices import flipPolarization as flip
from Source.Utility.SystematicsConstants import *
from Source.Utility.Constants import *

def main(omega_data, omega_error, initialPolVector, initial_pol_error, \
         omega_range, flipPolarization, t, start_params, DoPlot = False, \
         PlotResiduals = False, CustomTransmissionProbabilities =False, \
         Polarizers= [1, 0, 1, 0]):
    "Fits the polarization after beam-propagation through the detector and returns parameters"
    
    detEfficiencyDown = mut.linearFunction(detEfficiencyDown0, delta_detEfficiencyDown, t)
    detEfficiencyUp = mut.linearFunction(detEfficiencyUp0, delta_detEfficiencyUp, t)
    BackgroundDown = mut.linearFunction(BackgroundDown0, delta_BackgroundDown, t)
    BackgroundUp = mut.linearFunction(BackgroundUp0, delta_BackgroundUp, t)
    
    transmission, transmission_err, loss, loss_err = \
    setDetectorParameters(CustomTransmissionProbabilities, *Polarizers)
    
    #spinflip before the polarization enters the detector system
    if flipPolarization:
        initialPolVector = flip(initialPolVector)
    
    NUp, NDown, polarization, Tup, Tdown, Loss, polarizationVector = \
    run.systematics(initialPolVector, *loss, *transmission, N_ges, detEfficiencyDown, \
                    detEfficiencyUp, BackgroundDown, BackgroundUp)
    
    #calculate statistical error of raw polarization
    pol_error = DSut.getPolError(polarizationVector, initialPolVector, \
                                 initial_pol_error, NUp+NDown, N_ges, \
                                 detEfficiencyUp*Tup, \
                                 detEfficiencyDown*Tdown, Loss, \
                                 detEfficiencyUp, detEfficiencyDown, *loss, \
                                 *transmission_err, *loss_err)
    
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
    
    return fit_params, param_std_dev, result, polarization, polarization_fit, \
pol_error, polarizationVector