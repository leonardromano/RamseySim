#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 12:17:45 2019

@author: lenny
"""
import numpy as np
import Source.Utility.Ramsey as Ramsey
import Source.Utility.Math_Utility as mut
import Source.Utility.QAplots as QAplots
import Source.Fitting as fit
import Source.Utility.Plotting as plot
import Source.Utility.Simulation_Utility as simut

from Source.Utility.Constants import *

def and_fit(omega_larmor, start_params, \
            random_number = random_numb, delta_N_up = 0, \
            delta_N_down = 0, add_noise = True, \
            omega_errorbars = True, DoPlot = False, PlotResiduals = False, \
            ome_err = 0.1, rnd_w_y = random_w_y, rnd_w_x = random_w_x, \
            N_tot = N_ges):
        
    omega_range                 = np.arange(omega_start, omega_end, 0.000001)
    # Calculate ramsey test data
    polarization_fit, minima_i, maxima_i, res_i, \
    omega_steep_increasing, omega_steep_decreasing  = \
    simut.calculateTestData(omega_range, omega_larmor)

    # Calculate new random test data on steepest points of the fit
    width                     = np.abs(omega_steep_decreasing - \
                                       omega_steep_increasing) * rnd_w_y
    omega_data                = simut.getRandomRange\
    (width, omega_start, omega_end, omega_steep_increasing, \
     omega_steep_decreasing, random_w_x, random_number)
    
    delta_omega_data, t_p_data, omega_rabi_data, Omega_eff_data = \
    mut.calculate_effective_frequency(omega_data, omega_larmor)
    
    polarizationVector        = Ramsey.polarizationVector\
    (delta_omega_data, t_p_data, omega_rabi_data, Omega_eff_data, T = T, \
     T1 = T1, alpha_0 = alpha_0)
    
    polarization              = mut.getPolarization(polarizationVector[0], \
                                                    polarizationVector[1], 0, 1)
        
    #---- Plot distribution of frequencies
    """
    plot.randomGauss(omega_range, mut.gauss_random_density_fwhm\
                     (omega_range, omega_steep_decreasing, width, \
                      frac_width = random_w_x) + mut.gauss_random_density_fwhm\
                      (omega_range, omega_steep_increasing, width, \
                       frac_width = rnd_w_x)) 
    """
    # Systematic neutron count error polarization shift
    
    polarization, pol_error, \
    omega_data, omega_error = simut.includeErrors(polarization, delta_N_up, \
                                            delta_N_down, N_tot, add_noise, \
                                            omega_errorbars, omega_data, \
                                            ome_err)

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

    return fit_params, param_std_dev, result, omega_data, polarization, \
omega_range, polarization_fit, pol_error, omega_error, minima_i, maxima_i, \
res_i, polarizationVector