#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 13:57:32 2019

@author: lenny
"""

import numpy as np
import Source.Utility.Ramsey as Ramsey
import Source.Utility.Math_Utility as mut

from Source.Utility.Constants import *
from Source.Utility.RandomNumberGenerator import random_numbers


def calculateTestData(omega_range, omega_larmor):
    pol_fit = Ramsey.analytical(omega_range, omega_larmor, T, T1, 0, 1)
    min_i, max_i, res_i = mut.find_extrema(pol_fit)
    omega_steep_inc, omega_steep_dec = \
    mut.find_steepest_points(np.array([omega_range[i] for i in min_i]), \
                             np.array([omega_range[i] for i in max_i]))
    
    return pol_fit, min_i, max_i, res_i, omega_steep_inc, omega_steep_dec

def getRandomRange(width, X0, Xend, mu_inc, mu_dec, frac_width, random_numb):
    X_dec                  = random_numbers\
    (mut.gauss_random_density_fwhm, X0, Xend, \
     int(random_numb/2), mu_dec, width, frac_width)
    X_inc                  = random_numbers\
    (mut.gauss_random_density_fwhm, X0, Xend, \
     int(random_numb/2), mu_inc, width, frac_width)
    Xrange                = np.sort(np.append(X_dec, X_inc))
    return Xrange

def makeErrorPositiveAgain(error, length):
    "Dieser Code macht keinen Sinn"
    error_array = np.ones(length)*error
    if error != 0:
        if error < 0:
            error *= (-1)
        if error < 0:
            error_array *= (-1)
    return error, error_array

def includeErrors(polarization, delta_N_up, delta_N_down, N_tot, \
                  add_noise, omega_errorbars, omega_data, omega_noise):
    
    delta_N_up, delta_N_up_array = makeErrorPositiveAgain(delta_N_up, \
                                                          len(polarization))
    delta_N_down, delta_N_down_array = makeErrorPositiveAgain(delta_N_down, \
                                                              len(polarization))
        
    #systematic error of polarization    
    polarization += mut.ramsey_polarization_systematic_error\
    (polarization, N_ges, delta_N_up_array, delta_N_down_array)
        
    # Add noise
    pol_error = mut.ramsey_polarization_counting_error (polarization, N_tot)
    
    if add_noise == True:
        polarization      += np.random.normal(loc = 0, scale = pol_error)
        
    if omega_errorbars == True:
        omega_error        = np.array([random_numbers\
                                       (mut.gauss_random_density, -0.5, \
                                        0.5, 3, 0, omega_noise)[0] for i in \
                                        range(len(omega_data))])
        omega_data        += omega_error
        print("updated omega_data")
    else:
        omega_error = None
    return polarization, pol_error, omega_data, omega_error