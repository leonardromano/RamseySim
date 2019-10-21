#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 12:10:32 2019

@author: lenny
"""

import numpy as np
import Source.Utility.Ramsey as Ramsey

from Source.Utility.Constants import *
from lmfit import Parameters, minimize

def ramsey_fringes(omega, polarization, sigma_x = None, \
                       sigma_y = None, p0 = [omega_0, T, T1, 0, 1], \
                       xtol = 1e-12):
    # Initialize fit parameters with bounds
    params = Parameters()
    params.add('omega_larmor', value = p0[0], vary = True)#, min = p0[0]-1, max = p0[0]+1) 
    params.add('T', value = p0[1], vary = True)#, min = p0[1]-10, max = p0[1]+10)	
    params.add('T1', value = p0[2], vary = True)#, min = p0[2]-100, max = p0[2]+100)
    params.add('offset', value = p0[3], vary = False)
    params.add('scale', value = p0[4], vary = False)

    # Minimize chi-squared using the defined residual function 
    # with the initialized parameters and return the fit data
    result = minimize(Ramsey.residual, params, args = \
                      (omega, polarization, sigma_x, sigma_y), \
                      scale_covar = False, method = 'leastsq', \
                      xtol = xtol)
    popt   = np.array([result.params['omega_larmor'].value, \
                       result.params['T'].value, \
                       result.params['T1'].value, \
                       result.params['offset'].value, \
                       result.params['scale'].value])
    stderr = np.array([result.params['omega_larmor'].stderr, \
                       result.params['T'].stderr, \
                       result.params['T1'].stderr, \
                       result.params['offset'].stderr, \
                       result.params['scale'].value])

    return popt, stderr, result