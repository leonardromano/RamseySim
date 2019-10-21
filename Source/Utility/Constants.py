#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 11:46:37 2019

@author: lenny
"""
import numpy as np

#______________________________________________________________________________

# Neutron parameters and magnetic environment
#B1_amp      = 0.0018#np.sqrt(np.log(2)/np.pi) / (gamma_n*2) #0.0028/5           
# Pulse field amplitude in muT #0.0085733882#
hbar        = 6.58211951440 * 10**(-16)             # eVs
E           = 2.*10**4                               # V/cm
pulse_fwhm  = 2.                                     # Gauß pulse fwhm in s
sigma_fact  = 5.                                     # (2 * sigma_frac * sigma) pulse lenght for calculation of average rabi frequency
B_0         = 1.263                                 # Holding filed in muT
gamma_n     = 2. * np.pi * 29.1646943                # Hz/muT, (omega = gamma * B)
omega_0     = gamma_n*B_0                           # Neutron resonance frequency for the given B_0
T           = 250.                                   # Neutron free precession duration
T1          = 4400.                                  # Longitudinal spin relaxation time
alpha_0     = 0.90                                  # Visibility at t = 0 after filling
alpha       = alpha_0*np.exp(-T/T1) #0.85           # Visibility
N_ges       = 57150.                                 # N_up + N_down of the respective Ramsey cycle
omega_start = 231.17                                # Start frequency for the simulation
omega_end   = 231.71                                # End Frequency for the simulation
pulse_field = 'linear'                              # Shape of the pulse field. Options are 'linear' and 'circular'
random_w_x  = 0.15                                 # width fraction for the data points on the omega axis
random_w_y  = 0.125                               # width fraction for the data points on the polarization axis
random_numb = 21600.                                # Number of data points per Ramsey cycle

#__________________________________________________________________________