#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 11:38:47 2019

@author: lenny
"""
import numpy as np

from Source.Utility.Constants import *
from scipy import integrate

def poisson_distribution(k, lamb):
    return lamb**k/np.array([np.math.factorial(i) for i in k])*np.exp(-lamb)

def gauss_random_density_fwhm(omega, mu, fwhm, frac_width = random_w_x):
    prob = np.zeros(len(omega))
    for i in range(len(mu)):
        prob += np.sqrt(np.log(2))*gauss_random_density\
        (omega, mu[i], fwhm[i]/2./np.sqrt(2.*np.log(2)))
    prob /= max(prob)
    omega_center = (omega[-1] + omega[0])/2
    modulation = gauss_random_density(omega, omega_center, \
                                      (omega[-1]-omega[0])*frac_width/2./\
                                      np.sqrt(2.*np.log(2.)))
    modulation /= max(modulation)
    prob *= modulation
    return prob

def gauss_random_density(omega, mu, sigma):
    "returns a normalized gaussian function"
    return 1/np.sqrt(2*np.pi*sigma**2) * np.exp(-(((omega - mu)/sigma)**2)/2)

def calculate_effective_frequency(omega, omega_larmor):
    "return the parameters used for the caluclation of the ramsey polarization"
    delta_omega  = omega - omega_larmor
    sigma_pulse = pulse_fwhm / (2*np.sqrt(2*np.log(2)))
        
    #average rabi frequency for the gauss pulses
    t_p        = 2 * sigma_fact * sigma_pulse
    omega_rabi = (gamma_n / t_p) * integrate.quad(lambda t: \
                 gauss_flip_pulse_modulation\
                 (t, sigma_fact*sigma_pulse, pulse_fwhm), 0, t_p)[0]
      
    Omega_eff = np.sqrt(delta_omega**2 + omega_rabi**2)/2
    return delta_omega, t_p, omega_rabi, Omega_eff
    
def ramsey_polarization_counting_error(polarization, N_total):
    "return the statistical error of the ramsey polarization"
    return np.sqrt((1-polarization**2)/N_total)
    
def ramsey_polarization_systematic_error(polarization, N_total, \
                                         delta_N_up = 0, delta_N_down = 0):
    "return the systematic error of the ramsey polarization"
    return ((1-polarization)*delta_N_up - (1+polarization)*delta_N_down)/\
(N_total + delta_N_up + delta_N_down)
    
def gauss_flip_pulse_modulation(t, position, width, gamma = gamma_n, \
                                    flip_angle = np.pi/2):
    "return modulation function of the gaussian flip pulse"
    return (flip_angle/gamma) * gauss_random_density(t, position, \
           width/2./np.sqrt(2*np.log(2)))
    
def gauss_flip_pulse(t, position, width, omega, gamma = gamma_n, \
                         phase = 0, flip_angle = np.pi/2):
    "return gaussian flip pulse"
    return gauss_flip_pulse_modulation(t, position, width, gamma,flip_angle) \
 * np.sin(omega*t + phase)
    
def find_extrema(polarization):
    "Find extrema by comparison to neighbouring points"
    minima_indices = []
    maxima_indices = []
    for i in range(1, len(polarization)-1):
        if polarization[i+1] > polarization[i] and \
        polarization[i-1] > polarization[i]:
            minima_indices.append(i)
        if polarization[i+1] < polarization[i] and \
        polarization[i-1] < polarization[i]:
            maxima_indices.append(i)
    # Find index of resonance frequency
    res_index = [i for i in polarization].index\
    (min([polarization[i] for i in minima_indices]))
    return minima_indices, maxima_indices, res_index
    
def find_steepest_points(omega_minima, omega_maxima):
    "find the steepest points of a function by assuming them in the middle \
    between neighbouring extrema"
    if len(omega_maxima) > len(omega_minima):
        omega_steep_decreasing   = (omega_maxima[:-1] + omega_minima)/2
        omega_steep_increasing   = (omega_maxima[1:] + omega_minima)/2
    elif len(omega_maxima) < len(omega_minima):
        omega_steep_decreasing   = (omega_maxima + omega_minima[1:])/2
        omega_steep_increasing   = (omega_maxima + omega_minima[:-1])/2
    elif len(omega_maxima) == len(omega_minima) and \
    omega_maxima[0] < omega_minima[0]:
        omega_steep_decreasing   = (omega_maxima + omega_minima)/2
        omega_steep_increasing   = (omega_maxima[1:] + omega_minima[:-1])/2
    elif len(omega_maxima) == len(omega_minima) and \
    omega_maxima[0] > omega_minima[0]:
        omega_steep_decreasing   = (omega_maxima[:-1] + omega_minima[1:])/2
        omega_steep_increasing   = (omega_maxima + omega_minima)/2
    return omega_steep_increasing, omega_steep_decreasing

def getPolarization(up, down, offset = 0, scale = 1):
    "returns the polarization from the detector system"
    return (up-down+offset)*scale

def linearFunction(x0, Dx, t):
    "returns a linear function evaluated at t"
    return x0 + Dx*t

def alternatingTrace(matrix):
    "returns the alternating trace of a matrix"
    rowIndex = 0
    trace = 0
    for row in matrix:
        trace+=row[rowIndex]*(-1)**rowIndex
        rowIndex+=1
    return trace

def weightedAverage(ys, dy):
    "calculates the weighted average of y"
    mean = 0
    normalization = 0
    i=0
    for y in ys:
        mean+=y*dy[i]
        normalization+=dy[i]
        i+=1
    return mean/normalization, 1/np.sqrt(normalization)