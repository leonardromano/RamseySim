#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 12:38:29 2019

@author: lenny
"""
import numpy as np
import Source.Simulation as simulate
import Source.Utility.Ramsey as Ramsey
import Source.DetectorSystematics_main as systematics

from Source.Utility.Math_Utility import weightedAverage
from Source.Utility.Constants import *

def main(dipole_moment = 1e-26, Niter = 1, addNoise = False, \
         omega_errorbars = False, DoPlot = False, PlotResiduals = False, fp = 0.32, fomega = 0.011604, \
         flipPolarization = False, DoDetectorSystematics = False):
    
    omega_edm = 2*E/hbar * dipole_moment  #3e-7
    omega_larmor_test_data1 = omega_0 + omega_edm
    omega_larmor_test_data2 = omega_0 - omega_edm
    #np.random.seed(37)
    larmor_frequency1      = np.array([])
    larmor_frequency_std_err1 = np.array([])
    spin_down_sys_count_error1 = np.array([])
    spin_up_sys_count_error1 = np.array([])
    
    larmor_frequency2      = np.array([])
    larmor_frequency_std_err2 = np.array([])
    spin_down_sys_count_error2 = np.array([])
    spin_up_sys_count_error2 = np.array([])
    #random_w_x_array = np.array([])
    #random_w_y_array = np.array([])
    
    if(DoDetectorSystematics == True):
        larmor_frequency_det1        = np.array([])
        larmor_frequency_det_std_err1 = np.array([])
        spin_down_sys_count_error1 = np.array([])
        spin_up_sys_count_error1 = np.array([])
    
        larmor_frequency_det2       = np.array([])
        larmor_frequency_det_std_err2 = np.array([])
        spin_down_sys_count_error2 = np.array([])
        spin_up_sys_count_error2 = np.array([])
        BadFitCounter = 0
    
    

    for i in range(Niter):
        print("Performing Ramsey simulation {}/{}".format(i+1,Niter))
        random_w_x = fp#0.32*((i+1)*0.01)
        random_w_y = fomega#0.0041+6.7e-5*i#0.0125/3*(i+1)#0.0125#0.125*0.1*(i+1)
        
        dNd = 0
        dNu = 0
        dT = 0.01
        
        # Do the procedure for the E-up data
        param1, std1, result, omega, pol, omega_fit, pol_fit, pol_error, \
        omega_error, mini, maxi, res_i, polarizationVector = simulate.and_fit\
        (omega_larmor_test_data1, [omega_0, T + dT, T1, 0, 1], random_numb, \
         dNu, dNd, addNoise, omega_errorbars, DoPlot, PlotResiduals, 2e-13, \
         random_w_y, random_w_x, N_ges)
        
        larmor_frequency1         = np.append(larmor_frequency1, param1[0])
        larmor_frequency_std_err1 = np.append(larmor_frequency_std_err1, std1[0])
        spin_down_sys_count_error1 = np.append(spin_down_sys_count_error1, dNd)
        spin_up_sys_count_error1 = np.append(spin_up_sys_count_error1, dNu)
        print('Reduced chi-squared: {}'.format(result.redchi))
        
        if(DoDetectorSystematics == True):
            fit_params1, std_dev1, sys_result, polarization, polarization_fit, \
            polar_error, polVector = \
            systematics.main(omega, omega_error, polarizationVector, pol_error, \
                             omega_fit, flipPolarization, i*dT, param1, DoPlot, \
                             PlotResiduals)
            if(std_dev1[0] == None):
                "Avoids crashing on unrealistic uncertainties"
                std_dev1[0] = std1[0]
                print("Unrealistic uncertainty occured!")
                BadFitCounter+=1
            larmor_frequency_det1        = np.append(larmor_frequency_det1, fit_params1[0])
            larmor_frequency_det_std_err1 = np.append(larmor_frequency_det_std_err1, std_dev1[0])
            print('Reduced chi-squared: {}'.format(sys_result.redchi))
        
        # Do the procedure for the E-down data
        param2, std2, result, omega, pol, omega_fit, pol_fit, pol_error, \
        omega_error, mini, maxi, res_i, polarizationVector = simulate.and_fit\
        (omega_larmor_test_data2, [omega_0, T + dT, T1, 0, 1], random_numb, \
         dNu, dNd, addNoise, omega_errorbars, DoPlot, PlotResiduals, 2e-13, \
         random_w_y, random_w_x, N_ges)
        
        larmor_frequency2         = np.append(larmor_frequency2, param2[0])
        larmor_frequency_std_err2 = np.append(larmor_frequency_std_err2, std2[0])
        spin_down_sys_count_error2 = np.append(spin_down_sys_count_error2, dNd)
        spin_up_sys_count_error2 = np.append(spin_up_sys_count_error2, dNu)
        print('Reduced chi-squared: {}'.format(result.redchi))
        
        if(DoDetectorSystematics == True):
            fit_params2, std_dev2, sys_result, polarization, polarization_fit, \
            polar_error, polVector = \
            systematics.main(omega, omega_error, polarizationVector, pol_error, \
                             omega_fit, flipPolarization, i*dT, param2, DoPlot, \
                             PlotResiduals)
            if(std_dev2[0] == None):
                "Avoids crashing on unrealistic uncertainties"
                std_dev2[0] = std2[0]
                print("Unrealistic uncertainty occured!")
                BadFitCounter+=1
            larmor_frequency_det2       = np.append(larmor_frequency_det2, fit_params2[0])
            larmor_frequency_det_std_err2 = np.append(larmor_frequency_det_std_err2, std_dev2[0])
            print('Reduced chi-squared: {}'.format(sys_result.redchi))

    dist1 = np.abs(pol - Ramsey.analytical(omega, *param1))
    inside_count = sum([1 for i in range(len(dist1)) if \
                        dist1[i] < pol_error[i]])/random_numb * 100
    
    weights1 = 1./larmor_frequency_std_err1**2
    
    omega_mean1, omega_mean_err1 = weightedAverage(larmor_frequency1, weights1)
    
    weights2 = 1./larmor_frequency_std_err1**2
    
    omega_mean2, omega_mean_err2 = weightedAverage(larmor_frequency2, weights2)
    
    EDM = (omega_mean1-omega_mean2)*hbar/(4*E)

    print('')
    print('Estimated omega from edm: ' + str(abs(omega_mean1-omega_mean2)/2) \
          + ', True value: ' + str(omega_edm))
    print('')
    print('Estimated EDM: ' + str(EDM) \
          + ', True value: ' + str(dipole_moment))
    print('')
    print('Percentage of data points coincidal with the fit function: \
          (within error bar): ' + str(inside_count))
    print('')
    print('Larmor frequency 1:     ' + str(omega_mean1) + ' +- ' + str(omega_mean_err1))
    print('Larmor frequency 2:     ' + str(omega_mean2) + ' +- ' + str(omega_mean_err2))
    #print('Free precession time: ' + str(param[1]) + ' +- ' + str(std[1]))
    #print('T1 time:              ' + str(param[2]) + ' +- ' + str(std[2]))
    print('')
    print('EDM:                  ' + str(EDM) \
          + ' +- ' + str((omega_mean_err1 + omega_mean_err2) * hbar / (2*E)))
    
    if(DoDetectorSystematics == True):
        weights_det1 = 1./larmor_frequency_det_std_err1**2
        omega_mean_det1, omega_mean_det_err1 = weightedAverage(larmor_frequency_det1, weights_det1)
    
        weights_det2 = 1./larmor_frequency_det_std_err2**2
        omega_mean_det2, omega_mean_det_err2 = weightedAverage(larmor_frequency_det2, weights_det2)

        dist_det = np.abs(polarization - Ramsey.analytical(omega, *fit_params1))
        inside_count_det = sum([1 for i in range(len(dist_det)) if \
                        dist_det[i] < polar_error[i]])/random_numb * 100
    
        EDM_det = (omega_mean_det1-omega_mean_det2)*hbar/(4*E)

        print('')
        print('Estimated omega from edm after detector systematics: ' + \
              str(abs(omega_mean_det1-omega_mean_det2)) + ', True value: ' + str(omega_edm))
        print('')
        print('Estimated EDM: ' + str(EDM_det) \
          + ', True value: ' + str(dipole_moment))
        print('')
        print('Percentage of data points coincidal with the fit function: \
              (within error bar): ' + str(inside_count_det))
        print('')
        print('Larmor frequency 1:     ' + str(omega_mean_det1) + ' +- ' + str(omega_mean_det_err1))
        print('Larmor frequency 2:     ' + str(omega_mean_det2) + ' +- ' + str(omega_mean_det_err2))
        #print('Free precession time: ' + str(param[1]) + ' +- ' + str(std[1]))
        #print('T1 time:              ' + str(param[2]) + ' +- ' + str(std[2]))
        print('')
        print('EDM:                  ' + str(EDM_det) \
              + ' +- ' + str((omega_mean_det_err1 + omega_mean_det_err2) * hbar / (4*E)))
    
        print('Systematic Shift due to detector:       ' + str(EDM_det - EDM))
        
        print('Fraction of bad fits: ' + str(100*BadFitCounter/2/Niter) +' %')

