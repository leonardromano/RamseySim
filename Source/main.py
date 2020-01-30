#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 12:38:29 2019

@author: lenny
"""


import Source.Utility.main_Utility as ut

def main(dipole_moment = 1e-26, Niter = 1, addNoise = False, \
         omega_errorbars = False, DoPlot = False, PlotResiduals = False, \
         fp = 0.32, fomega = 0.011604, flipPolarization = False, \
         DoDetectorSystematics = False, CustomTransmissionProbabilities = False, \
         Polarizers_system1 = [1, 0, 1, 0], Polarizers_system2 = [1, 0, 1, 0], \
         sigma_gauge1 = 2.0, sigma_gauge2 = 2.0, keys = ["0", "0"]):
    #initialise initial values
    omegas = ut.initOmega(dipole_moment)
    #initialise Output Objects
    OutputObjects1 = [ut.initialiseOutputObjects()]
    OutputObjects2 = [ut.initialiseOutputObjects()]
    
    if(DoDetectorSystematics == True):
        OutputObjects1.append(ut.initialiseOutputObjects())
        OutputObjects2.append(ut.initialiseOutputObjects())

    for i in range(Niter):
        print("Performing Ramsey simulation {}/{}".format(i+1,Niter))
        #initialise time Coordinate and Systematic shifts
        dNd = 0
        dNu = 0
        dT = 0
        # Do the procedure for the E-up data
        omega, params1, pols1, pol_errors1 = \
        ut.mainProcedure(OutputObjects1, omegas[1], i*dT, dNu, dNd, addNoise, \
                         omega_errorbars, DoPlot, PlotResiduals, fomega, fp, \
                         DoDetectorSystematics, flipPolarization, \
                         CustomTransmissionProbabilities, Polarizers_system1, \
                         sigma_gauge1, keys[0])
        # Do the procedure for the E-down data
        omega, params2, pols2, pol_errors2 = \
        ut.mainProcedure(OutputObjects2, omegas[2], i*dT, dNu, dNd, addNoise, \
                         omega_errorbars, DoPlot, PlotResiduals, fomega, fp, \
                         DoDetectorSystematics, flipPolarization, \
                         CustomTransmissionProbabilities, Polarizers_system2, \
                         sigma_gauge2, keys[1])
    # Print the results
    EDM = ut.presentResults(omega, params1[0], pols1[0], pol_errors1[0], \
                            OutputObjects1[0], OutputObjects2[0], \
                            omegas[0], dipole_moment)
    
    if(DoDetectorSystematics == True):
        EDM_det = ut.presentResults(omega, params1[1], pols1[1], \
                                    pol_errors1[1], OutputObjects1[1], \
                                    OutputObjects2[1], omegas[0], \
                                    dipole_moment, True)
    
        print('Systematic Shift due to detector:       ' + str(EDM_det - EDM))

