#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 14:07:51 2019

@author: lenny
"""

from Source.main import main
import time

EDM = 3e-26
NumberOfIterations = 1
addNoise = True
omega_errorbars = False
DoPlot = False
PlotResiduals = False
fp = 0.15
fomega = 0.011604
flipPolarization = False
DoDetectorSystematics = True
CustomTransProb = False
Polarizers_system1 = [0.75, 0.05, 0.9, 0.2] #[upup, updown, downdown, downup]
Polarizers_system2 = [0.90, 0.05, 0.90, 0.0] #[upup, updown, downdown, downup]
keys = ["Asy1_20_inv", "Asy1_2_inv"]
sigma_gauge = [48.6, 23.05]

#Time measurement start
starttime = time.time()

main(EDM, NumberOfIterations, addNoise, omega_errorbars, DoPlot, PlotResiduals, \
     fp, fomega, flipPolarization, DoDetectorSystematics, CustomTransProb, \
     Polarizers_system1, Polarizers_system2, *sigma_gauge, keys)

#Time measurement stop
endtime = time.time()
print("Finished in " + str((endtime-starttime)/60.) + " minutes")