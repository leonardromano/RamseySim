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
DoPlot = True
PlotResiduals = True
fp = 0.15
fomega = 0.011604
flipPolarization = False
DoDetectorSystematics = True
CustomTransProb = True
Polarizers_system1 = [0.9, 0.05, 0.9, 0.05] #[upup, updown, downdown, downup]
Polarizers_system2 = [0.85, 0.05, 0.9, 0.05]
sigma_gauge = [0.391, 0.42]

starttime = time.time()

main(EDM, NumberOfIterations, addNoise, omega_errorbars, DoPlot, PlotResiduals, \
     fp, fomega, flipPolarization, DoDetectorSystematics, CustomTransProb, \
     Polarizers_system1, Polarizers_system2, *sigma_gauge)
endtime = time.time()
print("Finished in " + str((endtime-starttime)/60.) + " minutes")