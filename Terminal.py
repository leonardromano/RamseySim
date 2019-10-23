#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 14:07:51 2019

@author: lenny
"""

from Source.main import main
import time

EDM = 1e-26
NumberOfIterations = 100
addNoise = True
omega_errorbars = False
DoPlot = False
PlotResiduals = False
fp = 0.15
fomega = 0.011604
flipPolarization = False
DoDetectorSystematics = True

starttime = time.time()

main(EDM, NumberOfIterations, addNoise, omega_errorbars, DoPlot, PlotResiduals, \
     fp, fomega, flipPolarization, DoDetectorSystematics)
endtime = time.time()
print("Finished in " + str((endtime-starttime)/60.) + " minutes")