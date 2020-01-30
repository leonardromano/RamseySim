#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 14:24:44 2019

@author: universe
"""

import numpy as np
import Source.Utility.DetectorParameters as DP
import Source.Utility.DetectorSystematics_Utility as dSU

def setDetectorParameters(CustomTransmissionProbabilities, pupup, pupdown, \
                          pdowndown, pdownup, key = "0", pupup_err = 0, pupdown_err = 0, \
                          pdowndown_err = 0, pdownup_err = 0):
    "Decides the mode how to determine the detector Parameters and sets them"
    if (CustomTransmissionProbabilities == True):
        print(pupup, pdownup)
        transmission = \
        np.array([dSU.PrimaryTransmissionCoefficient(pupup, pdownup), \
                  dSU.SecondaryTransmissionCoefficient(pdowndown, pupdown), \
                  dSU.PrimaryTransmissionCoefficient(pdowndown, pupdown), \
                  dSU.SecondaryTransmissionCoefficient(pupup, pdownup)])
        print(transmission)
        transmission_err = \
        np.array([dSU.PrimaryTransmissionError(pupup, pdownup, pupup_err, pdownup_err), \
                  dSU.SecondaryTransmissionError(pdowndown, pupdown, pdowndown_err, pupdown_err), \
                  dSU.PrimaryTransmissionError(pdowndown, pupdown, pdowndown_err, pupdown_err), \
                  dSU.SecondaryTransmissionError(pupup, pdownup, pupup_err, pdownup_err)])
        loss = np.array([dSU.LossCoefficient(pupup, pdownup), \
                         dSU.LossCoefficient(pdowndown, pupdown)])
        print(loss)
        loss_err = \
        np.array([dSU.LossError(pupup, pdownup, pupup_err, pdownup_err), \
                  dSU.LossError(pdowndown, pupdown, pdowndown_err, pupdown_err)])
        return transmission, transmission_err, loss, loss_err
    else:
        DetectorParameters = DP.DetectorParameters()
        tdownup, tdownup_err, tdowndown, tdowndown_err, tupup, tupup_err, \
        tupdown, tupdown_err, lossUp, lossUp_err, lossDown, lossDown_err = \
        DetectorParameters.getItemWithKey(key)
        transmission = np.array([tupup, tupdown, tdowndown, tdownup])
        transmission_err = np.array([tupup_err, tupdown_err, \
                                     tdowndown_err, tdownup_err])
        loss = np.array([lossUp, lossDown])
        loss_err = np.array([lossUp_err, lossDown_err])
        return transmission, transmission_err, loss, loss_err