#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 11:57:36 2019

@author: universe
"""
import Source.Utility.DetectorParameters as DP

tdownup, tdownup_err, tdowndown, tdowndown_err, tupup, tupup_err, \
tupdown, tupdown_err, lossUp, lossUp_err, lossDown, lossDown_err = \
(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
DetectorParameters = DP.DetectorParameters()
key = "Ideal"
detEfficiencyDown0 = 1
delta_detEfficiencyDown = 0
detEfficiencyUp0 = 1
delta_detEfficiencyUp = 0
BackgroundDown0 = 0
delta_BackgroundDown = 0
BackgroundUp0 = 0
delta_BackgroundUp = 0
sigma_gauge = 0.0058

def setTransmissionParameters(key):
    global tdownup, tdownup_err, tdowndown, tdowndown_err, tupup, tupup_err, \
    tupdown, tupdown_err, lossUp, lossUp_err, lossDown, lossDown_err
    tdownup, tdownup_err, tdowndown, tdowndown_err, tupup, tupup_err, \
    tupdown, tupdown_err, lossUp, lossUp_err, lossDown, lossDown_err = \
    DetectorParameters.getItemWithKey(key)

setTransmissionParameters(key)