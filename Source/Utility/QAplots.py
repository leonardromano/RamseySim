#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 13:16:23 2019

@author: universe
"""
import numpy as np
import matplotlib.pyplot as plt

def plotResiduals(x, y, y_error, modelFunction):
    "plot the residual of the fit"
    #Plot fitted function
    plt.figure()
    plt.clf()
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Residual $(P_{sim} - P_{Ramsey})/\sigma$')
    plt.errorbar(x/(2*np.pi),(y-modelFunction)/y_error ,y_error, marker = "o", \
                 ecolor = "black", mfc = "blue")
    plt.show()
    plt.close()