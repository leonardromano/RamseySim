#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 12:27:38 2019

@author: lenny
"""

import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

def fit_and_data(omega_data, polarization_data, \
                      polarization_error, omega_error, omega_fit, \
                      polarization_fit):
    #Plot fitted function
    plt.figure('Ramsey fringes')
    plt.clf()
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Polarization')
    plt.plot(omega_fit/(2*np.pi), polarization_fit, \
             label = 'Fitted Ramsey pattern')

    # Plot calculated Ramsey pattern
    plt.errorbar(omega_data/(2*np.pi), polarization_data, \
                 polarization_error, omega_error, fmt = 'o', \
                 markersize = 3, capsize = 4, color = 'red')

    # Plot legend
    plt.legend()

    #----
    fig, ax = plt.subplots()
    ax.plot(omega_fit/(2*np.pi), polarization_fit, \
            label = 'Fitted Ramsey pattern')
    ax.errorbar(omega_data/(2*np.pi), polarization_data, \
                polarization_error, omega_error, fmt = 'o', \
                markersize = 3, capsize = 4, color = 'red')
    #ax.legend()
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Polarization')
    axins = inset_axes(ax, '25%','25%',loc = 1)
    axins.plot(omega_fit/(2*np.pi), polarization_fit)
    axins.errorbar(omega_data/(2*np.pi), polarization_data, \
                   polarization_error, omega_error, fmt = 'o', \
                   markersize = 3, capsize = 4, color = 'red')
    axins.set_xlim(36.8311,36.8331)
    axins.set_ylim(-0.75,1.03)
    mark_inset(ax, axins, loc1=2, loc2=4, fc='green', ec='black')
    axins.get_xaxis().get_major_formatter().set_useOffset(False)
    axins.get_yaxis().get_major_formatter().set_useOffset(False)
    plt.yticks(visible=False) 
    plt.xticks(visible=False)
    #----

    #Plot the frequency distribution of the test data
    plt.figure('Random numbers')
    plt.clf()
    plt.hist(omega_data/(2*np.pi), 200, facecolor = 'red', alpha = 0.75)
    
def randomGauss(Xrange, distribution):
    fig, ax = plt.subplots()
    ax.plot(Xrange/np.pi/2, distribution)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Probability')
    axins = inset_axes(ax, '25%','25%',loc = 1)
    axins.plot(Xrange/np.pi/2,distribution)
    axins.set_xlim(36.8311,36.8331)
    axins.set_ylim(-0.024,1.03)
    mark_inset(ax, axins, loc1=2, loc2=4, fc='green', ec='black')
    axins.get_xaxis().get_major_formatter().set_useOffset(False)
    axins.get_yaxis().get_major_formatter().set_useOffset(False)
    plt.yticks(visible=False) 
    plt.xticks(visible=False)