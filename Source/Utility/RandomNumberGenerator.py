#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 13:12:18 2019

@author: lenny
"""

import numpy as np

# _________________________________________________________________________________

def integrate_probability_density (density_function , x_start , \
                                   x_end , length , *params):
    x_list = np.arange(x_start , x_end , (x_end - x_start)/length , \
                       dtype=np.float64)
    F_list = np.cumsum(density_function(x_list , *params))
    F_list /= F_list [-1]
    return x_list , F_list

def find_interpolation_points (y, x_list , F_list):
    lower_point_indices = np.array ([np.where(F_list < i)[0][ -1] for i in y])
    upper_point_indices = lower_point_indices + 1
    x_low = np.array ([ x_list[i] for i in lower_point_indices ])
    x_high = np.array ([ x_list[i] for i in upper_point_indices ])
    y_low = np.array ([ F_list[i] for i in lower_point_indices ])
    y_high = np.array ([ F_list[i] for i in upper_point_indices ])
    return x_low , x_high , y_low , y_high

def interpolate_x(y, x_low , x_high , y_low , y_high):
    return (y - y_high)*( x_high - x_low)/( y_high - y_low) + x_high

def random_numbers(density_function , x_start , x_end , length , \
                   *params):
    random_uniform = np.random.uniform (0, 1, length)
    x_list, F_list = integrate_probability_density\
    (density_function, x_start, x_end, length, *params)
    x_low , x_high , y_low , y_high = \
    find_interpolation_points(random_uniform, x_list, F_list)
    return interpolate_x(random_uniform , x_low , x_high , y_low , y_high)

#if __name__ == '__main__ ':
#rnd = RandomNumberGenerator ()
#rand = rnd.random_numbers(<dens_function >, 0, 8, 20000 , 1)

# Integrate a function h(x,* params) between x1 and x2:
#def h(x,a,b,c):
#return <function evaluation code >
#params = (1,2,3) #(a,b,c) values
#integral = integrate.quad(lambda x : h(x, *params), x1 , x2)