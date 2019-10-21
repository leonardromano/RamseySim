#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 10:07:24 2019

@author: universe
"""

class DetectorParameters:
    def __init__(self):
        """
        Paramaters are sorted in the following order:
            tdownup, tdownup_err, 
            tdowndown, tdowndown_err, 
            tupup, tupup_err,
            tupdown, tupdown_err, 
            lossUp, lossUp_err, 
            lossDown, lossDown_err
        """
        self.parameters = {
                "Ideal": (0, 0, 1, 0, 1, 0, 0, 0, 0.1462, 0.0024, 0.1480, 0.0024),
                "Abs1SlImp": (0.0437, 0.0018, 0.9577, 0.0027, 0.9561, 0.0016, \
                              0.0442, 0.0018, 0.1720, 0.0014, 0.1741, 0.0014),
                "highAbsSlImp": (0.0451, 0.0019, 0.9511, 0.0024, 0.954,	0.003, \
                                 0.0438, 0.0013, 0.2930, 0.0029, 0.3000, 0.0029),
                "Abs1QuImp": (0.0556, 0.0002, 0.9432, 0.0012, 0.9445, 0.001, \
                              0.0561, 0.0016, 0.2127, 0.0008, 0.2135, 0.0008),
                "NoFilters": (0.4924, 0.0029, 0.5023, 0.0016, 0.505, 0.002, \
                              0.5, 0.0007, 0.0408, 0.0013, 0.0376, 0.0012),
                "Asy1_1": (0.0437, 0.0018, 0.9577, 0.0027, 0.9561, 0.0016, \
                              0.0442, 0.0018, 0.1720, 0.0014, 0.1741, 0.0014),
                "Asy1_2": (0.04419, 0.00026, 0.8863, 0.0012, 0.952, 0.0028, \
                           0.12467, 7E-05, 0.1754, 0.0003, 0.1569, 0.0003),
                "Asy1_3": (0.0441, 0.0003, 0.8224, 0.0026, 0.939, 0.004, \
                           0.1965, 0.0011, 0.1754, 0.0006, 0.1421, 0.0005),
                "Asy1_4": (0.0429, 0.0011, 0.746, 0.007, 0.9553, 0.0008, \
                           0.2631, 0.0029, 0.1741, 0.0006, 0.1318, 0.0023),
                "Asy1_5": (0.0442, 0.0002, 0.724, 0.007, 0.9396, 0.0029, \
                           0.3219, 0.0022, 0.1755, 0.0005, 0.1169, 0.0009),
                "Asy1_6": (0.0425, 0.0006, 0.69, 0.004, 0.924, 0.007, \
                           0.3466, 0.0016, 0.1753, 0.0006, 0.1103, 0.0005),
                "Asy2_1": (0.0538, 0.0002, 0.974, 0.004, 0.925, 0.004, \
                           0.046, 0.0003, 0.2132, 0.0006, 0.1745, 0.0006),
                "Asy2_2": (0.0538, 0.0012, 0.916, 0.008, 0.932, 0.004, \
                           0.1304, 0.0002, 0.2133, 0.0004, 0.1582, 0.0006),
                "Asy2_3": (0.056, 0.0009, 0.817, 0.007, 0.942, 0.004, \
                           0.195, 0.004, 0.2147, 0.0019, 0.1452, 0.0029),
                "Asy2_4": (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                "Asy2_5": (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                "Asy2_6": (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                "Asy3_1": (0.0776, 0.001, 1, 0.009, 0.85, 0.017, \
                           0.0492, 0.0005, 0.3052, 0.0011, 0.1740, 0.0008),
                "Asy3_2": (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                "Asy3_3": (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                "Asy3_4": (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                "Asy3_5": (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                "Asy3_6": (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
                }
        
    def getItemWithKey(self, key):
        "returns the value in the dictionary with some key"
        try:
            return self.parameters.get(key)
        except TypeError:
            print ("wrong key")
            return (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)