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
                "Asy1_1": (0.0444, 0.0007, 0.9554, 0.0007, 0.9555, 0.0004, \
                           0.0447, 0.0005, 0.1740, 0.0005, 0.17434, 0.00029),
                "Asy1_2": (0.04420,	0.00025, 0.8802, 0.0029, 0.9549, 0.0015, \
                           0.12464, 0.00012, 0.17537, 0.00029, 0.15728, 0.00027),
                "Asy1_2_inv": (0.12464, 0.00012, 0.9549, 0.0015, 0.8802, 0.0029,\
                           0.04420,	0.00025, 0.15728, 0.00027, 0.17537, 0.00029),
                "Asy1_3": (0.0441, 0.0003, 0.8224, 0.0026, 0.939, 0.004, \
                           0.1965, 0.0011, 0.1754, 0.0006, 0.1421, 0.0005),
                "Asy1_3_inv": (0.1965, 0.0011, 0.939, 0.004, 0.8224, 0.0026,\
                           0.0441, 0.0003, 0.1421, 0.0005, 0.1754, 0.0006),
                "Asy1_4": (0.0429, 0.0011, 0.74698, 0.00022, 0.9553, 0.0008, \
                           0.2532, 0.0008, 0.1741, 0.0006, 0.12944, 0.00028),
                "Asy1_5": (0.04436, 0.00015, 0.697, 0.008, 0.953, 0.004, \
                           0.312, 0.004, 0.1750, 0.0004, 0.1168, 0.0005),
                "Asy1_6": (0.0437, 0.0006, 0.673, 0.007, 0.953, 0.005, \
                           0.339, 0.004, 0.1747, 0.0004, 0.1108, 0.0006),
                "Asy1_00":(0.04447, 0.00016, 1, 0, 0.9555, 0.0003, 0, 0, \
                           0.1744, 0.0005, 0.1814, 0.0029),
                "Asy1_00_inv":(0, 0, 0.9555, 0.0003, 1, 0, 0.04447, 0.00016, \
                               0.1814, 0.0029, 0.1744, 0.0005),
                "Asy1_10":(0.0452, 0.0009, 0.9145, 0.0020, 0.9557, 0.0012, \
                           0.0854, 0.0015, 0.1751, 0.0029, 0.1659, 0.0029),
                "Asy1_10_inv":(0.0854, 0.0015, 0.9557, 0.0012, 0.9145, 0.0020, \
                           0.0452, 0.0009, 0.1659, 0.0029, 0.1751, 0.0029),
                "Asy1_20":(0.0433, 0.0006, 0.844, 0.007, 0.952, 0.004, \
                           0.1621, 0.0011, 0.1772, 0.0021, 0.1503, 0.0022),
                "Asy1_20_inv":(0.1621, 0.0011, 0.952, 0.004, 0.844, 0.007,\
                           0.0433, 0.0006, 0.1503, 0.0022, 0.1772, 0.0021),
                "Asy2_1": (0.05388, 0.00024, 0.9560, 0.0019, 0.9450, 0.0019, \
                           0.0458, 0.0004, 0.2127, 0.0004, 0.17443, 0.00026),
                "Asy2_2": (0.0539, 0.0011, 0.882, 0.009, 0.944, 0.003, \
                           0.1302, 0.0005, 0.2131, 0.0003, 0.1583, 0.0004),
                "Asy2_3": (0.056, 0.0009, 0.817, 0.007, 0.942, 0.004, \
                           0.195, 0.004, 0.2147, 0.0019, 0.1452, 0.0029),
                "Asy2_4": (0.0545, 0.0004, 0.7470, 0.0007, 0.9454, 0.001, \
                           0.2531, 0.0005, 0.2123, 0.0005, 0.12946, 0.00028),
                "Asy2_5": (0.05454, 0.00016, 0.6922, 0.0027, 0.9453, 0.0014, \
                           0.3084, 0.0011, 0.2124, 0.0005, 0.1168, 0.0007),
                "Asy2_6": (0.0535, 0.0010, 0.673, 0.006, 0.943, 0.007, \
                           0.3390, 0.0030, 0.2125, 0.0005, 0.1111, 0.0004),
                "Asy3_1": (0.0776, 0.001, 1, 0.009, 0.85, 0.017, \
                           0.0492, 0.0005, 0.3052, 0.0011, 0.1740, 0.0008),
                "Asy3_2": (0.0779, 0.0003, 0.949, 0.016, 0.839, 0.017, \
                           0.1343, 0.0017, 0.3044, 0.0010, 0.1600, 0.0007),
                "Asy3_3": (0.0806, 0.0006, 0.829, 0.019, 0.891, 0.022, \
                           0.209, 0.004, 0.298, 0.003, 0.1465, 0.0024),
                "Asy3_4": (0.08, 0.003, 0.7471, 0.0018, 0.894, 0.022, \
                           0.2531, 0.0011, 0.3055, 0.0026, 0.12940,	0.00028),
                "Asy3_5": (0.0816, 0.0022, 0.693, 0.005, 0.893, 0.024, \
                           0.3089, 0.0027, 0.3035, 0.0026, 0.1166, 0.0007),
                "Asy3_6": (0.079, 0.004, 0.666, 0.004, 0.894, 0.023, \
                           0.335, 0.003, 0.3011, 0.0026, 0.1105, 0.0007),
                "Asy4_1":(0.0484, 0.0004, 0.95552, 0.00025, 0.95187, 0.00018, \
                          0.0439, 0.0006, 0.1872, 0.0005, 0.17449, 0.00029),
                "Asy4_2":(0.04704, 0.00017, 0.883, 0.004, 0.9509, 0.0022, \
                          0.1257, 0.0003, 0.18661, 0.00026, 0.15858, 0.00025),
                "Asy4_3":(0.0463, 0.0002, 0.829, 0.003, 0.933, 0.007, 
                          0.1977, 0.0013, 0.1881, 0.0007, 0.1416, 0.0006),
                "Asy4_4":(0.048, 0.0003, 0.7470, 0.0007, 0.9518, 0.0007, \
                          0.2531, 0.0004, 0.1869, 0.0005, 0.12948, 0.00028),
                "Asy4_5":(0.04811, 0.00021, 0.6919, 0.0018, 0.9518, 0.0008, \
                          0.3086, 0.0013, 0.1871, 0.0005, 0.1168, 0.0007),
                "Asy4_6":(0.0476, 0.0004, 0.671, 0.008, 0.950, 0.004, \
                          0.338, 0.005, 0.1869, 0.0004, 0.1107, 0.0005),
                "Asy5_1":(0.0622, 0.0013, 0.9882, 0.017, 0.916, 0.044, \
                          0.0477, 0.0003, 0.2445, 0.0007, 0.1746, 0.0010),
                "Asy5_2":(0.0596, 0.0003, 0.899, 0.003, 0.88, 0.019, \
                          0.1272, 0.0006, 0.2461, 0.0011, 0.1582, 0.0005),
                "Asy5_3":(0.0645, 0.0018, 0.821, 0.01, 0.926, 0.013, \
                          0.195, 0.003, 0.240, 0.003, 0.1510, 0.0025),
                "Asy5_4":(0.0609, 0.0003, 0.794, 0.01, 0.877, 0.014, \
                          0.2712, 0.0022, 0.2461, 0.0006, 0.1290, 0.0006),
                "Asy5_5":(0.0603, 0.0006, 0.745, 0.007, 0.876, 0.011, \
                          0.333, 0.003, 0.24430, 0.00028, 0.1172, 0.0003),
                "Asy5_6":(0.0592, 0.0005, 0.717, 0.007, 0.872, 0.016, \
                          0.360, 0.004, 0.2449, 0.0007, 0.1113, 0.0007),
                "Asy6_1":(0.0716, 0.0008, 0.959, 0.008, 0.908, 0.014, \
                          0.0488, 0.0009, 0.2835, 0.0004, 0.17392, 0.00023),
                "Asy6_2":(0.072, 0.004, 0.8767, 0.0028, 0.9215, 0.0024, \
                          0.127, 0.004, 0.2829, 0.0006, 0.1584, 0.0005),
                "Asy6_3":(0, 0, 0, 0, 0, 0, \
                          0, 0, 0, 0, 0, 0),
                "Asy6_4":(0.0774, 0.0017, 0.7472, 0.0019, 0.920, 0.017, \
                          0.2531, 0.0007, 0.2827, 0.0006, 0.12946, 0.00028),
                "Asy6_5":(0.0776, 0.0015, 0.692, 0.005, 0.920, 0.007, \
                          0.3090, 0.0029, 0.2829, 0.0006, 0.1168, 0.0007),
                "Asy6_6":(0.072, 0.004, 0.677, 0.010, 0.907, 0.025, \
                          0.341, 0.005, 0.2828, 0.0005, 0.1106, 0.0005),
                "Asy7_1":(0, 0, 0, 0, 0, 0, \
                          0, 0, 0, 0, 0, 0),
                "Asy7_2":(0, 0, 0, 0, 0, 0, \
                          0, 0, 0, 0, 0, 0),
                "Asy7_3":(0, 0, 0, 0, 0, 0, \
                          0, 0, 0, 0, 0, 0),
                "Asy7_4":(0, 0, 0, 0, 0, 0, \
                          0, 0, 0, 0, 0, 0),
                "Asy7_5":(0, 0, 0, 0, 0, 0, \
                          0, 0, 0, 0, 0, 0),
                "Asy7_6":(0, 0, 0, 0, 0, 0, \
                          0, 0, 0, 0, 0, 0)
                }
        
    def getItemWithKey(self, key):
        "returns the value in the dictionary with some key"
        try:
            return self.parameters.get(key)
        except TypeError:
            print ("wrong key")
            return (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)