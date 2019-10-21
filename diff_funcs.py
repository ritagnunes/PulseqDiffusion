#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 14:09:08 2019

@author: rgnunes
"""

import math

import numpy as np

def get_dirs(ndirs):

    if ndirs==3:
        g= np.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
    elif ndirs==6:
        g= np.zeros((3, ndirs))        
        #obtained using gen_scheme (mrtrix)
        g= np.array([[-0.283341,  -0.893706,  -0.347862],
                     [-0.434044,   0.799575,  -0.415074],
                     [0.961905,    0.095774,   0.256058],
                     [-0.663896,   0.491506,   0.563616],
                     [-0.570757,  -0.554998,   0.605156],
                     [-0.198848,  -0.056534,  -0.978399]])                
    elif ndirs==12:
        g= np.zeros((3, ndirs))        
        #obtained using gen_scheme (mrtrix)
        g= np.array([[ 0.648514,   0.375307,   0.662249],
                     [-0.560493,   0.824711,   0.075496],
                     [-0.591977,  -0.304283,   0.746307],
                     [-0.084472,  -0.976168,   0.199902],
                     [-0.149626,  -0.006494,  -0.988721],
                     [-0.988211,  -0.056904,  -0.142130],
                     [0.864451,   -0.274379,  -0.421237],
                     [-0.173549,  -0.858586,  -0.482401],
                     [-0.039288,   0.644885,  -0.763269],
                     [0.729809,    0.673235,  -0.118882],
                     [0.698325,   -0.455759,   0.551929],
                     [-0.325340,   0.489774,   0.808873]])
    else:
        g = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    return g

def calc_bval(G, delta, Delta):
    bval= (2*math.pi*G*delta)**2*(Delta-delta/3)
    a=(Delta-delta/3)#+((gdiff_risetime)**3)/30-delta*((gdiff_risetime)**2)/6)

    return bval