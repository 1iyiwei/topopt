#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 17:07:13 2018

@author: bram
"""

import numpy as np
import cv2 as cv
import math

# filter
def filt(x, rmin, dc):
    rminf = math.floor(rmin)

    dcn = np.zeros(x.shape)
    nely, nelx = 11, 11

    for i in range(nelx):
        for j in range(nely):
            sum = 0.0
            for k in range(max(i-rminf, 0), min(i+rminf+1, nelx)):
                for l in range(max(j-rminf, 0), min(j+rminf+1, nely)):
                    weight = max(0, rmin - np.sqrt((i-k)**2+(j-l)**2));
                    sum = sum + weight;
                    dcn[j,i] = dcn[j,i] + weight*x[l,k]*dc[l,k];
        
            dcn[j,i] = dcn[j,i]/(x[j,i]*sum);

    return dcn

# new filter based upon C++ accelerated code
def filtNew(x, rmin, dc):
    rminf = math.floor(rmin)
    nely, nelx = x.shape

    # define normalized convolution kernel based upon rmin
    size = rminf*2+1
    kernel = np.zeros((size, size))
    for i in range(size):
        for j in range(size):
            dis = np.sqrt((rminf-i)**2 + (rminf-j)**2)
            kernel[i, j] = np.max((0, rmin - dis))
    kernel = kernel/np.sum(kernel)  # normalisation

    # elementwise multiplication of x and dc
    xdc = x*dc

    # execute convolution on xdc using opencv, with reflected 101 border
    xdcn = cv.filter2D(xdc, -1, kernel, borderType=cv.BORDER_REFLECT_101)

    # then dcn = xdcn/x (elementwise)
    dcn = xdcn/x

    return dcn


rmin = 1.5
x = np.ones((11, 11))
dc = np.ones((11, 11))
dc[4, 0] = 10
dc[5, 1] = 10
dc[3, 1] = 10

dcn = filt(x, rmin, dc)
dcn2 = filtNew(x, rmin, dc)

if np.allclose(dcn, dcn2, rtol=1e-10):
    print('succes')
    