#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
about calculating the next time slice using Crank-Nicholson scheme, 
cylindrical coordinate laplacian , points being numbered as follows:
 
 *    *    *    * ... *    *    *
 n(m-1)                         nm-1
 ............
 *    *    *    * ... *    *    *
           i+n 
 *    *    *    * ... *    *    *
      i-1  i   i+1
 *    *    *    * ... *    *    *
           i-n
..............
 *    *    *    * ... *    *    *
 0    1    2    3 ... n-3  n-2  n-1
@author: ariel
"""

import numpy as np
import scipy.sparse as ssp

def createA(h, tau, m, n):
    N=m*n
    ex = np.ones(N)
    offsets = np.array([-n, -1,0,1,n])
    ipos = np.repeat(np.array([0] + [(i + 0.5)/(2*h*h*(i)*(i)) for i in range(1,n-1)] + [0]),m)

    ineg = np.repeat(np.array([0] + [(i + 0.5)/(2*h*h*(i)*(i)) for i in range(1,n-1)] + [0]),m)

    data = np.array([-ex*1/(2*h*h), -ineg, ex*(1/tau + 2/(h*h)), -ipos, -ex*1/(2*h*h) ],dtype='float64')

    return ssp.dia_array((data, offsets), shape=(N, N))
   

def create_b(h, tau):
    return np.ones(h*tau)

def step(h, tau):
    return ssp.linalg.spsolve(createA(h, tau), create_b(h, tau))