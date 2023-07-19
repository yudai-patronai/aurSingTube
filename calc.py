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
    """
    returns left-side matrix A in Ax_new = Bx_old + ... in dia_array form
    """
    N=m*n
    ex = np.ones(N)
    offsets = np.array([-n, -1,0,1,n])
    diag = ex*(1/tau + 2/(h*h))
    # 1st is lost:
    right = -np.array([(i%n + 0.5)/(2*h*h*(i%n)*(i%n)) for i in range(N)])
    # last is lost:
    left = -np.array([(i%n - 0.5)/(2*h*h*(i%n)*(i%n)) for i in range(N)])
    lower = -ex*1/(2*h*h)
    # bottom boarder:
    lower[N-n:] += -np.ones(n)*1/(2*h*h)
    upper = -ex*1/(2*h*h)
    # top boarder:
    upper[n:2n] += -np.ones(n)*1/(2*h*h)
    data = np.array([lower, left, diag, right, upper ],dtype='float64')

    return ssp.dia_array((data, offsets), shape=(N, N))
   

def create_b(h, tau):
    return np.ones(h*tau)

def step(h, tau):
    return ssp.linalg.spsolve(createA(h, tau), create_b(h, tau))