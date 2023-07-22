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
    ======================================================================
    diag[0] right[1] 0     0   ... upper[n] 0 ...
    left[0] diag[1]  right[2] 0   ... 0     upper[n+1] ...
    0       left[1]  diag[2]  right[3]   .....
    ......................................
    ..........  lower[N-2n-2] 0     ....    left[N-3] diag[N-2] right[N-1]
    .............0           lower[N-2n-1] ...0       left[N-2] diag[N-1]   
    ======================================================================
    """
    N=m*n
    ex = np.ones(N)
    offsets = np.array([-n, -1,0,1,n])
    diag = ex*(1/tau + 2/(h*h))
    # 1st is lost:
    right = -np.array([0] + [(i%n + 0.5)/(2*h*h*(i%n)*(i%n)) if i%n else 0 for i in range(N-1)])
   
    # last is lost:
    left = -np.array([(i%n - 0.5)/(2*h*h*(i%n)*(i%n)) if i%n else 0 for i in range(1,N)] + [0])
    # leftmost boarder:
    lmost = np.arange(m)*n
    print(lmost)
    diag[lmost] = 1/tau + 3/(h*h)
    print(diag)
    left[lmost[1:]-1] = 0
    right[lmost+1] = -2/(h*h)
     # rightmost boarder:
    rmost = lmost + n - 1 
    left[rmost-1] *= 2
    right[rmost[:-1]+1] = 0
    # first n are lost
    lower = -ex*1/(2*h*h)
    # bottom boarder:
    lower[N-2*n:N-n] *= 2 
    #last n are lost
    upper = -ex*1/(2*h*h)
    # top boarder:
    upper[n:2*n] *= 2
    # first n are lost
    data = np.array([lower, left, diag, right, upper ],dtype='float64')

    return ssp.dia_array((data, offsets), shape=(N, N))
   

def create_b(h, tau, m, n):
    """
    returns right-side matrix B in Ax_new = Bx_old + ... in dia_array form
    ======================================================================
    diag[0] right[1] 0     0   ... upper[n] 0 ...
    left[0] diag[1]  right[2] 0   ... 0     upper[n+1] ...
    0       left[1]  diag[2]  right[3]   .....
    ......................................
    ..........  lower[N-2n-2] 0     ....    left[N-3] diag[N-2] right[N-1]
    .............0           lower[N-2n-1] ...0       left[N-2] diag[N-1]   
    ======================================================================
    """
    N=m*n
    ex = np.ones(N)
    offsets = np.array([-n, -1,0,1,n])
    diag = ex*(1/tau - 2/(h*h))
    # 1st is lost:
    right = np.array([0] + [(i%n + 0.5)/(2*h*h*(i%n)*(i%n)) if i%n else 0 for i in range(N-1)])
   
    # last is lost:
    left = np.array([(i%n - 0.5)/(2*h*h*(i%n)*(i%n)) if i%n else 0 for i in range(1,N)] + [0])
    # leftmost boarder:
    lmost = np.arange(m)*n
    diag[lmost] = 1/tau - 3/(h*h)
    left[lmost[1:]-1] = 0
    right[lmost+1] = -2/(h*h)
     # rightmost boarder:
    rmost = lmost + n - 1 
    left[rmost-1] *= 2
    right[rmost[:-1]+1] = 0
    # first n are lost
    lower = ex*1/(2*h*h)
    # bottom boarder:
    lower[N-2*n:N-n] *= 2 
    #last n are lost
    upper = ex*1/(2*h*h)
    # top boarder:
    upper[n:2*n] *= 2
    # first n are lost
    data = np.array([lower, left, diag, right, upper ],dtype='float64')
    return np.ones(h*tau)

def step(h, tau):
    return ssp.linalg.spsolve(createA(h, tau), create_b(h, tau))