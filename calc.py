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

class Diff():
    def __init__(self, h, tau, m,n):
        self.N = m*n
        self.m = m
        self.n = n
        self.h = h 
        self.tau = tau
        self.A = self.createA()
        print(self.A)
        self.B = self.createB()
        print(self.B)

    def createA(self):
        """
        returns left-side matrix A in Ax_new = Bx_old + ... in dia_array form
        ======================================================================
        diag[0] right[1] 0     0   ... upper[n] 0 ...
        left[0] diag[1]  right[2] 0   ... 0     upper[n+1] ...
        0       left[1]  diag[2]  right[3]   .....
        ......................................
        ..........  lower[N-n-2] 0     ....    left[N-3] diag[N-2] right[N-1]
        .............0           lower[N-n-1] ...0       left[N-2] diag[N-1]   
        ======================================================================
        Parameters
        ----------
        h : float
            space step.
        tau : float
            time step.
        m : int
            grid x size.
        n : int
            grid y size.
    
        """
        hh = self.h*self.h
        ex = np.ones(self.N)
        offsets = np.array([-self.n, -1,0,1,self.n])
        diag = ex*(1/self.tau + 2/hh)
        # 1st is lost:
        right = -np.array([0] + [(1 + 0.5/(i%self.n))/(2*hh) if i%self.n else 0 for i in range(self.N-1)])
       
        # last is lost:
        left = -np.array([(1 - 0.5/(i%self.n))/(2*hh) if i%self.n else 0 for i in range(1,self.N)] + [0])
        # leftmost boarder:
        lmost = np.arange(self.m)*self.n
        #print(lmost)
        diag[lmost] = 1/self.tau + 3/(hh)
        #print(diag)
        left[lmost[1:]-1] = 0
        right[lmost+1] = -2/(hh)
         # rightmost boarder:
        rmost = lmost + self.n - 1 
        left[rmost-1] *= 2
        right[rmost[:-1]+1] = 0
        # first n are lost
        lower = -ex*1/(2*hh)
        # bottom boarder:
        lower[self.N-2*self.n:self.N-self.n] *= 2 
        #last n are lost
        upper = -ex*1/(2*hh)
        # top boarder:
        upper[self.n:2*self.n] *= 2
        # first n are lost
        data = np.array([lower, left, diag, right, upper ],dtype='float64')
    
        return ssp.dia_array((data, offsets), shape=(self.N, self.N)).asformat('csr')
       
    
    def createB(self):
        """
        returns right-side matrix B in Ax_new = Bx_old + ... in dia_array form
        ======================================================================
        diag[0] right[1] 0     0   ... upper[n] 0 ...
        left[0] diag[1]  right[2] 0   ... 0     upper[n+1] ...
        0       left[1]  diag[2]  right[3]   .....
        ......................................
        ..........  lower[N-n-2] 0     ....    left[N-3] diag[N-2] right[N-1]
        .............0           lower[N-n-1] ...0       left[N-2] diag[N-1]   
        ======================================================================
        Parameters
        ----------
        h : float
            space step.
        tau : float
            time step.
        m : int
            grid x size.
        n : int
            grid y size.
    
        """
        hh = self.h*self.h
        ex = np.ones(self.N)
        offsets = np.array([-self.n, -1,0,1,self.n])
        diag = ex*(1/self.tau - 2/(hh))
        # 1st is lost:
        right = np.array([0] + [(1 + 0.5/(i%self.n))/(2*hh) if i%self.n else 0 for i in range(self.N-1)])
       
        # last is lost:
        left = np.array([(1 - 0.5/(i%self.n))/(2*hh) if i%self.n else 0 for i in range(1,self.N)] + [0])
        # leftmost boarder:
        lmost = np.arange(self.m)*self.n
        diag[lmost] = 1/self.tau - 3/(hh)
        left[lmost[1:]-1] = 0
        right[lmost+1] = -2/(hh)
         # rightmost boarder:
        rmost = lmost + self.n - 1 
        left[rmost-1] *= 2
        right[rmost[:-1]+1] = 0
        # first n are lost
        lower = ex*1/(2*hh)
        # bottom boarder:
        lower[self.N-2*self.n:self.N-self.n] *= 2 
        #last n are lost
        upper = ex*1/(2*hh)
        # top boarder:
        upper[self.n:2*self.n] *= 2
        # first n are lost
        data = np.array([lower, left, diag, right, upper ],dtype='float64')
        return ssp.dia_array((data, offsets), shape=(self.N, self.N)).asformat('csr')
    
    def step(self, xold ):
        """
        gives x_new known x_old with Ax_new = Bx_old + ... 
        
    
        Parameters
        ----------
        h : float
            space step.
        tau : float
            time step.
        m : int
            grid x size.
        n : int
            grid y size.
    
        Returns
        -------
        ndarray
            new step values.
    
        """
        return ssp.linalg.spsolve(self.A, self.B@xold)