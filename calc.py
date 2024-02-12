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
import sympy as smp

class Diff():
    """
    class, encompassing all the calculating algs, 
    init takes following args:
        h - space step
        tau - time step
        m - vertical size in number of nodes
        n - horizontal size in number of nodes
        a - array of source function coeffs
        D - nparray of diffusion coefficients 
        func - list of sympy function expressions
        coefs - list of coefficients
        
    """
    def __init__(self, h, tau, m, n, D = np.ones(1), a = np.ones(1), func = ['0']):
        if len(func) != len(D):
            print('dimensions do not match')
            raise TypeError
        self.nComps = len(D)
        self.N = m*n
        self.m = m
        self.n = n
        self.h = h 
        self.tau = tau
        self.crnt = tau/(h*h) * D  # Courant values
        self.A = self.createA()
        #print('A=', self.A)
        self.B = self.createB()
        #print('B=', self.B)
        self.a = a
        
        for _ in range(self.nComps):
            pass
        
        # to do !!!!preparing lambda functions for source evaluation
        #row = smp.lambdify((x, y), smp.Matrix((x, x + y, x*y)), modules='numpy')
        

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
        #hh = self.h*self.h
        ex = np.ones(self.N)
        offsets = np.array([-self.n, -1,0,1,self.n])
        diag = ex + 2*self.crnt
        # 1st is lost:
        right = -0.5*self.crnt*np.array([0] + [(1 + 0.5/(i%self.n)) if i%self.n else 0 for i in range(self.N-1)])
       
        # last is lost:
        left = -0.5*self.crnt*np.array([(1 - 0.5/(i%self.n)) if i%self.n else 0 for i in range(1,self.N)] + [0])
        # leftmost boarder:
        lmost = np.arange(self.m)*self.n
        #print(lmost)
        diag[lmost] = 1 + 3*self.crnt
        #print(diag)
        left[lmost[1:]-1] = 0
        right[lmost+1] = -2*self.crnt
         # rightmost boarder:
        rmost = lmost + self.n - 1 
        left[rmost-1] = -self.crnt
        right[rmost[:-1]+1] = 0
        # first n are lost
        lower = -ex*0.5*self.crnt
        # bottom boarder:
        lower[self.N-2*self.n : self.N-self.n] *= 2 
        #last n are lost
        upper = -ex*0.5*self.crnt
        # top boarder:
        upper[self.n : 2*self.n] *= 2
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
        #hh = self.h*self.h
        ex = np.ones(self.N)
        offsets = np.array([-self.n, -1,0,1,self.n])
        diag = ex - 2*self.crnt
        # 1st is lost:
        right = 0.5*self.crnt*np.array([0] + [(1 + 0.5/(i%self.n)) if i%self.n else 0 for i in range(self.N-1)])
       
        # last is lost:
        left = 0.5*self.crnt*np.array([(1 - 0.5/(i%self.n)) if i%self.n else 0 for i in range(1,self.N)] + [0])
        # leftmost boarder:
        lmost = np.arange(self.m)*self.n
        diag[lmost] = 1 - 3*self.crnt
        left[lmost[1:]-1] = 0
        right[lmost+1] = 2*self.crnt
         # rightmost boarder:
        rmost = lmost + self.n - 1 
        left[rmost-1] = self.crnt
        right[rmost[:-1]+1] = 0
        # first n are lost
        lower = ex*0.5*self.crnt
        # bottom boarder:
        lower[self.N-2*self.n : self.N-self.n] *= 2 
        #last n are lost
        upper = ex*0.5*self.crnt
        # top boarder:
        upper[self.n : 2*self.n] *= 2
        # first n are lost
        data = np.array([lower, left, diag, right, upper ],dtype='float64')
        return ssp.dia_array((data, offsets), shape=(self.N, self.N)).asformat('csr')

    def source(self, xold):
        """
        is responsible for adding source components to both sides of the difference equation
        \partial_t u = \partialD u + f(u)
        Au_new + f(u_new) = Bu_old + f(u_old)
        Au_new + (\p f / \p u_old) u_new = Bu_old + (\p f / \p u_old) u_old

        """
        
        partial = self.a - 2*self.a*xold    # to change!
        return (partial, func)
    
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
        #left = self.A + ssp.eye(self.m*self.n, dtype = 'float64', format = 'csr').multiply(self.source(xold))
        #right = self.B@xold + ssp.eye(self.m*self.n, dtype = 'float64', format = 'csr').multiply(self.source(xold))@xold
        
        left = self.A 
        right = self.B@xold + self.a * xold - self.a * np.power(xold, 2)
        return ssp.linalg.spsolve(left, right)
         

class Diff1d():
    def __init__(self, h, tau, n, fdata):


    def createA(self):
        pass

    def createB(self):
        pass

    def step(self):
        
