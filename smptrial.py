# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#%%
import numpy as np
import sympy as smp

frm = ['k1*u1', 'k2*u1*u2', 'k3*(u1-u3)']
args = smp.symbols('k1 k2 k3 u1 u2 u3 ')
x,y = smp.symbols('x y')
##row = smp.lambdify((x, y), smp.Matrix((x, x + y, x*y)), modules='numpy')
vec = smp.Matrix(frm).T
dvec = smp.diff(vec, x)
row = smp.lambdify(args, vec, modules='numpy')
drow = smp.lambdify(args, dvec, modules='numpy')

u1 = np.ones(3)
u2 = np.array([0,1,1], dtype='float')
u3 = np.array([1,1,0], dtype='float')

k1 = 1.0
k2 = 0.5
k3 = -1.0

unew = row(k1,k2,k3,u1,u2,u3)

#%%
import numpy as np
import sympy as smp


x = smp.symbols('x', real = True)
f = smp.simplify(smp.Piecewise((-x, x < 0), (x, x > 0), (0, x==0)))
func  = smp.lambdify(x, f, modules='numpy')

args = np.arange(-1, 1, 0.1)

print(func(args))

#print(f.subs(x, -1))  # Output: 1
#print(f.subs(x, 1))   # Output: 2


