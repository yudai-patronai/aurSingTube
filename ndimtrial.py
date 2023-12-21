import numpy as np
import total as tt

m = 10
n = 10
u = [np.zeros((m*n)),   np.ones((m*n))]
mesh = tt.CalcMesh(m, n, 0.1, 0.1, 'u1 u2')
mesh.smth = u
mesh.snapshot(0)