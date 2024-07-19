import numpy as np
import vtk
from os import mkdir, chdir
from time import asctime
import calc as cc 
from math import floor, ceil
from collections import namedtuple as nmtuple
import sympy as smp

# Класс расчётной сетки
class CalcMesh:

    # grid vsize x hsize with h space step
    
    def __init__(self, vsize, hsize, step, tau, names = None, ncomps = 1):
        """

        Parameters
        ----------
        vsize : int
            
        hsize : int
            
        step : float
            
        tau : float
            
        names : str
            names of components sep-d by spaces.
        ncomps : int
            Number of components. The default is 1.

        Returns
        -------
        CalcMesh class instance.

        """
        # 2D-сетка из расчётных точек, у каждой из которых, тем не менее, 3 координаты
        self.hsize = hsize
        self.vsize = vsize
        self.step = step
        self.tau = tau
        if names:
            self.ncomps = len(names.split())
        else:
            self.ncomps = ncomps
        xys = np.array(np.meshgrid(np.arange(hsize),np.arange(vsize)), dtype='float64')
        
        
        xys = np.multiply( xys, step, casting="unsafe")
        self.cds=xys.reshape(xys.shape[0],xys.shape[1]*xys.shape[2])
        
        self.nodes = np.vstack((self.cds, np.zeros(shape=self.cds.shape[1:]))).T
        
        if names:
            self.names = names.split()
        else:
            self.names = ' '.join(['u' + str(i+1) for i in range(ncomps)])
        #self.comp = cc.Diff(self.step, self.tau, self.vsize, self.hsize)
        #-------------------------------------------------------------
    
    def startdistrib(self, funcs):
        if len(funcs)!=self.ncomps:
            raise ValueError('number of funcs should match num of comps')
            exit()
        crds = smp.symbols('r z')
        fs = smp.Matrix(funcs).T
        u0 = smp.lambdify(crds, fs, modules = 'numpy')
        return u0(self.cds.T)
        
    #  tau size time step for the whole grid 
    def move(self, num, tau):
    
        calc = cc.Diff(self.step, tau, self.vsize, self.hsize, )
        #print('left=', calc.A.nonzero())
        #print('right=', calc.B.nonzero())
        
        self.smth = calc.step(self.smth)
    
    def snapshot(self, snap_number):
        # data arrays
        smth = []
        for name in self.names:
            arr = vtk.vtkDoubleArray()
            arr.SetName(name)
            smth.append(arr)
        #-------------------------    
        
        #geometric points
        points = vtk.vtkPoints()
       
        for x in self.nodes:
            points.InsertNextPoint(x)
        #-------------------------------
        
        structuredGrid = vtk.vtkStructuredGrid()
        structuredGrid.SetDimensions(self.hsize, self.vsize, 1)
        structuredGrid.SetPoints(points)
        
        # filling vtk data arrays
        n = len(self.smth)
        print(n)
        for i in range(n):
            values = np.reshape(self.smth[i], (self.smth[i].size), order = 'C')     # self.smth[i] must be of size n-by-1
            print(i)
            print(type(smth[i]))
            smth[i].SetNumberOfComponents(1)
            smth[i].SetNumberOfTuples(values.shape[0])
            smth[i].SetVoidArray(np.ascontiguousarray(values), values.size, 1)
            structuredGrid.GetPointData().AddArray(smth[i])     # appending data arrays
        #-----------------------------------------
        
        #--------------------------------------------
        
        #print(structuredGrid)
        
        # Присоединяем векторное и скалярное поля к точкам
        
        #self.structuredGrid.GetPointData().AddArray(vel)

        #------------------------------------------------------------

        # Создаём снапшот в файле с заданным именем

        writer = vtk.vtkXMLStructuredGridWriter()
        writer.SetInputDataObject(structuredGrid)
        writer.SetFileName("cubic2d-step-" + str(snap_number) + ".vts")
        writer.Write()

# setting up starting distribution       
def startdistrib(m, n, names = 'u1'):
    """

    Parameters
    ----------
    m : vertical size, int.
    n : horizontal size, int
    names : str of matter names sep-ted by spaces
    
    Returns
    -------
    matts : namedtuple of numpy ndarrays of size (m*n, ) , floats

    """
    samp = np.concatenate((np.zeros(floor(n * 0.75)), np.ones(ceil(n * 0.25))))
    x = np.tile(samp, m)
    data = [x]
    M = ntup('u', names)
    us = M(data)
    return us

if __name__=="__main__":
    # Размер расчётной сетки, точек на сторону

    vsize = 13
    hsize = 13
    # Шаг точек по пространству
    h = 0.1
    # Шаг по времени
    tau = 0.002

    

    nsteps = 20

    #init distributions
    names = 'a1 a2'
    Dist = nmtuple('dist', 'exprs conds')
    a1 = Dist(['0'], [True])
    a2 = Dist(['1'], [True])

    Distrs = nmtuple('u0', names)
    dist = Distrs(*[a1, a2]) #make a tuple of all dists to pass to startdistrib



    #--------------------------------

    # Создаём сетку заданного размера
    m = CalcMesh(vsize, hsize, h, tau, names)
    #---------------------------------
    
    # setting up starting distibutions
    #at this point the mesh is introduced with its data
    #startdistrib returns a namedtuple 
    m.smth = m.startdistrib(dist)
    #---------------------------------------------

    # create dir for specific cartoon, move to it, leave it at the end

    dirname = 'cartoon' + '-'.join(asctime()[3:-5].replace(':','=').split())
    mkdir(dirname)
    chdir(dirname)

    # Пишем её начальное состояние в VTK
    m.snapshot(0)
    np.savetxt('nodes',m.nodes)



    # Делаем шаги по времени,
    # на каждом шаге считаем новое состояние и пишем его в VTK

    for i in range(1, nsteps):

        m.move(i, tau)
        m.snapshot(i)

    chdir("..")
    #---------------------------------------------------------------
