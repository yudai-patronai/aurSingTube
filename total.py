import numpy as np
import vtk
from os import mkdir, chdir
from time import asctime
import calc as cc 
from math import floor, ceil

# Класс расчётной сетки
class CalcMesh:

    # grid vsize x hsize with h space step
    
    def __init__(self, vsize, hsize, step, tau, ncomps = 1):
        # 2D-сетка из расчётных точек, у каждой из которых, тем не менее, 3 координаты
        self.hsize = hsize
        self.vsize = vsize
        self.step = step
        self.tau = tau
        self.ncomps = ncomps
        xys = np.array(np.meshgrid(np.arange(hsize),np.arange(vsize)), dtype='float64')
        
        
        xys = np.multiply( xys, step, casting="unsafe")
        cds=xys.reshape(xys.shape[0],xys.shape[1]*xys.shape[2])
        
        self.nodes = np.vstack((cds, np.zeros(shape=cds.shape[1:]))).T
        print(self.nodes.shape)
        print(self.nodes.dtype)

      
        print('coords=\n ', self.nodes)
        
        #self.comp = cc.Diff(self.step, self.tau, self.vsize, self.hsize)
        #-------------------------------------------------------------
        
 
    #  tau size time step for the whole grid 
    def move(self, num, tau):
    
        calc = cc.Diff(self.step, tau, self.vsize, self.hsize)
        #print('left=', calc.A.nonzero())
        #print('right=', calc.B.nonzero())
        
        self.smth = calc.step(self.smth)
    
    def snapshot(self, snap_number):
       
        smth = vtk.vtkDoubleArray()
        smth.SetName("smth")
        
        #geometric points
        points = vtk.vtkPoints()
       
        for x in self.nodes:
            points.InsertNextPoint(x)
        #-------------------------------
        
        values = np.reshape(self.smth,(self.smth.shape[0], 1))
        smth.SetNumberOfComponents(values.shape[1])
        smth.SetNumberOfTuples(values.shape[0])
        smth.SetVoidArray(values, values.size, 1)
        
        structuredGrid = vtk.vtkStructuredGrid()
        structuredGrid.SetDimensions(self.hsize, self.vsize, 1)
        structuredGrid.SetPoints(points)
        
        
        #print(structuredGrid)
        
        # Присоединяем векторное и скалярное поля к точкам
        structuredGrid.GetPointData().AddArray(smth)
        #self.structuredGrid.GetPointData().AddArray(vel)

        #------------------------------------------------------------

        # Создаём снапшот в файле с заданным именем

        writer = vtk.vtkXMLStructuredGridWriter()
        writer.SetInputDataObject(structuredGrid)
        writer.SetFileName("cubic2d-step-" + str(snap_number) + ".vts")
        writer.Write()
# setting up starting distribution       
def startdistrib(m, n, ncomps = 1):
    """

    Parameters
    ----------
    m : vertical size, int.
    n : horizontal size, int

    Returns
    -------
    x : numpy ndarray of size (m*n, ) , floats

    """
    samp = np.concatenate((np.zeros(floor(n * 0.75)), np.ones(ceil(n * 0.25))))
    x = np.tile(samp, m)
    return x

if __name__=="__main__":
    # Размер расчётной сетки, точек на сторону

    vsize = 13
    hsize = 13
    # Шаг точек по пространству
    h = 0.1
    # Шаг по времени
    tau = 0.002

    

    nsteps = 50


    # Создаём сетку заданного размера
    m = CalcMesh(vsize, hsize, h, tau)
    #---------------------------------
    
    # setting up starting distibutions
    
    
    m.smth = startdistrib(vsize, hsize)
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
