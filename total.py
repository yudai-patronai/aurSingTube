import numpy as np
import vtk
from os import mkdir, chdir
from time import asctime
import calc as cc 

# Класс расчётной сетки
class CalcMesh:

    # Конструктор сетки size x size точек с шагом h по пространству
    def __init__(self, vsize, hsize, step):
        # 2D-сетка из расчётных точек, у каждой из которых, тем не менее, 3 координаты
        self.hsize = hsize
        self.vsize = vsize
        self.step = step
        xys = np.array(np.meshgrid(np.arange(hsize),np.arange(vsize)), dtype='float64')
        
        
        xys = np.multiply( xys, step, casting="unsafe")
        cds=xys.reshape(xys.shape[0],xys.shape[1]*xys.shape[2])
        
        self.nodes = np.vstack((cds, np.zeros(shape=cds.shape[1:]))).T
        print(self.nodes.shape)
        print(self.nodes.dtype)

      
        print('coords=\n ', self.nodes)
        #-------------------------------------------------------------
        
        
        
       

      

        # Модельная скалярная величина распределена как-то вот так
        self.smth = np.ones((self.hsize*self.vsize,1))
        #self.smth[0] = 1
        
        
 
    # Метод отвечает за выполнение для всей сетки шага по времени величиной tau
    def move(self, num, tau):
        # По сути метод просто двигает все точки c их текущими скоростями
        
        calc = cc.Diff(self.step, tau, self.vsize, self.hsize)
        #print('left=', calc.A.nonzero())
        #print('right=', calc.B.nonzero())
        
        self.smth = calc.step(self.smth)
    
    def snapshot(self, snap_number):
       
        smth = vtk.vtkDoubleArray()
        smth.SetName("smth")

        
        
        points = vtk.vtkPoints()
        #points.SetData(vtkarr)
        for x in self.nodes:
            points.InsertNextPoint(x)
        #print(points)
        
        #print(self.smth.shape)
        #print(self.smth)
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

if __name__=="__main__":
    # Размер расчётной сетки, точек на сторону
    vsize = 50
    hsize = 100
    # Шаг точек по пространству
    h = 0.1
    # Шаг по времени
    tau = 0.1

    # Создаём сетку заданного размера
    m = CalcMesh(vsize, hsize, h)

    # create dir for specific cartoon, move to it, leave it at the end

    dirname = 'cartoon' + '-'.join(asctime()[3:-5].split())
    mkdir(dirname)
    chdir(dirname)

    # Пишем её начальное состояние в VTK
    m.snapshot(0)
    np.savetxt('nodes',m.nodes)



    # Делаем шаги по времени,
    # на каждом шаге считаем новое состояние и пишем его в VTK
    for i in range(1, 50):
        m.move(i, tau)
        m.snapshot(i)

    chdir("..")
    #---------------------------------------------------------------
