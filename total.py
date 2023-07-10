import numpy as np
import vtk
from os import mkdir, chdir
from time import asctime


# Класс расчётной сетки
class CalcMesh:

    # Конструктор сетки size x size точек с шагом h по пространству
    def __init__(self, vsize, hsize, step):
        # 2D-сетка из расчётных точек, у каждой из которых, тем не менее, 3 координаты
        self.hsize = hsize
        self.vsize = vsize
        xys = np.array(np.meshgrid(np.arange(hsize),np.arange(vsize)), dtype='float64')
        
        
        xys = np.multiply( xys, step, casting="unsafe")
        cds=xys.reshape(xys.shape[0],xys.shape[1]*xys.shape[2])
        
        self.nodes = np.vstack((cds, np.zeros(shape=cds.shape[1:]))).T
        print(self.nodes.shape)
        print(self.nodes.dtype)

      
        print('coords=\n ', self.nodes)
        #-------------------------------------------------------------
        
        
        
       

      

        # Модельная скалярная величина распределена как-то вот так
        self.smth = np.zeros((self.hsize, self.vsize))
        
        
        # Профиль скорости взят какой-нибудь с потолка, есть только компонента по Z
        #self.velocity = np.zeros(shape=(3, size, size), dtype=np.double)
        self.velocity = np.power(self.nodes[0] - self.nodes[1], 2) # (x-y)^2
        
        
    # Метод отвечает за выполнение для всей сетки шага по времени величиной tau
    def move(self, tau):
        # По сути метод просто двигает все точки c их текущими скоростями
        self.smth += self.velocity[:,:] * tau

    # Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    def snapshot(self, snap_number):
        # Сетка в терминах VTK

        # Скалярное поле на точках сетки
        smth = vtk.vtkDoubleArray()
        smth.SetName("smth")

        # Векторное поле на точках сетки
        """
        vel = vtk.vtkDoubleArray()
        vel.SetNumberOfComponents(4)
        vel.SetName("vel")
        """
        
        # Обходим все точки нашей расчётной сетки
        # Делаем это максимально неэффективным, зато наглядным образом
        # nodes - point coords, 3-by-height-by-width, vtkPoints
        # smth - scalar values, 1-by-height-by-width, vtkDoubleArray
        # vel - vector values, 3-by-height-by-width, vtkDoubleArray
        
        # the following is all about attaching data to points    
        
        #for i in range(0, number):
        #for j in range(0, self.nodes.shape[0]):
                # Вставляем новую точку в сетку VTK-снапшота

         #       points.InsertNextPoint(self.nodes[i,0],self.nodes[i,1],self.nodes[i,2])
                # Добавляем значение скалярного поля в этой точке
                #smth.InsertNextValue(self.smth[i,j])
                # Добавляем значение векторного поля в этой точке
                #vel.InsertNextTuple((self.velocity[0][i,j], self.velocity[1][i,j], self.velocity[2][i,j], self.velocity[3][i,j]))
        vtkarr = vtk.vtkDoubleArray()
        vtkarr.SetNumberOfComponents(self.nodes.shape[1])
        vtkarr.SetNumberOfTuples(self.nodes.shape[0])
        vtkarr.SetVoidArray(self.nodes.copy(order='C'), self.nodes.size, 1)
        print(vtkarr)
        
        points = vtk.vtkPoints()
        points.SetData(vtkarr)
        print(points)
        
        print(self.smth.shape)
        print(self.smth)
        values = np.reshape(self.smth,(self.smth.shape[0]*self.smth.shape[1], 1))
        smth.SetNumberOfComponents(values.shape[1])
        smth.SetNumberOfTuples(values.shape[0])
        smth.SetVoidArray(values, values.size, 1)
        
        structuredGrid = vtk.vtkStructuredGrid()
        structuredGrid.SetDimensions(self.hsize, self.vsize, 1)
        structuredGrid.SetPoints(points)
        
        
        print(structuredGrid)
        
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
    vsize = 120
    hsize = 50
    # Шаг точек по пространству
    h = 0.1
    # Шаг по времени
    tau = 0.01

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
    #for i in range(1, 100):
    #    m.move(tau)
    #    m.snapshot(i)

    chdir("..")
    #---------------------------------------------------------------
