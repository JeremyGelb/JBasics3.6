# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 18:57:36 2019

@author: gelbj
"""

from osgeo import gdal, osr
import shapely
import numpy as np
import GeomTools as GT
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection, LineCollection
from collections import defaultdict

from libs.QuadTree import Index as QdIndex

##############################################################
## Definition des erreurs que l'on peut rencontrer
##############################################################
class ExtentError(Exception):
    """
    Raise me when we ask somethong out of the bound of an object
    """
    pass


##############################################################
## Classe centrale SpArray
#etend le classique np.array pour lui donner une dimension spatiale
##############################################################    

class SpArray(np.ndarray):
    """
        Extent : une etendue specifie ainsi : (minX,minY,MaxX,MaxY)
        src : un wkt de projection tel que renvoye par osr
        NB : shape=(Row,Col)
        NB : PixelShape = (Width,Height)
    """
    def __new__(cls, data, Src,Extent):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        if type(data)!=np.ndarray :
            input_array = np.array(data)
        else : 
            input_array = data
        obj = np.asarray(input_array).view(cls)
        # add the new attribute to the created instance
        obj.Src = Src
        obj.Extent = Extent
        obj.PixelShape = ((Extent[2]-Extent[0])/input_array.shape[1],(Extent[3]-Extent[1])/input_array.shape[0])
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None: return
        self.Src = getattr(obj, 'Src', None)
        self.Extent = getattr(obj, 'Extent', None)
        #defining the other elements
        self.PixelShape = getattr(obj,'PixelShape', None)
    
    def Sample(self,Coords) : 
        """
        Renvoit la valeur du pixel au coordonnes indiquees si le pixel 
        se trouve dans l'etendue du raster, sinon, souleve une erreur
        """
        if GT.InExtent(shapely.geometry.Point(Coords),self.Extent)==False : 
            raise ExtentError("The point "+str(Coords)+" is not in the boundary of the raster : "+str(self.Extent))
        #step1 : calculer la colonne dans laquelle on se trouve
        Col =  int((Coords[0]-self.Extent[0])/self.PixelShape[0])
        Row = int((Coords[1]-self.Extent[1])/self.PixelShape[1])
        return self[Row,Col]
    
    def Plot(self,Color="viridis",Vlim=None) :
        if Vlim is None : 
            plt.imshow(self,extent = (self.Extent[0],self.Extent[2],self.Extent[1],self.Extent[3]),cmap=Color)
        else : 
            plt.imshow(self,extent = (self.Extent[0],self.Extent[2],self.Extent[1],self.Extent[3]),cmap=Color,vmin=Vlim[0],vmax=Vlim[1])


##############################################################
## Classe centrale MultiSpArray
#Permet la lecture et l'ecriture de gros raster spatiaux compose de tuiles de sparray
##############################################################  
class MultiSpArray(object) : 
    
    def __init__(self,MaxSize,Source=None,Extent=None,PixelShape=None,Crs=None,Default=None) : 
        """
        Pour creer un Multi Array, deux option possibles : 
            -Definir une source, et alors toutes les meta donnees proviennent de cette source
            -Definir les meta donnees a la main ainsi que la valeur de defaut
        
        Extent : une etendue specifie ainsi : (minX,minY,MaxX,MaxY)
        Crs : un wkt de projection tel que renvoye par osr
        NB : Shape=(Width,Height)
        NB : PixelShape = (Width,Height)
        NB : TileIndex[Colonne,Row] = {'Extent': [MinX, MinY, MaxX, MaxY], 'Pxs': [StartCol, EndCol, StartRow, EndRow]}
        NB : Shape2=(Cols,Rows)
        """
        #1) setting des parametres
        if Source is None :
            self.Extent = Extent
            self.Shape = [Extent[2]-Extent[0],Extent[3]-Extent[1]]
            self.PixelShape = PixelShape
            self.Crs=Crs
            self.Default=Default
        else : 
            DataSource = gdal.Open(Source)
            self.Crs = DataSource.GetProjection()
            ulx, xres, xskew, uly, yskew, yres  = DataSource.GetGeoTransform()
            lrx = ulx + (DataSource.RasterXSize * xres)
            lry = uly + (DataSource.RasterYSize * yres)
            self.Extent = [ulx,lry,lrx,uly]
            self.Shape = [self.Extent[2]-self.Extent[0],self.Extent[3]-self.Extent[1]]
            self.PixelShape = [xres,-yres]
            self.Shape2 = [DataSource.RasterXSize,DataSource.RasterYSize]
        
        self.Source = Source
        self.MaxSize=MaxSize
        self.__Arrays = defaultdict(lambda : None)
        #2) construire la structures des tuiles
        NbWidth = int(self.Shape[0]/MaxSize)
        if NbWidth*MaxSize<self.Shape[0] : 
            NbWidth+=1
        NbHeight = int(self.Shape[1]/MaxSize)
        if NbHeight*MaxSize<self.Shape[1] : 
            NbHeight+=1
        TotCol=0
        #iterer sur les colonnes pour delimiter les tuiles
        self.TileIndex = {}
        for i in range(NbWidth) : 
            if TotCol+MaxSize<=self.Shape[0] : 
                AddCol = MaxSize
            else : 
                AddCol = self.Shape[0]-TotCol
            StartCol = TotCol
            EndCol = TotCol + AddCol
            MinX = self.Extent[0]+self.PixelShape[0]*StartCol
            MaxX = self.Extent[0]+self.PixelShape[0]*EndCol
            #iterer sur les lignes pour delimiter les tuiles
            TotRow = 0
            for j in range(NbHeight) : 
                if TotRow+MaxSize<=self.Shape[1] : 
                    AddRow = MaxSize
                else : 
                    AddRow = self.Shape[1]-TotRow
                StartRow = TotRow
                EndRow = TotRow+AddRow
                MaxY = self.Extent[3]-self.PixelShape[1]*StartRow
                MinY = self.Extent[3]-self.PixelShape[1]*EndRow
                self.TileIndex[(i,j)] = {"Extent":[MinX,MinY,MaxX,MaxY],"Pxs":[int(StartCol),int(EndCol),int(StartRow),int(EndRow)]}
                TotRow+=AddRow
            TotCol+=AddCol
        #Creer finalement l'index spatial
        #BBox demande par QIndex : (xmin,ymin,xmax,ymax)
        self.SpIndex = QdIndex(self.Extent)
        for key,value in self.TileIndex.items() : 
            print(key)
            self.SpIndex.insert(key,value["Extent"])
            
    def GetRaster(self,Key,FromSource=True) : 
        """
        Permet d'acceder aux donnees presentes dans le raster source
        """
        if self.__Arrays[Key] is None or FromSource :
            Coordinates = self.TileIndex[Key]
            Window = Coordinates["Pxs"]
            Extent = Coordinates["Extent"]
            DataSource = gdal.Open(self.Source)
            print("key : "+str(Key))
            print("Window : "+str(Window))
            print("xoff="+str(Window[0]))
            print("yoff="+str(Window[2]))
            print("xsize="+str(Window[1]- Window[0]))
            print("ysize="+str(Window[3]-Window[2]))
            Data = DataSource.ReadAsArray(xoff=Window[0],yoff=Window[2], xsize = Window[1]- Window[0], ysize = Window[3]-Window[2])
            print(Data)
            Array = SpArray(Data,self.Crs,Extent)
            del DataSource
        else : 
            Array = self.__Arrays[Key]
        return Array
    
    
    def PlotGrid(self,LineColor="g",LineWidth=2) : 
        """
        Permet de dessiner les emprises des tuiles du raster multiple
        """
        Lines = []
        for Key,Value in self.TileIndex.items() : 
            Extent = Value["Extent"]
            Lines.append([(Extent[0],Extent[1]),(Extent[1],Extent[1]),(Extent[1],Extent[3]),(Extent[0],Extent[3]),(Extent[0],Extent[1])])
        lc = LineCollection(Lines,colors=LineColor,linewidths=LineWidth)
        fig = plt.gcf()
        ax = fig.axes[0]
        ax.add_collection(lc)
        ax.set_xlim(self.Extent[0],self.Extent[2])
        ax.set_ylim(self.Extent[1],self.Extent[3])
        plt.show()
            
            
            
        
    def Plot(self,Extent=None,FromSource=True) :
        """
        Permet de dessiner tout ou partie du raster
        """
        if Extent is None : 
            Extent = self.Extent
        
        Keys = self.SpIndex.intersect(Extent)
        f = plt.figure()
        Arrays=[]
        Mins=[]
        Maxs=[]
        for Key in Keys : 
            Raster = self.GetRaster(Key,FromSource)
            Arrays.append(Raster)
            Mins.append(np.min(Raster))
            Maxs.append(np.max(Raster))
        for Array in Arrays : 
            Array.Plot(Vlim=[np.min(Mins),np.max(Maxs)])
        a = f.axes[0]
        a.set_xlim(Extent[0],Extent[2])
        a.set_ylim(Extent[1],Extent[3])
        plt.show()
        plt.colorbar()
   
################################
#Tests !
        
if __name__=="__main__" : 
    Link = r"L:\TEMP\BatiParis_TestRaster.tif"
    Raster1 = MultiSpArray(100,Source=Link)
    Raster1.Plot()
    Raster1.PlotGrid(LineColor="r")
    