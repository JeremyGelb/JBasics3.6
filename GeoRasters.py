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
    
    def Plot(self,Color="viridis") : 
        plt.imshow(self,extent = (self.Extent[0],self.Extent[2],self.Extent[1],self.Extent[3]),cmap=Color)



##############################################################
## Classe centrale MultiSpArray
#Permet l'ecriture de gros raster spatiaux compose de tuiles de sparray
##############################################################  
class MultiSpArray(object) : 
    
    def __init__(self,Extent,MaxSize,Shape,PixelShape,Src=None,Source=None) : 
        """
        Extent : une etendue specifie ainsi : (minX,minY,MaxX,MaxY)
        src : un wkt de projection tel que renvoye par osr
        NB : Shape=(Row,Col)
        NB : PixelShape = (Width,Height)
        """
        #1) setting des parametres
        self.Extent = Extent
        self.Shape = Shape
        self.PixelShape = PixelShape
        self.Source = Source
        self.Src=Src
        self.MaxSize=MaxSize
        #2) construire la structures des tuiles
        NbWidth = int(self.Shape[1]/MaxSize)
        if NbWidth*MaxSize<self.Shape[1] : 
            NbWidth+=1
        NbHeight = int(Shape[0]/MaxSize)
        if NbHeight*MaxSize<self.Shape[0] : 
            NbHeight+=1
        TotCol=0
        #iterer sur les colonnes pour delimiter les tuiles
        self.TileIndex = {}
        for i in range(NbWidth) : 
            if TotCol+MaxSize<=self.Shape[1] : 
                AddCol = MaxSize
            else : 
                AddCol = self.Shape[1]-TotCol
            StartCol = TotCol
            EndCol = TotCol + AddCol
            MinX = self.Extent[0]+self+PixelShape[0]*StartCol
            MaxX = self.Extent[0]+self+PixelShape[0]*EndCol
            #iterer sur les lignes pour delimiter les tuiles
            TotRow = 0
            for j in range(NbHeight) : 
                if TotRow+MaxSize<=self.Shape[0] : 
                    AddRow = MaxSize
                else : 
                    AddRow = self.Shape[0]-TotRow
                StartRow = TotRow
                EndRow = TotRow+AddRow
                MinY = self.Extent[3]-self.PixelShape[1]*StartRow
                MaxY = self.Extent[3]-self.PixelShape[1]*EndRow
                self.TileIndex[(i,j)] = {"Extent":[MinX,MinY,MaxX,MaxY],"Pxs":[StartRow,EndRow,StartCol,EndCol]}
                TotRow+=AddRow
            TotCol+=AddCol
        #Creer finalement l'index spatial
        #BBox demande par QIndex : (xmin,ymin,xmax,ymax)
        self.SpIndex = QdIndex(self.Extent)
        for key,value in self.TileIndex.values() : 
            self.SpIndex.insert(key,value["Extent"])
        
        
        

   
################################
#Tests !
        
if __name__=="__main__" : 
    Values = [[0,1,2,3],[4,5,6,7]]
    Extent = [0,0,20,20]
    Src = "blabla"
    Raster = SpArray(Values,Src,Extent)
    print(Raster.Sample([1,1]))
    print(Raster.Sample([19,19]))
    Raster.Plot("magma")