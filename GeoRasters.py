# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 18:57:36 2019

@author: gelbj
"""

from osgeo import gdal, osr
import shapely
from shapely import wkt
import numpy as np
import GeomTools as GT
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection, LineCollection
from collections import defaultdict
import copy

from libs.QuadTree import Index as QdIndex

##############################################################
## Definition des erreurs que l'on peut rencontrer
##############################################################
class ExtentError(Exception):
    """
    Raise me when we ask something out of the bound of an object
    """
    pass


##############################################################
## Definition d'une classe qui va gerer les extent
##############################################################
class GeoExtent(object) :
    
    def __init__(self,Xmin,Ymin,Xmax,Ymax):
        self.Xmin = Xmin
        self.Xmax = Xmax
        self.Ymin = Ymin
        self.Ymax = Ymax
        self.Width = abs(self.Xmax-self.Xmin)
        self.Height = abs(self.Ymax-self.Ymin)
        self.Geom = shapely.geometry.Polygon([(Xmin,Ymin),(Xmax,Ymin),(Xmax,Ymax),(Xmin,Ymax),(Xmin,Ymin)])
        
    def __repr__(self) : 
        return "Xmin : {} - Xmax : {} - Ymin : {} - Ymax : ".format(self.Xmin,self.Xmax,self.Ymin,self.Ymax)
        
    def List(self,Format=1) :
        """
        Format 1 : lower corner ; upper corner
        Format 2 : XX ; YY
        """
        if Format==1 : 
            return (self.Xmin,self.Ymin,self.Xmax,self.Ymax)
        elif Format == 2 : 
            return (self.Xmin,self.Xmax,self.Ymin,self.Ymax)
        else : 
            raise ValueError("This format is not supported by the List methode")
    
    def Plot(self,LineColor="r",LineWidth=1) :
        """
        Plot the extent as a line in the current plot or create a new one
        """
        try :
            fig = plt.gcf()
            ax = fig.axes[0]
        except IndexError : 
            fig,ax = plt.subplots()
        lc = LineCollection([list(self.Geom.boundary.coords)],colors=LineColor,linewidths=LineWidth)
        ax.add_collection(lc)
        plt.show()
        
        
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
        Extent = GeoExtent(*Extent)
        # We first cast to be our class type
        if type(data)!=np.ndarray :
            input_array = np.array(data)
        else : 
            input_array = data
        obj = np.asarray(input_array).view(cls)
        # add the new attribute to the created instance
        obj.Src = Src
        obj.Extent = Extent
        obj.PixelShape = (Extent.Width/input_array.shape[1],Extent.Height/input_array.shape[0])
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None: return
        self.Src = getattr(obj, 'Src', None)
        self.Extent = getattr(obj, 'Extent', None)
        #defining the other elements
        self.PixelShape = getattr(obj,'PixelShape', None)
        
    ########################################
    ## methodes d'acces aux pixels individuels
    ########################################
    
    def Sample(self,Coords) : 
        """
        Renvoit la valeur du pixel au coordonnes indiquees si le pixel 
        se trouve dans l'etendue du raster, sinon, souleve une erreur
        """
        if self.Extent.Geom.intersects(shapely.geometry.Point(Coords))==False : 
            raise ExtentError("The point "+str(Coords)+" is not in the boundary of the raster : "+str(self.Extent))
        #step1 : calculer la colonne dans laquelle on se trouve
        Col =  int((Coords[0]-self.Extent[0])/self.PixelShape[0])
        if Col == self.shape[0] : 
            Col-=1
        Row = self.shape[1]-int((Coords[1]-self.Extent[1])/self.PixelShape[1])
        if Row == self.shape[1] : 
            Row-=1
        return self[Row,Col]
    
    def PxCoords(self,Coords) : 
        """
        Renvoit la valeur du pixel au coordonnes indiquees si le pixel 
        se trouve dans l'etendue du raster, sinon, souleve une erreur
        """
        if self.Extent.Geom.intersects(shapely.geometry.Point(Coords))==False :
            raise ExtentError("The point "+str(Coords)+" is not in the boundary of the raster : "+str(self.Extent))
        #step1 : calculer la colonne dans laquelle on se trouve
        Col = (Coords[0]-self.Extent.Xmin)/self.PixelShape[0]
#        if int(Col)<Col : 
#            Col = int(Col)+1
#        else : 
#            Col = int(Col)
        Col = int(Col)
        if Col == self.shape[1] : 
            Col-=1
        Row = self.shape[0]-(Coords[1]-self.Extent.Ymin)/self.PixelShape[1]
        Row = int(Row)
        if Row == self.shape[0] : 
            Row-=1
        return [Row,Col]
    
    def PxGeom(self,PxCoords,Type="Point") : 
        """
        renvoit une geometrie du pixel, peut etre son centroid : Point, ou son polygone : Polygone
        PxCoords Row,Col
        """
        if Type=="Point" : 
            X = self.Extent.Xmin+((PxCoords[1])*self.PixelShape[0])+0.5*self.PixelShape[0]
            Y = self.Extent.Ymax-((PxCoords[0])*self.PixelShape[1])-0.5*self.PixelShape[1]
            return shapely.geometry.Point((X,Y))
        elif Type=="Polygone" : 
            x1 = self.Extent.Xmin+((PxCoords[0]-1)*self.PixelShape[0])
            x2 = self.Extent.Xmin+((PxCoords[0]-1)*self.PixelShape[0])+self.PixelShape[0]
            y1 = self.Extent.Ymax-((PxCoords[1]-1)*self.PixelShape[1])
            y2 = self.Extent.Ymax-((PxCoords[1]-1)*self.PixelShape[1])-self.PixelShape[1]
            return shapely.geometry.Polygon([[x1,y1],[x2,y1],[x2,y2],[x1,y2]])
        else : 
            raise ValueError("This type : "+Type+" is not supported")
    
    ########################################
    ## methodes d'acces aux pixels par groupe
    ########################################
    
    def ScanLine(self,Geom,Default=np.nan,Reverse=False) :
        """
        Fais passer une ligne qui remonte ligne par ligne le raster
        NB : Geom doit etre de type shapely
        retourne une selection du raster, soit un raster avec une valeur par defaut assignee a tout ce qui est au dehors
        """
        Step = self.PixelShape[1]
        #Intersected={}
        MaxIter = self.shape[0]
        StartX = self.Extent.Xmin
        EndX = self.Extent.Xmax
        StartY = self.Extent.Ymin+Step/2.0
        Copied = copy.deepcopy(self)
        As = []
        Bs = []
        #debut des iterations
        for e in range(MaxIter) :
            #creation de la ligne
            Y = StartY + Step*e
            TXT = "LINESTRING ("+str(StartX)+" "+str(Y)+","+str(EndX)+" "+str(Y)+")"
            Line = shapely.wkt.loads(TXT)
            #recuperation des points (intersections avec la geoemtries)
            Inter = Geom.intersection(Line)
            Name = Inter.geom_type
            #print(Name)
            if Name=="GeometryCollection" and len(Inter.geoms)==0 :
                continue
            elif Name=="MultiLineString" :
                Lines = list(Inter.geoms)
            elif Name=="LineString" :
                Lines = [Inter]
            elif Name=="Point" :
                Coords = Inter.coords[0]
                PxCoords = self.PxCoords(Coords)
                As.append(PxCoords[0])
                Bs.append(PxCoords[1])
                #Intersected[PxCoords] = self.Sample(Coords)
                continue
            #iteration sur ces points
            for line in Lines :
                #Bars.append(list(line.coords))
                Start,End = GT.GetExtremites(line)
                Coords1 = self.PxCoords(Start.coords[0])
                Coords2 = self.PxCoords(End.coords[0])
                #verfication pour les extremites
                Center1 = self.PxGeom(Coords1)
                Center2 = self.PxGeom(Coords2)
                if Center1.intersects(Geom) : 
                    StartCount=0
                else : 
                    StartCount=1
                if Center2.intersects(Geom) : 
                    EndCount=1
                else : 
                    EndCount=0
                
                LastPx = int(Coords2[1])-int(Coords1[1])
                for i in range(StartCount,LastPx+EndCount,1):
                    X = Coords1[1]+i
                    As.append(Coords2[0])
                    Bs.append(X)
                    #Intersected[(Coords2[0],X)]=self[Coords2[0]][X]
        mask = np.ones(self.shape,dtype=bool)
        Idxs = (tuple(As),tuple(Bs))
        mask[Idxs] = False
        if Reverse :
            Copied[~mask] = Default
        else :
            Copied[mask] = Default
        return Copied
    
    
    def SetPixels(self,PixelSelection) : 
        """
        Permet d'assigner des valeurs aux pixels
        """
        Xs = []
        Ys = []
        Values = []
        for key,value in PixelSelection : 
            Xs.append(key[0])
            Ys.append(key[1])
            Values.append(value)
        Idxs = (tuple(Xs),tuple(Ys))
        self[Idxs] = Values
    
    ########################################
    ## methode de dessin
    ########################################
    
    def Plot(self,Color="viridis",Vlim=None) :
        if Vlim is None : 
            plt.imshow(self,extent = self.Extent.List(Format=2),cmap=Color)
        else : 
            plt.imshow(self,extent = self.Extent.List(Format=2),cmap=Color,vmin=Vlim[0],vmax=Vlim[1])


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
            self.Extent = GeoExtent(*Extent)
            self.Shape = [self.Extent.Width,self.Extent.Height]
            self.PixelShape = PixelShape
            self.Crs=Crs
            self.Default=Default
        else : 
            DataSource = gdal.Open(Source)
            self.Crs = DataSource.GetProjection()
            ulx, xres, xskew, uly, yskew, yres  = DataSource.GetGeoTransform()
            lrx = ulx + (DataSource.RasterXSize * xres)
            lry = uly + (DataSource.RasterYSize * yres)
            self.Extent = GeoExtent(ulx,lry,lrx,uly)
            self.Shape = [self.Extent.Width,self.Extent.Height]
            self.PixelShape = [xres,-yres]
            self.Shape2 = [DataSource.RasterXSize,DataSource.RasterYSize]
        
        self.Source = Source
        self.MaxSize=MaxSize
        self.__Arrays = defaultdict(lambda : None)
        #2) construire la structures des tuiles
        NbWidth = int(self.Shape2[0]/MaxSize)
        if NbWidth*MaxSize<self.Shape2[0] : 
            NbWidth+=1
        NbHeight = int(self.Shape2[1]/MaxSize)
        if NbHeight*MaxSize<self.Shape2[1] : 
            NbHeight+=1
        TotCol=0
        #iterer sur les colonnes pour delimiter les tuiles
        self.TileIndex = {}
        for i in range(NbWidth) : 
            if TotCol+MaxSize<=self.Shape2[0] : 
                AddCol = MaxSize
            else : 
                AddCol = self.Shape2[0]-TotCol
            StartCol = TotCol
            EndCol = TotCol + AddCol
            MinX = self.Extent.Xmin+self.PixelShape[0]*StartCol
            MaxX = self.Extent.Xmin+self.PixelShape[0]*EndCol
            #iterer sur les lignes pour delimiter les tuiles
            TotRow = 0
            for j in range(NbHeight) : 
                if TotRow+MaxSize<=self.Shape2[1] : 
                    AddRow = MaxSize
                else : 
                    AddRow = self.Shape2[1]-TotRow
                StartRow = TotRow
                EndRow = TotRow+AddRow
                MaxY = self.Extent.Ymax-self.PixelShape[1]*StartRow
                MinY = self.Extent.Ymax-self.PixelShape[1]*EndRow
                self.TileIndex[(i,j)] = {"Extent":GeoExtent(MinX,MinY,MaxX,MaxY),"Pxs":[int(StartCol),int(EndCol),int(StartRow),int(EndRow)]}
                TotRow+=AddRow
            TotCol+=AddCol
        #Creer finalement l'index spatial
        #BBox demande par QIndex : (xmin,ymin,xmax,ymax)
        self.SpIndex = QdIndex(self.Extent.List())
        for key,value in self.TileIndex.items() : 
            self.SpIndex.insert(key,value["Extent"].List())
            
            
    ########################################
    ## methode d'access au rasters
    ########################################
    
    def GetRaster(self,Key,FromSource=True) : 
        """
        Permet d'acceder aux donnees presentes dans le raster source
        """
        if self.__Arrays[Key] is None or FromSource :
            Coordinates = self.TileIndex[Key]
            Window = Coordinates["Pxs"]
            Extent = Coordinates["Extent"]
            DataSource = gdal.Open(self.Source)
            Data = DataSource.ReadAsArray(xoff=Window[0],yoff=Window[2], xsize = Window[1]- Window[0], ysize = Window[3]-Window[2])
            #print(Data)
            Array = SpArray(Data,self.Crs,Extent.List())
            del DataSource
        else : 
            Array = self.__Arrays[Key]
        return Array
    
    def Iterate(self,FromSource=True) : 
        """
        Iterateur permettant de boucler sur toutes les tuils
        """
        for Key in self.TileIndex.keys() : 
            yield self.GetRaster(Key,FromSource=FromSource)
    
    ########################################
    ## methode d'access au pixels
    ########################################
    
    def IntersectedPixels(self,Geom,FromSource=True) : 
        """
        Geom doit etre une geometrie valide de type shapely
        renvoit une selection de pixel triee par tuile :  CoordTile : SpArray
        """
        Extent = Geom.bounds
        Dico = {}
        for key in self.SpIndex.intersect(Extent) : 
            Raster = self.GetRaster(key,FromSource=FromSource)
            if Raster.Extent.Geom.intersects(Geom)==True :
                Pxs = Raster.ScanLine(Geom)
                Dico[key] = Pxs
        return Dico
    
    def SetPixels(self,PixelSelection,FromSource=True) : 
        """
        Permet d'assigner des valeurs aux pixels
        """
        for key,Dico in PixelSelection.items() : 
            Raster = self.GetRaster(key,FromSource=FromSource)
            Raster.SetPixels(Dico)
            
    def Sample(self,Coords,FromSource=True) : 
        """
        Permet d'extraire la valeur d'un pixel a certaine coordonnees
        """
        Extent = [Coords[0],Coords[1],Coords[0],Coords[1]]
        for key in self.SpIndex.intersect(Extent) : 
            Raster = self.GetRaster(key,FromSource=FromSource)
            return Raster.Sample(Coords)
    
    
    ########################################
    ## methodes de calcul de base
    ########################################
    
    def Mean(self) : 
        Sum = 0
        N = 0
        for Raster in self.Iterate() : 
            Sum+=np.sum(Raster)
            N += Raster.shape[0] * Raster.shape[1]
        return Sum/N
    
    def Max(self) : 
        Max = float("-inf")
        for Raster in self.Iterate() :
            m = np.max(Raster)
            if m>Max : 
                Max = m
        return Max
    
    def Min(self) : 
        Min = float("inf")
        for Raster in self.Iterate() :
            m = np.min(Raster)
            if m<Min : 
                Min = m
        return Min
    
    def Std(self) : 
        Mean = self.Mean()
        Sum=0
        N = 0
        for Raster in self.Iterate() : 
            Sum+=np.sum((Raster-Mean)**2)
            N += Raster.shape[0] * Raster.shape[1]
        return Sum/N
        
    
    ########################################
    ## methode de dessin
    ########################################
    
    def PlotGrid(self,LineColor="g",LineWidth=2) : 
        """
        Permet de dessiner les emprises des tuiles du raster multiple
        """
        Lines = []
        for Key,Value in self.TileIndex.items() : 
            Extent = Value["Extent"]
            Lines.append(list(Extent.Geom.boundary.coords))
        lc = LineCollection(Lines,colors=LineColor,linewidths=LineWidth)
        try :
            fig = plt.gcf()
            ax = fig.axes[0]
        except IndexError : 
            fig,ax = plt.subplots()
        ax.add_collection(lc)
        ax.set_xlim(self.Extent.Xmin,self.Extent.Xmax)
        ax.set_ylim(self.Extent.Ymin,self.Extent.Ymax)
        plt.show()
            
        
    def Plot(self,Extent=None,FromSource=True) :
        """
        Permet de dessiner tout ou partie du raster
        """
        if Extent is None : 
            Extent = self.Extent.List()
        
        Keys = self.SpIndex.intersect(Extent)
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
        try :
            fig = plt.gcf()
            ax = fig.axes[0]
        except IndexError : 
            fig,ax = plt.subplots()
        
        ax.set_xlim(Extent[0],Extent[2])
        ax.set_ylim(Extent[1],Extent[3])
        plt.show()
        plt.colorbar()
   
################################
#Tests !


if __name__=="__main__" : 
    
    import time
    
    def timing(f):
        def wrap(*args):
            time1 = time.time()
            ret = f(*args)
            time2 = time.time()
            print('{:s} function took {:.3f} ms'.format(f.__name__, (time2-time1)*1000.0))
    
            return ret
        return wrap
    
    def PlotPoly(Poly,LineColor="r",LineWidth=1) :
        Lines=[list(Poly.boundary.coords)]
        lc = LineCollection(Lines,colors=LineColor,linewidths=LineWidth)
        fig = plt.gcf()
        ax = fig.axes[0]
        ax.add_collection(lc)
        plt.show()
          
    
    @timing
    def test() : 
        Raster1.IntersectedPixels(Poly)
    
    Arr1=np.array([range(10) for e in range(10)])
    Raster1 = SpArray(data=Arr1,Extent=(0,0,100,100),Src="")
    Poly = shapely.wkt.loads("Polygon ((34 20.5, 150 20.5, 75.2 30, 34 50.2, 34 20.5))")
    Pxs1 = Raster1.ScanLine(Poly,Default=0)
    Pxs2 = Raster1.ScanLine(Poly,Default=0,Reverse=True)
    Pxs1.Plot()
#    Raster1.Plot()
#    Xs=[]
#    Ys=[]
#    for Coord in Pxs.keys() : 
#        Pt = Raster1.PxGeom(Coord,Type="Point")
#        Xs.append(Pt.coords[0][0])
#        Ys.append(Pt.coords[0][1])
#    fig = plt.gcf()git 
#    plt.scatter(Xs,Ys,c="r")
#    PlotPoly(Poly)
    
#    Poly = shapely.wkt.loads("Polygon ((654666.14606649207416922 6861235.97581147402524948, 654588.73020110744982958 6861178.4277345510199666, 654517.48020110744982958 6861069.49744608905166388, 654597.63645110744982958 6861114.71379224304109812, 654636.00183572282548994 6861177.05754224304109812, 654692.1797203382011503 6861096.90129224304109812, 654781.92731649207416922 6861160.6152345510199666, 654808.64606649207416922 6861243.511869166046381, 654716.84318187669850886 6861301.74504224304109812, 654663.40568187669850886 6861260.63927301205694675, 654666.14606649207416922 6861235.97581147402524948))")
#    Link = r"L:\TEMP\BatiParis_TestRaster.tif"
#    Raster1 = MultiSpArray(400,Source=Link)
#    Raster1.Plot()
#    Raster1.Plot()
#    Raster1.PlotGrid(LineColor="r")
#    
#    print(Raster1.Sample((654673,6861062)))
#    Key = list(Raster1.SpIndex.intersect([654673,6861062,654673,6861062]))[0]
#    Raster = Raster1.GetRaster(Key)
    
#    test()
    
#    Selection = Raster1.IntersectedPixels(Poly)
#    Xs=[]
#    Ys=[]
#    for key,Dico in Selection.items() : 
#        Raster = Raster1.GetRaster(key)
#        for PxCoords in Dico.keys() : 
#            Pt = Raster.PxGeom(PxCoords,Type="Point")
#            Xs.append(Pt.coords[0][0])
#            Ys.append(Pt.coords[0][1])
#    fig = plt.gcf()
#    plt.scatter(Xs,Ys)
#    PlotPoly(Poly)
    
    