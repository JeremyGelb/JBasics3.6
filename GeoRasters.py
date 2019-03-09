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
import math
from libs.QuadTree import Index as QdIndex
import GeomTools as GT
from GeomTools import GeoExtent
from GeoVectors import gpd

##############################################################
## Definition des erreurs que l'on peut rencontrer
##############################################################
class ExtentError(Exception):
    """
    Raise me when we ask something out of the bound of an object
    """
    pass


class RasterIntegrityError(Exception):
    """
    Raise me when operations are performed on rasters that don't share there
        -pixelshape
        -crs
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
        
        
    def Clone(self) : 
        return copy.deepcopy(self)
        
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
    
    
    def SetValues(self,Target,TempValue=-999) : 
        """
        Permet d'assigner les valeurs du raster provenant d'un autre raster dont les pixels se superposent parfaitement
        """
        Inter = self.Extent.Geom.intersection(Target.Extent.Geom)
        
        OldPixels = self.ScanLine(Inter,TempValue,Reverse=True)
        NewPixels = Target.ScanLine(Inter,TempValue,Reverse=True)
        Idx1 = np.where(OldPixels==TempValue)
        Idx2 = np.where(NewPixels==TempValue)
        print(Idx1)
        print(Idx2)
        OldPixels[Idx1[0],Idx1[1]] = Target[Idx2[0],Idx2[1]]
        return OldPixels
        
    
    def ArrayCoords(self) : 
        """
        Permet de recuperer deux Array, contenant les coordonnes X et Y des centres des pixels de cet Array
        """
        
        def GetX(i,j) : 
            return self.Extent.Xmin+((j)*self.PixelShape[0])+0.5*self.PixelShape[0]
        def GetY(i,j) : 
            return self.Extent.Ymax-((i)*self.PixelShape[1])-0.5*self.PixelShape[1]

        X = np.fromfunction(GetX,self.shape,dtype=float)
        Y = np.fromfunction(GetY,self.shape,dtype=float)
        return SpArray(X,self.Src,self.Extent.List()),SpArray(Y,self.Src,self.Extent.List())
        
    ########################################
    ## Methode d'analyses spatiales
    ########################################
    def KernelDensity(self,Points,Radius,Function) : 
        """
        Calcule la densite kernel
        """
        X,Y = self.ArrayCoords()
        Kernels = self.Clone()
        Kernels.fill(0)
        for Pt in Points :
            Distances = np.sqrt((X-Pt.x)**2+(Y-Pt.y)**2)
            Kernel = Function(Distances,Radius)
            Kernels = Kernels+Kernel
        return Kernels
    
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
            NB : Maxsize is th emaximum possible size for a tile in pixels
        
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
            self.Shape2 = [self.Extent.Width/PixelShape[0],self.Extent.Height/PixelShape[1]]
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
        self._MemArray = defaultdict(lambda : None)
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
        self.TileShape = (NbWidth,NbHeight)
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
            
    
    def EmptyCopy(self,Default = 0) : 
        """
        retourne un multi array d'un format identique, sans source, avec une valeur par defaut
        """
        NewArray = MultiSpArray(self.MaxSize,Source=None,Extent=self.Extent.List(),PixelShape=self.PixelShape,Crs=self.Crs,Default=Default)
        return NewArray
       
    ########################################
    ## methode d'access au rasters
    ########################################
    
    def Merge(self,FromSource=True) : 
        """
        Permet de fusionner les differentes tuiles pour ne former qu'un spArray
        NB : attention a ne pas depasser la memoire disponible
        """
        Cols = []
        for i in range(self.TileShape[0]) : 
            Col = []
            for j in range(self.TileShape[1]) :
                Raster = self.GetRaster((i,j),FromSource)
                Col.append(Raster)
            Cols.append(np.concatenate(Col,axis=0))
        NewArray = np.concatenate(Cols,axis=1)
        NewSpArray = SpArray(NewArray,self.Crs,self.Extent.List())
        return NewSpArray
        
    
    def GetRaster(self,Key,FromSource=True) : 
        """
        Permet d'acceder aux donnees presentes dans le raster source
        """
        if self._MemArray[Key] is None and self.Source is None : 
            Coordinates = self.TileIndex[Key]
            Window = Coordinates["Pxs"]
            Extent = Coordinates["Extent"]
            xsize = Window[1]- Window[0]
            ysize = Window[3]-Window[2]
            Arr = np.full((ysize,xsize),self.Default)
            SpArr = SpArray(Arr,self.Crs,Extent.List())
            self._MemArray[Key] = SpArr
            return SpArr
        
        if self.Source is None and FromSource : 
            raise ValueError("Impossible to use fromsource here because there is no source...")
            
        if self._MemArray[Key] is None or FromSource :
            Coordinates = self.TileIndex[Key]
            Window = Coordinates["Pxs"]
            Extent = Coordinates["Extent"]
            DataSource = gdal.Open(self.Source)
            Data = DataSource.ReadAsArray(xoff=Window[0],yoff=Window[2], xsize = Window[1]- Window[0], ysize = Window[3]-Window[2])
            #print(Data)
            Array = SpArray(Data,self.Crs,Extent.List())
            del DataSource
        else : 
            Array = self._MemArray[Key]
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
    
    def IntersectedPixels(self,Geom,Default=np.nan,Reverse=False,FromSource=True) : 
        """
        Geom doit etre une geometrie valide de type shapely
        renvoit une selection de pixel triee par tuile :  CoordTile : SpArray
        """
        Extent = Geom.bounds
        Dico = {}
        Extents = []
        Is=[]
        Js=[]
        #recuperation des raster intersectes
        for key in self.SpIndex.intersect(Extent) : 
            Is.append(key[0])
            Js.append(key[1])
            Raster = self.GetRaster(key,FromSource=FromSource)
            if Raster.Extent.Geom.intersects(Geom)==True :
                Pxs = Raster.ScanLine(Geom,Default,Reverse)
                Dico[key] = Pxs
                print(Pxs.Extent.List())
                Extents.append(Pxs.Extent)
        #calcul de l'etendue globale de l'intersection
        print("step2")
        MinI = min(Is)
        MinJ = min(Js)
        Ex1 = Extents.pop(0)
        TotExt = Ex1.Merge(Extents)
        print(TotExt.List())
        NewArray = MultiSpArray(self.MaxSize,Source=None,Extent=TotExt.List(),PixelShape=self.PixelShape,Crs=self.Crs,Default=Default)
        for key,Raster in Dico.items() : 
            Newkey = (key[0]-MinI,key[1]-MinJ)
            NewArray._MemArray[Newkey]=Raster
        
        return NewArray
    
#    def SetPixels(self,PixelSelection,FromSource=True) : 
#        """
#        Permet d'assigner des valeurs aux pixels
#        """
#        for key,Dico in PixelSelection.items() : 
#            Raster = self.GetRaster(key,FromSource=FromSource)
#            Raster.SetPixels(Dico)
            
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
        return np.sqrt(Sum/N)
        
    ########################################
    ## methode d'analyse spatiale
    ########################################
    ########################################
    ## Methode d'analyses spatiales
    ########################################
    def KernelDensity(self,Points,Radius,Function) : 
        """
        Calcule la densite kernel
        """
        #Etape1 : faire un arbre pour acceder plus vite aux points
            
        Copy = self.EmptyCopy(0)
        for key in self.TileIndex.keys() : 
            Raster = Copy.GetRaster(key)
            OkPts = [Pt for Pt in Points if Pt.distance(Raster.Extent.Geom)<Radius]
            Kernel = Raster.KernelDensity(OkPts,Radius,Function)
            Copy._MemArray[key] = Kernel
        return Copy
    
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
            print(Key)
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
        
    ########################################
    ## methode d'export
    ########################################
    def Save(self,File,FromSource=True,Default=None) : 
        """
        Export the full raster to a file
        """
        if Default is None : 
            Default = self.Default
        Format = File.split(".")[1]
        if Format == "tif" : 
            self._TiffSave(File,Default,FromSource)
        else : 
            raise ValueError("The supported format are only : .tif")
        

    def _TiffSave(self,File,Default,FromSource=True) : 
        """
        Methode pour exporter un gros raster avec Gdal en ecrivant chunk par chunk
        """
        #Creation de la source de donnees
        RType = self.GetRaster((0,0)).dtype
        if RType is np.dtype('float64') : 
            Type = gdal.GDT_Float64
        elif RType is np.dtype('float32') : 
            Type = gdal.GDT_Float32
        elif RType is np.dtype('int32') :
            Type = gdal.GDT_Integer32
        elif RType is np.dtype('int64') :
            Type = gdal.GDT_Integer64
        else : 
            raise ValueError("The type od the raster is not supported : "+str(RType))
        
        driver = gdal.GetDriverByName('GTiff')
        DataSource = driver.Create(
            File,
            self.Shape2[0],  # xsize
            self.Shape2[1],  # ysize
            1,    # number of bands
            Type,  # The datatype of the pixel values
            options=[  # Format-specific creation options.
                'TILED=YES',
                'BIGTIFF=IF_SAFER',
                'BLOCKXSIZE=256',  # must be a power of 2
                'BLOCKYSIZE=256'   # also power of 2, need not match BLOCKXSIZE
            ])
        DataSource.SetGeoTransform((self.Extent.List()[0],self.PixelShape[0],0,self.Extent.List()[3],0,-self.PixelShape[1]))
        DataSource.SetProjection(self.Crs)        
        #DataSource.SetNoDataValue(-1)
        #DataSource.fill(-1)
        Band = DataSource.GetRasterBand(1)
        Band.SetNoDataValue(Default)
        #print("filling with default (may take a while)")
        #Band.fill(Default)
        #ecriture dans la source de donnees
        #Window=((WindowMaxY,WindowMinY),(WindowMaxX,WindowMinX))
        for key,window in self.TileIndex.items() :
            Array = self.GetRaster(key,FromSource)
            #print(Array)
            Window = window["Pxs"]
            #"Pxs":[int(StartCol),int(EndCol),int(StartRow),int(EndRow)]
            #print(Array.shape)
            #print((Window[0],Window[3]))
            Band.WriteArray(Array,xoff=Window[0],yoff=Window[2])
        Band.FlushCache()
        Band=None
        DataSource=None

##############################################################
## Fonctions de generations
##############################################################
        
###########################
##Fonction pour combiner plusieurs Array en 1 MultiArray    
def CombineArray(Arrays, Default) : 
    """
    This function will combine single Arrays in a MultiArray that will be irregular
    The missing part will be fill with the default value
    NB : all the rasters must have the same projection and the same pixel size
    if rasters are overlapping, the last value will be keeped
    """
    #### Recuperation de l'etendue totale


###########################
##Fonction pour generer une densite de type Kernel
def InverseLogDistance(Dists,BandWidth) : 
    Values =  1.0/np.log(Dists)
    Values[Dists>BandWidth] = 0
    return Values

def GaussianDistance(Dists,BandWidth) : 
    Z = Dists / float(BandWidth)
    Values = (2*math.pi)**(0.5)*np.exp(-1*(Z)**2)
    Values[Dists>BandWidth] = 0
    return Values

def QuadraticDistance(Dists,BandWidth) : 
    Z = Dists / float(BandWidth)
    Values =  (15.0/16.0)*(1-Z**2)**2
    Values[Dists>BandWidth] = 0
    return Values
    
        
    
    
################################
#Tests !


if __name__=="__main__" : 
    
#    import time
#    
#    def timing(f):
#        def wrap(*args):
#            time1 = time.time()
#            ret = f(*args)
#            time2 = time.time()
#            print('{:s} function took {:.3f} ms'.format(f.__name__, (time2-time1)*1000.0))
#    
#            return ret
#        return wrap
#    
#    def PlotPoly(Poly,LineColor="r",LineWidth=1) :
#        Lines=[list(Poly.boundary.coords)]
#        lc = LineCollection(Lines,colors=LineColor,linewidths=LineWidth)
#        fig = plt.gcf()
#        ax = fig.axes[0]
#        ax.add_collection(lc)
#        plt.show()
#          
#    
#    @timing
#    def test() : 
#        Raster1.IntersectedPixels(Poly)
#        
#        
#    Arr1 = np.array([range(10) for e in range(10)])
#    Raster1 = SpArray(data=Arr1,Extent=(0,0,100,100),Src="")
#    
#    Arr2 = np.array([range(7) for e in range(7)])
#    Raster2 = SpArray(data=Arr2,Extent=(50,50,120,120),Src="")
#    Raster1.Plot()
#    Raster1.Extent.Plot()
#    Raster2.Plot()
#    Raster2.Extent.Plot()
#    Raster3 = Raster1.SetValues(Raster2)
#    Raster3.Plot()
#    Poly = shapely.wkt.loads("Polygon ((34 20.5, 150 20.5, 75.2 30, 34 50.2, 34 20.5))")
#    Pxs1 = Raster1.ScanLine(Poly,Default=0)
#    Pxs2 = Raster1.ScanLine(Poly,Default=0,Reverse=True)
#    Pxs1.Plot()
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
    
    Poly = shapely.wkt.loads("Polygon ((654666.14606649207416922 6861235.97581147402524948, 654588.73020110744982958 6861178.4277345510199666, 654517.48020110744982958 6861069.49744608905166388, 654597.63645110744982958 6861114.71379224304109812, 654636.00183572282548994 6861177.05754224304109812, 654692.1797203382011503 6861096.90129224304109812, 654781.92731649207416922 6861160.6152345510199666, 654808.64606649207416922 6861243.511869166046381, 654716.84318187669850886 6861301.74504224304109812, 654663.40568187669850886 6861260.63927301205694675, 654666.14606649207416922 6861235.97581147402524948))")
    Link = r"I:\TEMP\BatiParis_TestRaster.tif"
    Raster1 = MultiSpArray(100,Source=Link)
    Raster1.Save(r"I:\TEMP\BatiParis_TestRaster_copied.tif",FromSource=True,Default=-999)
#    Merged = Raster1.Merge()
#    Raster1.Plot()
#    Merged.Plot()
#    X,Y = Merged.ArrayCoords()
#    X.Plot()
#    Y.Plot()
#    Center = Poly.centroid
#    Centercoords = (Center.x,Center.y)
#    Distance = np.sqrt((X-Centercoords[0])**2+(Y-Centercoords[1])**2)
#    
#    #K1 = KernelDensity(Merged,[Poly.centroid],100,GaussianDistance)
#    #K2 = KernelDensity(Merged,[Poly.centroid],100,InverseLogDistance)
#    K3 = Merged.KernelDensity([Poly.centroid],100,QuadraticDistance)
#    K4 = Raster1.KernelDensity([Poly.centroid],100,QuadraticDistance)
##    Raster1.Plot()
##    Raster1.PlotGrid(LineColor="r")
##    
###    print(Raster1.Sample((654673,6861062)))
###    Key = list(Raster1.SpIndex.intersect([654673,6861062,654673,6861062]))[0]
###    Raster = Raster1.GetRaster(Key)
##
##    
#    Selection = Raster1.IntersectedPixels(Poly,Default=0,Reverse=False)
#    Selection.Plot(FromSource=False)
#    for key,Raster in Selection.items() : 
#        Raster.Plot()
    #Selection.Mean()
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
    
    