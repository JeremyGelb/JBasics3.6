# -*- coding: utf-8 -*-
"""
Created on Mon Oct 01 10:48:15 2018

@author: GelbJ
"""

import geopandas as gpd
from osgeo import ogr
import shapely
from shapely import wkt
import numpy as np
import copy
import threading
import queue as Queue

###################################################
###################################################
#_______________________________fonctions de convertion
def ToShapely(Geom) :
    WKT= Geom.ExportToWkt()
    try :
        ShGeom = shapely.wkt.loads(WKT)
        return ShGeom
    except :
        print(WKT)

def ToOgr(Geom) :
    WKT = shapely.wkt.dumps(Geom)
    OgrGeom = ogr.CreateGeometryFromWkt(WKT)
    return OgrGeom


    
###################################################
###################################################
#_______________________________fonctions utilitaires
def Chunks(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))


###################################################
###################################################
#_______________________________fonctions de geometrie vectorielle (shapely based)
    

def GetExtent(Geom) : 
    x,y = zip(*list(Geom.envelope.boundary.coords))
    return min(x), min(y), max(x), max(y)
    
    
###################################################
###################################################

###################################################
###################################################
#_______________________________Monkey patching pour le Geodataframe
    
FileTest = "C:/Users/GelbJ/Desktop/Projets/___Data These/Paris/Paris/GIS Datas/chantiers-perturbants.shp"


######
# Fonction pour calculer l'etendue
#####

def get_extent(self) : 
    Bounds = self.bounds
    MinX = np.min(Bounds["minx"])
    MaxX = np.max(Bounds["maxx"])
    MinY = np.min(Bounds["miny"])
    MaxY = np.max(Bounds["maxy"])
    return (MinX, MaxX, MinY, MaxY)

gpd.geodataframe.GeoDataFrame.Extent = property(get_extent)
#############################################


######
# Fonction pour creer un clone
#####

def clone_me(self) : 
    return copy.deepcopy(self)

gpd.geodataframe.GeoDataFrame.Clone = clone_me
#############################################


#####
# Fonction pour creer une replique vide
#####

def make_me_empty(self) : 
    ToCopy = ["columns","crs","dtypes"]
    New = gpd.GeoDataFrame()
    for Prop in ToCopy : 
        New.__dict__[Prop] = self.__dict__[Prop]
    return New

gpd.geodataframe.GeoDataFrame.MakeOneEmpty = make_me_empty
#############################################


#####
# Fonction pour creer une replique vide
#####

def make_oid(self) : 
    self["OID"] = range(len(self))
    New = self.set_index("OID",drop=True)
    return New
    

gpd.geodataframe.GeoDataFrame.MakeOID = make_oid
#############################################


#####
# Fonction ajouter des lignes dans le DF
#####

def add_rows(self,Rows) :
    """
    NB : rows doit etre une liste de dictionnaire contenant tous les champs necessaires pour
    constituer une ligne, sinon une erreure sera renvoyee
    """
    Gpd2 = gpd.GeoDataFrame(Rows,columns=self.columns)
    return self.append(Gpd2)
        
gpd.geodataframe.GeoDataFrame.AddRows = add_rows


#############################################


#####
# Fonction pour Iterer sur le DF
#####

def iterate(self) : 
    for i in range(self.shape[0]) : 
        yield self.iloc[i,:]
        

gpd.geodataframe.GeoDataFrame.Iterate = iterate
#############################################



#####
# Fonction pour calculer un champs avec une fonction complexe
#####

def calculate_field(self,FieldName,Func,**kwargs) : 
    """
    Permet d'utiliser des fonctions complexes pour calculer un nouveau champs
    self referre a geodataframe sur lequel est calcule le champs
    Feat referre a la feature en cours d'iteration
    i referre au rang de cette feature dans le geodataframe
    la fonction doit donc ressembler a : Func(Feat,i,**kwargs)
    """
    Values = []
    for i in range(self.shape[0]) : 
        Feat = self.iloc[i,:]
        Values.append(Func(Feat,i,**kwargs))
    self[FieldName] = Values
        

gpd.geodataframe.GeoDataFrame.CalculateField = calculate_field
#############################################


#####
# Fonction pour calculer un champs avec une fonction complexe en mode multi threading
#####

def calculate_field_mc(self,FieldName,Func,Cores=2,**kwargs) : 
    """
    Permet d'utiliser des fonctions complexes pour calculer un nouveau champs en mode multithreading
    NB : accelere le calcul du champ dans une certaine mesure
    self referre a geodataframe sur lequel est calcule le champs
    Feat referre a la feature en cours d'iteration
    i referre au rang de cette feature dans le geodataframe
    la fonction doit donc ressembler a : Func(Feat,i,**kwargs)
    """
    def Worker(Args) : 
        Part,Func,DF,kwargs,results = Args
        for Id in Part :
            element = DF.iloc[Id,:]
            results.put((Id,Func(element,Id,**kwargs)))
            
    
    idx = range(self.shape[0])
    Parts = Chunks(idx,Cores)
    jobs=[]
    results = Queue.Queue()
    for i,Part in enumerate(Parts) : 
        print("Starting job : "+str(i))
        p = threading.Thread(target = Worker, args=((Part,Func,self,kwargs,results),))
        p.start()
        jobs.append(p)
    
    for job in jobs : 
        job.join()
    Results = []
    while not results.empty() :
        Results.append(results.get())
    Results.sort(key=lambda x : x[0])
    self[FieldName] = [R[1] for R in Results]
        

gpd.geodataframe.GeoDataFrame.CalculateFieldMC = calculate_field_mc
#############################################



#####
# Fonction d'indexion spatiale plus rapide
#####

def spatial_filter(self,Geom,Methode="intersect") : 
    """
    Geom doit etre de type shapely
    methode acceptee : intersect, dont intersect, inside, touche
    """
    if Methode in ["intersect", "inside", "touche"] :
        #Partie 1 : prefilter avec cx 
        xmin, ymin, xmax, ymax = GetExtent(Geom)
        indexed = self.cx[xmin:xmax, ymin:ymax]
        if Methode == "intersect" : 
            return indexed[indexed.intersects(Geom)]
        elif Methode == "touche" : 
            return indexed[indexed.touches(Geom)]
        elif Methode == "inside" : 
            Sub1 = indexed[indexed.intersects(Geom)]
            return Sub1[Sub1.intersects(Geom.boundary)==False]
    
    elif Methode == "not intersect" : 
        return self[self.intersects(Geom)==False]
    else : 
        raise ValueError("The method must be one of these : "+str(["intersect", "inside", "touche", "not intersect"]))

gpd.geodataframe.GeoDataFrame.SpatialFilter = spatial_filter        
#############################################




###############
## Afficher des fonctions utilitaires
###############
def help_me(self) : 
    Txt = """
    To reproject a dataframe, use this command : 
        world.to_crs({'init': 'epsg:3395'}) # world.to_crs(epsg=3395) would also work
        
    To access a specific Row in a DF, use this : 
        world.loc[id,:]
        
    To change the index column, use this command : 
        DF = DF.set_index("objectid",drop=True)
        
    To access to specifc row number use this command : 
        DF.iloc[id,:]
    """
    print(Txt)
    
gpd.geodataframe.GeoDataFrame.HelpMe = help_me



if __name__=="__main__" : 
    DF = gpd.read_file(FileTest)
    DF = DF.set_index("objectid",drop=True)
    DF2 = DF.to_crs(epsg=3857)
    G = DF2.iloc[1,]["geometry"]
    DF3 = DF2.SpatialFilter(DF2.iloc[1,]["geometry"].buffer(1500),"intersect")
    def GetPreviousArea(Feat,i,ThisDF = DF3) : 
        if i > 0 :
            PrevFeat = ThisDF.iloc[i-1,:]
            return PrevFeat["geometry"].area
        else : 
            return -999
    DF3.CalculateFieldMC("PrevArea",GetPreviousArea,Cores=2)
    DF4 = DF3[["PrevArea","cp_arrondis","geometry"]]
    DF5 = DF4.AddRows([{"PrevArea":1002,"cp_arrondis":74150,"geometry":G},{"PrevArea":1002,"cp_arrondis":74150,"geometry":G}])
    DF6 = DF5.MakeOID()
    