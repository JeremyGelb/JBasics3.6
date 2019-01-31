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
# Fonction d'indexion spatiale plus rapide
#####

def spatial_filter(self,Geom,Methode) : 
    """
    Geom doit etre de type shapely
    methode acceptee : intersect, dont intersect, inside, touche
    """
    if Methode in ["intersect", "inside", "touche"] :
        #Partie 1 : prefilter avec cx 
        xmin, ymin, xmax, ymax = GetExtent(Geom)
        indexed = self.cx[xmin:xmax, ymin:ymax]
        if Methode == "intersect" : 
            return indexed[indexed.intersect(Geom)]
        elif Methode == "touche" : 
            return indexed[indexed.touches(Geom)]
        elif Methode == "inside" : 
            Sub1 = indexed[indexed.intersect(Geom)]
            return Sub1[Sub1.intersect(Geom.boundary)==False]
    
    elif Methode == "dont intersect" : 
        return self[self.intersect(Geom)==False]
        


DF = gpd.read_file(FileTest)
DF2 = DF.Clone()