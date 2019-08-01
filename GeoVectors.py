# -*- coding: utf-8 -*-
"""
Created on Mon Oct 01 10:48:15 2018

@author: GelbJ
"""

import geopandas as gpd
from osgeo import ogr,osr
import shapely
from shapely import wkt
from shapely import wkb
import numpy as np
import copy
import threading
import queue as Queue
import sqlite3
import os
from GeomTools import GeoExtent
import itertools
import pandas as pd
import geojson

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
    return GeoExtent(MinX,MinY,MaxX, MaxY)

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
# Fonction pour creer un champs OID
#####

def make_oid(self) : 
    self["OID"] = range(len(self))
    New = self.set_index("OID",drop=False)
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
        

gpd.geodataframe.GeoDataFrame.Apply = calculate_field
#############################################


######
## Fonction pour calculer un champs avec une fonction complexe en mode multi threading
######
#
#
#def calculate_field_mc(self,Func,Cores=2,**kwargs) : 
#    """
#    Permet d'utiliser des fonctions complexes pour calculer un nouveau champs en mode multithreading
#    NB : accelere le calcul du champ dans une certaine mesure
#    self referre a geodataframe sur lequel est calcule le champs
#    Feat referre a la feature en cours d'iteration
#    i referre au rang de cette feature dans le geodataframe
#    la fonction doit donc ressembler a : Func(Feat,i,**kwargs)
#    """
#    def Worker(Args) : 
#        Part,Func,DF,kwargs,results = Args
#        for Id in Part :
#            element = DF.iloc[Id,:]
#            results.put((Id,Func(element,Id,**kwargs)))
#            
#    
#    idx = range(self.shape[0])
#    Parts = Chunks(idx,Cores)
#    jobs=[]
#    results = Queue.Queue()
#    for i,Part in enumerate(Parts) : 
#        print("Starting job : "+str(i))
#        p = threading.Thread(target = Worker, args=((Part,Func,self,kwargs,results),))
#        p.start()
#        jobs.append(p)
#    
#    for job in jobs : 
#        job.join()
#    Results = []
#    while not results.empty() :
#        Results.append(results.get())
#    Results.sort(key=lambda x : x[0])
#    return [R[1] for R in Results]
#        
#
#gpd.geodataframe.GeoDataFrame.ApplyMC = calculate_field_mc
##############################################


#####
# Fonction pour calculer un champs avec une fonction complexe en mode multi threading
#####


def calculate_field_mc(self,Func,Cores=2,**kwargs) : 
    """
    Permet d'utiliser des fonctions complexes pour calculer un nouveau champs en mode multithreading
    NB : accelere le calcul du champ dans une certaine mesure
    self referre a geodataframe sur lequel est calcule le champs
    Feat referre a la feature en cours d'iteration
    NB : keep in mind the dataframe is splitted... This will fail is you try to reach features in an other chunk
    i referre au rang de cette feature dans le geodataframe
    la fonction doit donc ressembler a : Func(Feat,i,**kwargs)
    """
    from multiprocessing import Process

    def Worker(Args) : 
        ID,Part,Func,kwargs,results = Args
        Values = Part.apply(Func,axis=1,**kwargs)
        results.put((ID,Values))
            
    
    chunks = np.array_split(self, Cores)
    jobs=[]
    results = Queue.Queue()
    for i,Part in enumerate(chunks) : 
        print("Starting job : "+str(i))
        p = threading.Thread(target = Worker, args=((i,Part,Func,kwargs,results),))
        p.start()
        jobs.append(p)
    
    for job in jobs : 
        job.join()
    Results = []
    print("End of all jobs !")
    while not results.empty() :
        Results.append(results.get())
    Results.sort(key=lambda x : x[0])
    return pd.Series(itertools.chain(Results))
    #return [R[1] for R in Results]
        

gpd.geodataframe.GeoDataFrame.ApplyMC = calculate_field_mc
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
        Index = self.sindex
        indexed = self.iloc[list(Index.intersection(Geom.bounds))]
        #xmin, ymin, xmax, ymax = GetExtent(Geom)
        #indexed = self.cx[xmin:xmax, ymin:ymax]
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
## Enregistrer dans une bd sqlite
###############

def save_sdb(self,File,TableName,GeomField="geometry",Id="OID") : 
    """
    Fonction permettant d'enregistrer le geodataframe dans une BD sqlite
    NB : pour cela on passe comme des brutes par OGR... donc il faut avoir les dependance necessaires
    """
    #0) verifier que la BD existe, sinon la creer
    DB_PATH = File
    Driver=ogr.GetDriverByName('SQLite')
    if os.path.isfile(File)==False : 
        Database = Driver.CreateDataSource(File,options=['SPATIALITE=yes'])
    Source = Driver.Open(File,1)
    GeomType = self.geom_type[0].lower()
        
    if GeomType == "linestring" :
        geom_type = ogr.wkbLineString25D
    elif GeomType == "multilinestring" : 
        geom_type = ogr.wkbMultiLineString
    if GeomType == "point" :
        geom_type = ogr.wkbPoint
    elif GeomType == "multipoint" : 
        geom_type = ogr.wkbMultiPoint
    if GeomType == "polygon" :
        geom_type = ogr.wkbPOLYGON
    elif GeomType == "multipolygon" : 
        geom_type = ogr.wkbMultiPolygon

    #si le layer existe deja, il faut le supprimer
    if Source.GetLayerByName(TableName) is not None :
        Source.DeleteLayer(TableName.lower())
    #recuperer le Src
    EPSG = self.crs["init"].split(':')[1]
    SpatialRef = osr.SpatialReference()
    SpatialRef.ImportFromEPSG(int(EPSG))
    
    SortieLayer = Source.CreateLayer(TableName, srs = SpatialRef,geom_type = geom_type ,options=['FORMAT=SPATIALITE'])
    #creation des champs du layer
    for Name,Value in zip(self.columns,self.iloc[0]) :
        if Name != GeomField :
            Type = type(Value)
            if Type is str : 
                LengthMax =  np.max(list(map(len,self[Name])))
                NewField = ogr.FieldDefn(Name, ogr.OFTString)
                NewField.width=int(LengthMax)
            elif Type is float : 
                NewField = ogr.FieldDefn(Name, ogr.OFTReal)
            elif Type is int : 
                NewField = ogr.FieldDefn(Name, ogr.OFTInteger)
            else : 
                raise ValueError("This Field ("+Name+") type is currently not supported : "+str(Type))
            SortieLayer.CreateField(NewField)
    SortieLayer.StartTransaction()     
    #remplissage du Layer
    FDefn = SortieLayer.GetLayerDefn()
    for row in self.Iterate() :
        Feature = ogr.Feature(FDefn)
        for Field in self.columns :
            if Field != GeomField :
                value = row[Field]
                Feature.SetField(Field,value)
                    
        Feature.SetGeometry(ToOgr(row[GeomField]))
        SortieLayer.CreateFeature(Feature)
    SortieLayer.CommitTransaction()
    del SortieLayer
    Source.Destroy()
        

gpd.geodataframe.GeoDataFrame.SaveSdb = save_sdb



###############
## convertir en objet geojson
###############

def to_geojson(self) :
    """
    Fonction pour convertir notre joli geodataframe en un Geojson object from geojson library
    """
    GeomConverter = {"Point" : geojson.Point,
                     "LineString" : geojson.LineString,
                     "Polygon" : geojson.Polygon,
                     "MultiPoint":geojson.MultiPoint,
                     "MultiLineString":geojson.MultiLineString,
                     "MultiPolygon":geojson.MultiPolygon}
    Features = []
    for Feat in self.iterfeatures() :
        Dico = {Key:Value for Key,Value in Feat.items() if Key!="geometry"}
        Geom = GeomConverter[Feat["geometry"]["type"]](Feat["geometry"]["coordinates"])
        Features.append(geojson.Feature(geometry=Geom,properties=Dico))

    CRS = {"type": "name","properties": {"name": self.crs["init"]}}
    Final = geojson.FeatureCollection(Features,crs=CRS)
    return Final

gpd.geodataframe.GeoDataFrame.ToGeojson= to_geojson


###############
## convertir en objet topojson avec js2py
###############

def to_topojson(self,Name="foo") :
    import js2py
    from path import Path
    Root = str(Path(__file__).parent)
    GeoJson = self.ToGeojson()
    JSFile = open(Root+"/libs/topojson.js","r")
    JSTring = JSFile.read()
    JSTring += "topojson.topology({"+Name+": "+str(GeoJson)+"})"
    try :
        result = js2py.eval_js(JSTring)
    except Exception as e :
        Out = open(str(Path(__file__).parent)+"/logtopojson.txt","w")
        Out.write(JSTring)
        Out.close()
        raise e
    result = result.to_dict()
    result["crs"] = GeoJson["crs"]
    result = str(result).replace("'",'"')
    return result

gpd.geodataframe.GeoDataFrame.ToTopojson= to_topojson

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
        
    To get the projection : 
        DF.crs
    """
    print(Txt)
    
gpd.geodataframe.GeoDataFrame.HelpMe = help_me


###################################################
###################################################
#_______________________________fonctions de lecture de fichiers externes

def ReadSdb(DataBase,LayerName) :
    #ouverture de la couche
    DB = ogr.Open(DataBase)
    Layer = DB.GetLayerByName(LayerName)
    #lecture des champs
    # LayerDef = Layer.GetLayerDefn()
    # Fields = []
    # for i in range(LayerDef.GetGeomFieldCount()) :
    #     Field = LayerDef.GetFieldDefn(i)
    #     Fields.append({
    #             "Name":Field.GetName(),
    #             "Type" : Field.GetTypeName()
    #             })
    #lecture des donnees
    Data = []
    Geoms = []
    for i in range(Layer.GetFeatureCount()) :
        Feat = Layer.GetNextFeature()
        Fields = Feat.items()
        Geoms.append(ToShapely(Feat.GetGeometryRef()))
        Data.append(Fields)
    SpRef = Layer.GetSpatialRef()
    EPSG = SpRef.GetAttrValue("AUTHORITY",1)
    DF = pd.DataFrame.from_dict(Data)
    crs = {'init': 'epsg:'+str(EPSG)}
    GDF = gpd.geodataframe.GeoDataFrame(DF,crs=crs,geometry=Geoms)
    GDF.SpatialFilter = spatial_filter
    return (GDF)



#from multiprocessing import  Pool
#from functools import partial
#import numpy as np
#
#def parallelize(data, func, num_of_processes=8):
#    data_split = np.array_split(data, num_of_processes)
#    pool = Pool(num_of_processes)
#    data = pd.concat(pool.map(func, data_split))
#    pool.close()
#    pool.join()
#    return data
#
#def run_on_subset(func, data_subset):
#    return data_subset.apply(func, axis=1)
#
#def parallelize_on_rows(data, func, num_of_processes=8):
#    return parallelize(data, partial(run_on_subset, func), num_of_processes)
#




if __name__=="__main__" : 
    DF = gpd.read_file(gpd.datasets.get_path("naturalearth_cities"))
    import time
    
    def CostlyFunction(Feat) : 
        for e in range(5000) : 
            Buff = Feat["geometry"].buffer(20)
        return Buff
    
    
    print("duree classique : ")
    Start = time.time()
    V2 = DF.apply(CostlyFunction,axis=1)
    End = time.time()
    print("duration : "+str(End-Start))
    print("")
    print("duree booste 6 : ")
    Start = time.time()
    V1 = DF.ApplyMC(CostlyFunction,Cores=10)
    End = time.time()
    print("duration : "+str(End-Start))    
    
    #DF.SaveSdb("C:/Users/gelbj/OneDrive/Bureau/TEMP/BDTest.sdb",'Cities2')
#    DF = gpd.read_file(FileTest)
#    DF = DF.set_index("objectid",drop=True)
#    DF2 = DF.to_crs(epsg=3857)
#    G = DF2.iloc[1,]["geometry"]
#    DF3 = DF2.SpatialFilter(DF2.iloc[1,]["geometry"].buffer(1500),"intersect")
#    def GetPreviousArea(Feat,i,ThisDF = DF3) : 
#        if i > 0 :
#            PrevFeat = ThisDF.iloc[i-1,:]
#            return PrevFeat["geometry"].area
#        else : 
#            return -999
#    DF3.CalculateFieldMC("PrevArea",GetPreviousArea,Cores=2)
#    DF4 = DF3[["PrevArea","cp_arrondis","geometry"]]
#    DF5 = DF4.AddRows([{"PrevArea":1002,"cp_arrondis":74150,"geometry":G},{"PrevArea":1002,"cp_arrondis":74150,"geometry":G}])
#    DF6 = DF5.MakeOID()
    