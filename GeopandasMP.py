import geopandas as gpd
import shapely
from joblib import Parallel, delayed
import numpy as np
import itertools
from datetime import datetime
from tqdm.auto import tqdm

def Function(Row,BuffSize1,BuffSize2,allGDF) : 
    B1 = Row["geometry"].buffer(BuffSize1)
    B2 = Row["geometry"].buffer(BuffSize2)
    Diff = B2.difference(B1)
    
    index = allGDF.sindex
    BBox = Diff.bounds
    Candidates = list(index.intersection(BBox))
    if len (Candidates) == 0 : 
        return 0
    Candidates = allGDF.loc[Candidates]
    OkInter = Candidates[Candidates["geometry"].intersects(Diff)]
    print(str(Row["ORIG_FID"])+" / 13312")
    return len(OkInter["geometry"])

def FunctionNP(Row,BuffSize1,BuffSize2,allGDF) : 
    B1 = Row["geometry"].buffer(BuffSize1)
    B2 = Row["geometry"].buffer(BuffSize2)
    Diff = B2.difference(B1)
    
    index = allGDF.sindex
    BBox = Diff.bounds
    Candidates = list(index.intersection(BBox))
    if len (Candidates) == 0 : 
        return 0
    Candidates = allGDF.loc[Candidates]
    OkInter = Candidates[Candidates["geometry"].intersects(Diff)]
    return len(OkInter["geometry"])

def __apply(gdf,fun,kwargs) :
    return(gdf.apply(fun,axis=1,**kwargs))

def Parallelize(gdf,fun,ncores,**kwargs) : 
    chunks = np.array_split(gdf, ncores)
    AllValues = Parallel(n_jobs=ncores)(delayed(__apply)(chunk,fun,kwargs) for chunk in chunks)
    return list(itertools.chain.from_iterable(AllValues))
    
    
    
Data = gpd.read_file(r"H:\Cours\Analyse spatiale\Donnees et exemples\Seance07-09\Datas\Ilots_avec_population_Points.gpkg")


### Test de la fonction en mode uniprocess
tqdm.pandas(desc="my bar!")
T1 = datetime.now()
Polys = Data.apply(Function,axis=1,BuffSize1=500,BuffSize2=1000,allGDF=Data)
T2 = datetime.now()
Diff = (T2-T1).total_seconds()
print("Seconds ellapsed for the uniprocessing: "+str(round(Diff)))

### Test de la fonction en mode multiprocess 2 coeurs
T1 = datetime.now()
Polys = Parallelize(Data,FunctionNP,2,BuffSize1=500,BuffSize2=1000,allGDF=Data)
T2 = datetime.now()
Diff = (T2-T1).total_seconds()
print("Seconds ellapsed for the multiprocessing (2cores): "+str(round(Diff)))

### Test de la fonction en mode multiprocess 8 coeurs
T1 = datetime.now()
Polys = Parallelize(Data,FunctionNP,8,BuffSize1=500,BuffSize2=1000,allGDF=Data)
T2 = datetime.now()
Diff = (T2-T1).total_seconds()
print("Seconds ellapsed for the multiprocessing (8cores): "+str(round(Diff)))
