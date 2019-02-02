# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 19:20:30 2019

@author: gelbj
"""

import shapely
from shapely import geometry,wkt

"""
Ensemble de fonction de geometrie basee sur shapely
"""

##################################################
##fonctions utilitaires liees a l'entendue d'objet
##################################################

def InExtent(Point,Extent) : 
    """
    Extent must be defined as : (minX,minY,MaxX,MaxY)
    """
    X,Y = Point.coords[0]
    if Extent[2]>=X and Extent[0]<=X and Extent[3]>=Y and Extent[1]<=Y : 
        return True
    else :
        return False