# -*- coding: utf-8 -*-
"""
Created on Tue May 14 09:48:50 2019

@author: GelbJ
"""

import os





class Path(str) : 
    """
    Classe permettant de gerer les path avec des strings
    sans avoir de mauvaises surprises
    """
    def __new__(cls, content) : 
        return str.__new__(cls,content.replace("//","/").replace("\\","/"))
    
    @property
    def parent(self) : 
        Parts = self.split("/")
        return Path("/".join(Parts[0:-1]))
    
    @property
    def name(self) : 
        return self.split("/")[-1]
    
    def joinpath(self,String) : 
        return Path(self+"/"+String)
    
    def isfile(self) : 
        return os.path.isfile(self)
    
    def isdir(self) : 
        return os.path.isdir(self)
    
    def walkfiles(self) : 
        if self.isdir() == False : 
            raise ValueError("I am not a folder : "+self)
        else :
            for root, dirs, files in os.walk(self, topdown=False) : 
                for file in files : 
                    yield Path(root).joinpath(file)
    
    def walkdirs(self) : 
        if self.isdir() == False : 
            raise ValueError("I am not a folder : "+self)
        else :
            for root, dirs, files in os.walk(self, topdown=False) : 
                for direct in dirs : 
                    yield Path(root).joinpath(direct)
        
        
    