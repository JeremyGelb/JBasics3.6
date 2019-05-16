# -*- coding: utf-8 -*-
"""
Created on Tue May 14 09:48:50 2019

@author: GelbJ
"""

import os
import fnmatch
import re




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
    
    def walkfiles(self,include=None,exclude=None) : 
        # transform glob patterns to regular expressions
        if type(include) is str : 
            include = [include]
        if type(exclude) is str : 
            exclude = [exclude]
        if include is not None :
            includes = r'|'.join([fnmatch.translate(x) for x in include])
        if exclude is not None :
            excludes = r'|'.join([fnmatch.translate(x) for x in exclude]) or r'$.'
        if self.isdir() == False : 
            raise ValueError("I am not a folder : "+self)
        else :
            for root, dirs, files in os.walk(self, topdown=False) : 
                for file in files : 
                    if include is None and exclude is None : 
                        yield Path(root).joinpath(file)
                    elif include is not None and exclude is None : 
                        if re.match(includes,file) :
                            yield Path(root).joinpath(file)
                    elif include is None and exclude is not None : 
                        if not re.match(excludes,file) : 
                            yield Path(root).joinpath(file)
                    else : 
                        if re.match(includes,file) and not re.match(excludes,file) : 
                            yield Path(root).joinpath(file)
    
    def walkdirs(self) : 
        if self.isdir() == False : 
            raise ValueError("I am not a folder : "+self)
        else :
            for root, dirs, files in os.walk(self, topdown=False) : 
                for direct in dirs : 
                    yield Path(root).joinpath(direct)
        
        
    
titi = Path('G:/___Data These')
E = titi.walkfiles("*.shp","*ID1_PA*")
print(next(E))