# JBasics3.6
Some tools to work with spatial datas, based on geopandas, gdal, ogr and shapely

Actually the project is not well advanced" I am porting some old stuff from python 2.7 and I decided to go with geopandas instead of continuing with base OGR.
Basicaly, for Vectors dataset, I rely on Geopandas and extend the geodataframe class with some Monkey Patching
For Rasters, I build a class that help to deal with large raster dataset and to find pixels intersected with geometries (not ready actually)
