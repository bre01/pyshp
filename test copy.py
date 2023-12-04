import shp
import json

""" sf = shapefile.Reader("shapefiles/blockgroups")
geoj=sf.__geo_interface__
with open('a.json','w') as fp:
    json.dump(geoj,fp);
print("done") """


# sf = shp.Reader("shapefiles/test/line.shp")
#sf = shp1.Shp("/Users/bre/source/Open-source-GIS-Course/KunmingProjectedUtm.shp")
sf = shp.Shp("/Users/bre/rgeo-shapefile/test/shapelib_testcases/test.shp")
geoj = sf.__json__
with open("originFinal.json", "w") as fp:
    json.dump(geoj, fp)
print("done")
