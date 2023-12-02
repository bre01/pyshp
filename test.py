import shp
import json

""" sf = shapefile.Reader("shapefiles/blockgroups")
geoj=sf.__geo_interface__
with open('a.json','w') as fp:
    json.dump(geoj,fp);
print("done") """


#sf = shp.Reader("shapefiles/test/line.shp")
sf = shp.Reader("/Users/bre/source/Open-source-GIS-Course/KunmingProjectedUtm.shp")
geoj=sf.__json__
with open('v1.json','w') as fp:
    json.dump(geoj,fp);
print("done")
