from IPython.core.debugger import Tracer
trace = Tracer()
import ogr
# Get the driver
driver = ogr.GetDriverByName('ESRI Shapefile')
# Open a shapefile
shapefileName = "../../data/ne_10m_admin_0_countries.shp"
dataset = driver.Open(shapefileName, 0)

layer = dataset.GetLayer()
#trace()
for index in xrange(layer.GetFeatureCount()):
    feature = layer.GetFeature(index)
    print "%s;" % feature['NAME_LONG'],
    geometry = feature.GetGeometryRef()
    #geometry for polygon as WKT, inner rings, outer rings etc. 
    print geometry
