#!/usr/bin/python
import itertools
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches
import matplotlib.collections
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as ar3
import mpl_toolkits.mplot3d.axes3d as ax3
from IPython.core.debugger import Tracer
trace = Tracer()

def generate_cpp_readable_file(filename):
    import shapefile as S
    with open(filename, "w") as f:
        countries = S.Reader("../data/ne_10m_admin_0_countries")
        oceans = S.Reader("../data/ne_10m_ocean")
        cntocean = 0
        for r, c in zip(oceans.iterRecords(), oceans.iterShapes()):
            f.write("%d " % len(c.parts))
            s = np.array(np.array(c.points))
            c.parts.append(len(s))
            for a, b in zip(c.parts, c.parts[1:]):
                cntocean += 1
                f.write("%d " % (b - a))
                f.write(" ".join("%2.20f" % x for x in s[a:b].flatten()))
                f.write("\n")
        print "Ocean Polygons: ", cntocean
        return
        for r, c in zip(countries.iterRecords(), countries.iterShapes()):
            f.write("%d " % len(c.parts))
            s = np.array(np.array(c.points))
            c.parts.append(len(s))
            # shapefiles are in cw order, but
            # google s2 wants ccw order
            for a, b in reversed(zip(c.parts, c.parts[1:])):
                f.write("%d " % (b - a))
                f.write(" ".join("%2.20f" % x for x in s[a:b].flatten()))
                f.write("\n")

def generate_wkt(filename):
    from osgeo import ogr
    # Get the driver
    driver = ogr.GetDriverByName('ESRI Shapefile')
    # Open a shapefile
    countries_shp = "../data/ne_10m_admin_0_countries.shp"
    #countries_shp = "../data/ne_10m_admin_0_boundary_lines_land.shp"
    ocean_shp = "../data/ne_10m_ocean.shp"
    country_ds = driver.Open(countries_shp, 0)
    ocean_ds = driver.Open(ocean_shp, 0)

    L = []
    with open(filename, "w") as f:
        for ds in [country_ds]:
            layer = ds.GetLayer()
            for index in xrange(layer.GetFeatureCount()):
                feature = layer.GetFeature(index)
                try:
                    #print "writing %d, `%s'" % (index, feature['type'])
                    #f.write("%s;%d;%d;" % (feature['adm0_a3_l'], feature['labelrank'], feature['labelrank']))
                    print "writing `%s'" % feature['NAME_LONG']
                    f.write("%s;%d;%d;" % (feature['NAME_LONG'], feature['MAPCOLOR13'], feature['scalerank']))
                except:
                    f.write("OCEAN;")
                geometry = feature.GetGeometryRef()
                L.append(geometry.Clone())
                #geometry for polygon as WKT, inner rings, outer rings etc. 
                f.write(str(geometry))
                f.write('\n')
    return L

def read_triangles(filename):
    import ogr
    L = []
    with open(filename) as f:
        for line in f.readlines():
            geom = ogr.CreateGeometryFromWkt(line)
            L.append(geom)
    return L

def intersections(triangles, countries):
    trace()
    L = []
    for c in countries[:3]:
        cpolys = []
        for t in triangles:
            isec = c.Intersection(t)
            if isec.IsEmpty():
                continue
            cpolys.append(isec)
        L.append(cpolys)
    return L


if __name__ == "__main__":
    #generate_cpp_readable_file("easy.txt")
    #triangles = read_triangles("triangles.txt")
    countries = generate_wkt("easy2.txt")
    #intersect = intersections(triangles, countries)
