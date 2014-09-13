#!/usr/bin/python
import itertools
from collections import defaultdict
import shapefile as S
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches
import matplotlib.collections
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as ar3
import mpl_toolkits.mplot3d.axes3d as ax3
from IPython.core.debugger import Tracer
trace = Tracer()

class SPHTriang:
    def __init__(self, filename):
        self.data = np.loadtxt(filename).reshape(-1, 5)
        self.latlng = self.data[:, :2]
        self.cart = self.data[:, 2:]

    def deg(self):
        return self.latlng * (180 / np.pi)

class FinElem:
    def __init__(self, filename, latlng=True):
        with open(filename) as f:
            L, colors, centroids = [], [], []
            for line in f:
                line = line.split()
                colors.append(float(line[1]))
                #centroids.append(map(float, line[2:4]))
                if latlng:
                    p = np.array(map(float, line[4:])).reshape(-1, 2)
                else:
                    p = np.array(map(float, line[4:])).reshape(-1, 3)
                L.append(p)
        self.centroids = np.array(centroids)
        self.colors = np.array(colors)
        print "Color Range: ", np.ptp(self.colors)
        self.polys = L

def plot_polys_2d(fn):
    fig, ax = plt.subplots(1, 1, figsize=(9, 16))

    #lines = np.loadtxt("src/dual_graph.txt").reshape(-1, 5)
    lines = np.loadtxt(fn).reshape(-1, 5)
    line_colors = lines[:, 0]
    lines = lines[:, 1:].reshape(-1, 2, 2)
    linecol = matplotlib.collections.LineCollection(lines, cmap="jet", lw=3)
    linecol.set_array(line_colors)
    #trace()
    #ax.add_collection(linecol)

    ##fe = FinElem("src/tbits.txt", latlng=True)
    #fe = FinElem("triagrid_poly.txt", latlng=True)
    fe = FinElem("flattened.txt", latlng=True)
    verts = fe.polys
    colors = fe.colors
    def peter(c):
        if c < 0: return 0
        else:     return 0.8
    #colors = map(peter, colors)
    colors = np.array(colors)
    #colors = np.random.uniform(size=len(verts))
    col = matplotlib.collections.PolyCollection(verts,
        array=colors, closed=False, antialiased=2, edgecolor="face", linewidths=1, alpha=1., cmap="Paired")
    #col.set_clim(0, 1)
    ax.add_collection(col)

    ax.set_xlim(-4, 4)
    ax.set_ylim(-3, 1.5)
    ax.set_axis_bgcolor('#909090')
    if False:
        lines = lines.reshape(-1, 4)
        ax.set_xlim(np.min(lines[:,2]), np.max(lines[:,2]))
        ax.set_ylim(np.min(lines[:,3]), np.max(lines[:,3]))
    #ax.xaxis.set_major_locator(plt.NullLocator())
    #ax.yaxis.set_major_locator(plt.NullLocator())
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.show()

def plot_polys_3d():
    fe = FinElem("fe.txt", latlng=False)
    fig = plt.figure(figsize=plt.figaspect(1),)
    ax = fig.add_subplot(111, projection="3d")

    verts = fe.polys
    triCol = ar3.Poly3DCollection(verts)
    z_color = np.random.uniform(size=len(verts))
    triCol.set_array(z_color)
    triCol.set_edgecolor('k')
    triCol.set_linewidth(2.0)
    ax.add_collection(triCol)
    ax.set_zlabel("z")
    ax.set_ylabel("y")
    ax.set_xlabel("x")
    plt.show()

def plot_trias_3d():
    triang = SPHTriang("sphtr.latlng")
    x, y, z = triang.cart.T

    fig = plt.figure(figsize=plt.figaspect(1),)
    ax = fig.add_subplot(111, projection="3d")

    # from (xyz, triangle, point) to (triangle, point, xyz)
    verts = np.array((x, y, z)).reshape(3, -1, 3).transpose(1, 2, 0)
    triCol = ar3.Poly3DCollection(verts)

    z_color = np.random.uniform(size=verts.shape[0])
    triCol.set_array(z_color)
    triCol.set_edgecolor('k')
    triCol.set_linewidth(2.0)
    ax.add_collection(triCol)
    ax.set_zlabel("z")
    ax.set_ylabel("y")
    ax.set_xlabel("x")
    plt.show()

def main():
    countries = S.Reader("../data/ne_10m_ocean")
    boundaries = S.Reader("../data/ne_10m_admin_0_boundary_lines_land")
    popp = S.Reader("../data/ne_10m_populated_places")
    #triang = SPHTriang("sphtr.latlng")
    #fe = FinElem("fe.txt")

    fig = plt.figure()
    ax = fig.add_subplot(111)

    verts, colors = [], []
    for r, c in zip(countries.iterRecords(), countries.iterShapes()):
       s = np.array(np.array(c.points))
       c.parts.append(len(s))
       color = np.random.uniform()
       for a, b in zip(c.parts, c.parts[1:]): 
           verts.append(s[a:b])
           colors.append(color)
    ax.add_collection(matplotlib.collections.PolyCollection(verts,
        array=np.array(colors), edgecolor="black", lw=1))

    #ax.add_collection(matplotlib.collections.PolyCollection(fe.polys,
        #edgecolor="green", facecolor=None, lw=1, alpha=0.2))

    rankdot = defaultdict(list)
    for r, s in zip(popp.iterRecords(), popp.iterShapes()):
       rankdot[r[0]].append((s.points[0][0], s.points[0][1]))
       if r[0] > 2:  # rank
           continue
       ax.text(s.points[0][0], s.points[0][1], r[4], fontsize=7)
    for rank, items in rankdot.iteritems():
        lat, lng = np.array(items).T
        ax.scatter(lat, lng, s=np.exp(-.5*rank)*20, edgecolor=None)

    #deg = triang.deg().reshape(-1, 3, 2)[:, :, ::-1]
    #ax.add_collection(matplotlib.collections.PolyCollection(deg,
    #    edgecolor="black", lw=1, facecolor=None, alpha=.1))

    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    plt.show()

def teststreamplot():
    from mpl_toolkits.basemap import Basemap
    m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
    m.drawmapboundary(fill_color='aqua')
    m.drawcoastlines()
    m.fillcontinents(color='coral',lake_color='aqua')
    m.drawparallels(np.arange(10,90,20))
    m.drawmeridians(np.arange(-180,180,30))
    plt.show()

if __name__ == "__main__":
    plot_polys_2d("triagrid.txt")
    #main()
    #teststreamplot()
