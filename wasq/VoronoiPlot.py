
import matplotlib.pyplot as plt

import numpy as np


def periodic_voronoi(X, size=360):
    from scipy.spatial import Voronoi
    from itertools import product
    ps, hs = [], []
    directions = [-size,0,size]
    for dxnX, dxnY in product(directions, directions):
        ps.append(X[:,0]+dxnX)
        hs.append(X[:,1]+dxnY)
    X2 = np.vstack([np.hstack(ps),np.hstack(hs)]).T
    return X2, Voronoi(X2)

def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    """
    from collections import defaultdict

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = defaultdict(list)
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges[p1].append((p2,v1,v2))
        all_ridges[p2].append((p1,v1,v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)


def plot_step(C):
    _,dim = C.shape
    subplots = dim / 2
    plt.figure(figsize=plt.figaspect(1/float(subplots)))
    for j,i in enumerate(xrange(0, dim-1, 2)):
        plt.subplot(1,subplots,j)
        phipsi = C[:,[i,i+1]]
        _, v = periodic_voronoi(phipsi)
        regions, verts = voronoi_finite_polygons_2d(v)
        for r in regions:
            poly = verts[r]
            plt.fill(*zip(*poly), linestyle='solid', fill=False)
        plt.scatter(phipsi[:,0], phipsi[:,1], marker='.')
        plt.axis([-180,180,-180,180])
        plt.xticks(range(-180, 180+60, 60))
        plt.yticks(range(-180, 180+60, 60))
    plt.tight_layout()


if __name__ == '__main__':
    import sys

    Cpath = sys.argv[1]
    Opath = sys.argv[2]

    C = np.loadtxt(Cpath)
    # numpy.savetxt flattens a (1,D) shaped array into a (D,) shaped
    # array, so fix it on the reload
    if len(C.shape) == 1:
        C = C.reshape((1,len(C)))
    plot_step(C)
    plt.savefig(Opath)
