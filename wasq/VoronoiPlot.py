import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import itertools
import collections


def make_periodic_2d(A, box, labels=None):
    """
    Return a periodic tiling of array A

    Parameters:
      A :: Nx2 array : the values to periodicify
      box :: |box| == 2 : the size to shift in the (+/-) directions
      labels :: maybe N array of `a`: label for each point in A
    """

    if labels is not None:
        print labels
        assert len(A) == len(labels), \
          'Size mismatch: number of labels ({}) should equal number of elements ({})'.format(len(A), len(labels))

    L = []
    xs, ys = [], []
    dx, dy = box
    dxs = [-dx, 0, dx]
    dys = [-dy, 0, dy]
    for dx, dy in itertools.product(dxs, dys):
        xs.append(A[:,0] + dx)
        ys.append(A[:,1] + dy)

        if labels is not None:
            L.extend(labels)

    B = np.vstack([np.hstack(xs),
                   np.hstack(ys)]).T

    if labels is None:
        return B, None
    else:
        assert len(B) == len(L), 'Size mismatch |B|={} =/= |L|={}'.format(len(B), len(L))
        return B, L

def periodic_voronoi(X, size=[360, 360], labels=None):
    from scipy.spatial import Voronoi
    X2, L = make_periodic_2d(X, size, labels=labels)
    return X2, L, Voronoi(X2)

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


def plot_cells(C, L=None, get_colors=lambda _:lambda _:None,
               linestyle='solid', fill=False, alpha=None, edgecolor='black',
               marker='.', markercolor='black'):
    """
    C :: NxD array of float: centroids
    L :: maybe N array of a: labels
    """

    _, labels, v = periodic_voronoi(C, labels=L)
    regions, verts = voronoi_finite_polygons_2d(v)
    colors = get_colors(labels)
    for i, r in enumerate(regions):
        color = colors(i)
        poly = verts[r]
        plt.fill(*zip(*poly), linestyle=linestyle, fill=fill, color=color, alpha=alpha, ec=edgecolor)
    plt.scatter(C[:,0], C[:,1], marker=marker, color=markercolor)

    # set axis to dihedral box
    plt.axis([-180,180,-180,180])
    # x/y ticks to from -180 to 180
    for axis in 'x y'.split():
        attr = axis + 'ticks'
        set_tick = getattr(plt, attr)
        set_tick(range(-180,180+60,60))


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
