#!/usr/bin/env python

from wasq.AdaptiveSampling import *
from wasq.VoronoiPlot import periodic_voronoi, voronoi_finite_polygons_2d

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import argparse


def getopts():
    a = argparse.ArgumentParser()
    a.add_argument('-c', '--cells', default='AS/cells')
    a.add_argument('-o', '--outfile', default='cells.png')

    opts = a.parse_args()
    return opts


def main(opts):
    cells = Cells.load_from_dir(opts.cells)
    I = len(cells._cells)
    cm = mpl.colors.LinearSegmentedColormap.from_list('mycolors', ['blue','red'])

    Z = [[0,0],[0,0]]
    levels = range(0, I-1, 1)
    colorbar = plt.contourf(Z, levels, cmap=cm)
    plt.clf()

    for i in xrange(len(cells._cells) - 1):
        print 'Iteration', i
        phipsi = cells._cells[i]
        # phipsi = np.vstack(cells._cells[:i+1])
        for phi, psi in phipsi:
            n = plt.cm.jet.N / len(cells._cells)
            k = i + i * n
            c = plt.cm.jet(k)
            plt.text(phi, psi, s=str(i), horizontalalignment='center', verticalalignment='center', color=c, fontsize=8)
        plt.axis([-180,180,-180,180])
        plt.xticks(range(-180, 180+60, 60))
        plt.yticks(range(-180, 180+60, 60))
        plt.tight_layout()

    cb = plt.colorbar(colorbar)
    cb.ax.set_ylabel('Iteration')

    # _, v = periodic_voronoi(cells.C)
    # regions, verts = voronoi_finite_polygons_2d(v)
    # for r in regions:
    #     poly = verts[r]
    #     plt.fill(*zip(*poly), linestyle='solid', fill=False)
        

    plt.savefig(opts.outfile)


if __name__ == '__main__':
    opts = getopts()
    main(opts)
