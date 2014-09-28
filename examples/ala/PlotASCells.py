#!/usr/bin/env python

from wasq.AdaptiveSampling import *
import wasq.VoronoiPlot as V

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
    I = cells.chunks
    cm = mpl.colors.LinearSegmentedColormap.from_list('mycolors', ['blue','cyan','green','yellow','orange','red'], N=I)

    # create a dummy figure so that we can show the colorbar
    # otherwise 'no mappeable' error on plt.colorbar

    alpha=0.5
    Z = [[0,0],[0,0]]
    levels = range(0, I-1, 1)
    colorbar = plt.contourf(Z, levels, cmap=cm, alpha=alpha)
    plt.clf()

    C, vor = V.periodic_voronoi(cells.C)
    regions, verts = V.voronoi_finite_polygons_2d(vor)

    for i, r in enumerate(regions):
        c = cells.chunk_of(i % len(cells.C))
        color = mpl.colors.colorConverter.to_rgba(cm(c), alpha=alpha)
        poly = verts[r]
        plt.fill(*zip(*poly), lw=0, ec=color, fc=color, fill=True)

    plt.scatter(cells.C[:,0], cells.C[:,1], marker='.', color='black')

    plt.axis([-180,180,-180,180])
    plt.xticks(range(-180, 180+60, 60))
    plt.yticks(range(-180, 180+60, 60))

    cb = plt.colorbar(colorbar)
    cb.ax.set_ylabel('Iteration')

    plt.xlabel('$\Phi$', fontsize=20)
    plt.ylabel('$\Psi$', fontsize=20)
    plt.tight_layout()
    plt.savefig(opts.outfile)


if __name__ == '__main__':
    opts = getopts()
    main(opts)
