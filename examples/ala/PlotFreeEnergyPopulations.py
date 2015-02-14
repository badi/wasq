#!/usr/bin/env python

from wasq.AdaptiveSampling import *
import wasq.VoronoiPlot as V
import trax

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import argparse


def getopts():
    a = argparse.ArgumentParser()

    a.add_argument('-o', '--outfile', default='FE.png')
    a.add_argument('-c', '--cells-dir', default='AS/cells')
    a.add_argument('-p', '--populations', default='populations.dat')

    opts = a.parse_args()
    return opts


def main(opts):
    cells = Cells.load_from_dir(opts.cells_dir)
    P = np.loadtxt(opts.populations, dtype=np.int)
    E = -np.log(P)

    print 'Plotting to', opts.outfile

    # dummy figure to show colorbar
    alpha = 0.3
    sE = sorted(set(E))
    energies = dict(zip(sE, range(len(sE))))
    colors = list(reversed('darkblue blue cyan green yellow orange red maroon'.split()))
    colormap = mpl.colors.LinearSegmentedColormap.from_list('mycolors', colors, N=len(sE))
    Z = [[0,0],[0,0]]
    levels = np.array(sE)

    sublevels = levels[range(0, len(levels), len(levels)/20)]
    colorbar = plt.contourf(Z, levels=sublevels, cmap=colormap, alpha=alpha)
    plt.clf()

    def get_colors(labels):
        def color(i):
            j = i % len(P)
            e = E[j]
            k = energies[e]
            c = colormap(k)
            return c
        return color

    V.plot_cells(cells.C, get_colors=get_colors, fill=True, alpha=alpha)
    cb = plt.colorbar(colorbar, format='%0.1f')
    cb.ax.set_ylabel('Free Energy')
    plt.savefig(opts.outfile)



if __name__ == '__main__':
    opts = getopts()
    main(opts)
