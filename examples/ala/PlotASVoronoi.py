#!/usr/bin/env python

from wasq.AdaptiveSampling import *
import wasq.VoronoiPlot as V

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import argparse


def getopts():
    a = argparse.ArgumentParser()
    a.add_argument('-c', '--cells',  default='AS/cells')
    a.add_argument('-o', '--outfile', default='cells.png')

    opts = a.parse_args()
    return opts


def main(opts):
    cells = Cells.load_from_dir(opts.cells)
    V.plot_step(cells.C)
    plt.savefig(opts.outfile)


if __name__ == '__main__':
    opts = getopts()
    main(opts)
