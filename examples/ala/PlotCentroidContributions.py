#!/usr/bin/env python

from wasq.AdaptiveSampling import *
import wasq.PoissonCover as PC
import wasq.VoronoiPlot as V

from scipy.spatial import Voronoi
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import cPickle as pickle

import argparse


def getopts():
    a = argparse.ArgumentParser()
    a.add_argument('-r', '--radius', type=float, default=10)
    a.add_argument('-c', '--cells', default='AS/cells')
    a.add_argument('-o', '--outfile', default='cells.png')
    a.add_argument('-m', '--md-cover', required=True)

    opts = a.parse_args()
    return opts


def main(opts):
    cells = Cells.load_from_dir(opts.cells)
    Cmd = pickle.load(open(opts.md_cover, 'rb'))

    Cnew, _ = PC.online_poisson_cover(cells.C, opts.radius, Cprev=Cmd)

    import pdb;pdb.set_trace()



if __name__ == '__main__':
    opts = getopts()
    main(opts)
