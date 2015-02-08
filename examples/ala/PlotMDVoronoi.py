#!/usr/bin/env python

from wasq.AdaptiveSampling import *
import wasq.VoronoiPlot as V

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import argparse


def getopts():
    a = argparse.ArgumentParser()
    a.add_argument('-x', '--trajs', nargs='+',  help='Trajectory files')
    a.add_argument('-t', '--topology', help='Topology file for loading trajectory')
    a.add_argument('-r', '--radius', type=float, default=10)
    a.add_argument('-o', '--outfile', default='cells.png')

    opts = a.parse_args()
    return opts


def main(opts):
    traj = mdtraj.load(opts.trajs, top=opts.topology)
    phipsi = calc_phipsi(traj)
    _, C, _ = PC.poisson_cover(phipsi, opts.radius, metric=dihedral_rmsd)
    V.plot_step(C)
    plt.savefig(opts.outfile)


if __name__ == '__main__':
    opts = getopts()
    main(opts)
