#!/usr/bin/env python

from wasq.AdaptiveSampling import *
import wasq.VoronoiPlot as V

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import argparse


def getopts():
    a = argparse.ArgumentParser()

    a.add_argument('-o', '--outfile', default='cells.png')

    ps = a.add_subparsers()

    pAS = ps.add_parser('AS', help='Adaptive Sampling')
    pAS.add_argument('-c', '--cells',  default='AS/cells')
    pAS.set_defaults(func=main_AS)

    pMD = ps.add_parser('MD', help='Molecular Dynamics')
    pMD.add_argument('-x', '--trajs', nargs='+',  help='Trajectory files')
    pMD.add_argument('-t', '--topology', help='Topology file for loading trajectory')
    pMD.add_argument('-r', '--radius', type=float, default=10)
    pMD.set_defaults(func=main_MD)

    opts = a.parse_args()
    return opts


def main_AS(opts):
    cells = Cells.load_from_dir(opts.cells)
    return cells.C

def main_MD(opts):
    traj = mdtraj.load(opts.trajs, top=opts.topology)
    phipsi = calc_phipsi(traj)
    _, C, _ = PC.poisson_cover(phipsi, opts.radius, metric=dihedral_rmsd)
    return C

def main(opts):
    C = opts.func(opts)
    V.plot_step(C)
    plt.savefig(opts.outfile)


if __name__ == '__main__':
    opts = getopts()
    main(opts)
