#!/usr/bin/env python

from wasq.AdaptiveSampling import *
import wasq.VoronoiPlot as V

import mdtraj

from scipy.spatial import distance

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import argparse
import os
import glob
import copy


def getopts():
    a = argparse.ArgumentParser()
    a.add_argument('-c', '--cells-dir', default='AS/cells')
    a.add_argument('-t', '--topol', default='topol.pdb')
    a.add_argument('-n', '--ndx', default='index.mndx')
    a.add_argument('-b', '--bins', type=int, default=50)
    a.add_argument('-o', '--outfile', default='cells.png')

    opts = a.parse_args()
    return opts

def plot_rmsd_distribution(cells, topol, atom_indices, bins=50):
    assert type(topol) is mdtraj.Trajectory, 'Expected Trajectory but got {}'.format(type(topotl))

    trajs = []
    for state in cells.L:
        t = copy.deepcopy(topol)
        t.xyz = state.x
        trajs.append(t)

    traj = trajs[0]
    traj = traj.join(trajs[1:])


    rmsds = []
    for frame in xrange(len(traj)):
        r = mdtraj.rmsd(traj, traj, frame=frame, atom_indices=atom_indices)
        rmsds.append(r)
    rmsds = np.vstack(rmsds)

    triu = np.triu_indices(len(rmsds))
    rmsds[triu] = -1
    np.fill_diagonal(rmsds, -1)
    rmsds = rmsds[np.where(rmsds >= 0)]

    plt.hist(rmsds, bins=bins)


if __name__ == '__main__':
    opts = getopts()

    cells = Cells.load_from_dir(opts.cells_dir)
    topol = mdtraj.load(opts.topol)
    idx   = np.loadtxt(opts.ndx, dtype=int)

    plot_rmsd_distribution(cells, topol, idx, bins=opts.bins)

    plt.xlabel('RMSD (nm)')
    plt.ylabel('Count')
    plt.savefig(opts.outfile)
