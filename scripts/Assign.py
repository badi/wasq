#!/usr/bin/env python

from wasq.util import load_cells_gps

import mdtraj
import numpy as np
import argparse
import re
import copy

def getopts():
    a = argparse.ArgumentParser()

    a.add_argument('-c', '--cells', default='cells.dat', help='AWE-format cells definition')
    a.add_argument('-x', '--trajs', nargs='+')
    a.add_argument('-t', '--topol', default='topol.pdb')
    a.add_argument('-r', '--radius', help='In Angstroms')
    a.add_argument('-n', '--ndx', default='atom.ndx')
    a.add_argument('-a', '--assignments', default='assignments.dat')
    a.add_argument('-p', '--populations', default='populations.dat')

    return a.parse_args()


def main(opts):
    print 'Loading atom indices file for trajectories', opts.ndx
    ndx = np.loadtxt(opts.ndx, dtype=np.int)

    print 'Loading cells from', opts.cells
    cells = mdtraj.load(opts.topol, atom_indices=ndx)
    cells.xyz = load_cells_gps(opts.cells)

    print 'Loading trajectories', ' '.join(opts.trajs)
    traj = mdtraj.load(opts.trajs, top=opts.topol, atom_indices=ndx)

    print 'Assigning to {} cells'.format(len(cells))
    rmsds = -np.ones((len(cells), len(traj)))
    for i in xrange(len(cells)):
        rmsds[i] = mdtraj.rmsd(traj, cells, frame=i)
    rmsds = rmsds.T
    A = -np.ones((len(traj),), dtype=np.int)
    for f in xrange(len(traj)):
        A[f] = rmsds[f].argmin()

    np.savetxt(opts.assignments, A, fmt='%d')

    print 'Computing populations'
    P = np.bincount(A)
    np.savetxt(opts.populations, P, fmt='%d')
        


if __name__ == '__main__':
    main(getopts())
