#!/usr/bin/env python

"""
Find the conformation furthest from a given reference
"""

from wasq.AdaptiveSampling import *

import mdtraj
import numpy as np
import matplotlib.pyplot as plt

import argparse
import copy

def getopts():
    a = argparse.ArgumentParser()

    a.add_argument('-o', '--outfile')
    a.add_argument('-n', '--ndx')
    a.add_argument('-t', '--topol')

    ps = a.add_subparsers()

    pAS = ps.add_parser('AS', help='Adaptive Sampling')
    pAS.add_argument('-c', '--cells', default='AS/cells')
    pAS.set_defaults(func=main_AS)

    opts = a.parse_args()
    return opts

def find_furthest(reference, conformations, atom_indices):
    rmsds = mdtraj.rmsd(conformations, reference, atom_indices=atom_indices)
    i     = np.argmax(rmsds)
    return rmsds[i], conformations[i]

def main_AS(opts):
    cells = Cells.load_from_dir(opts.cells)
    ndx   = np.loadtxt(opts.ndx, dtype=np.int)
    topol = mdtraj.load(opts.topol)

    trajs = []
    for state in cells.L:
        t = copy.deepcopy(topol)
        t.xyz = state.x
        trajs.append(t)

    traj = trajs[0]
    traj = traj.join(trajs[1:])

    return find_furthest(topol, traj, ndx)


def main(opts):
    rmsd, conf = opts.func(opts)
    print 'Larget rmsd:', rmsd
    conf.save(opts.outfile)

if __name__ == '__main__':
    opts = getopts()
    main(opts)
