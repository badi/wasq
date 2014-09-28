#!/usr/bin/env python

from wasq.AdaptiveSampling import *
import wasq.PoissonCover as PC
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
    a.add_argument('-a', '--as-dir', default='AS')
    a.add_argument('-t', '--walker-time', type=float, default=1)
    a.add_argument('-r', '--reference', default='folded.tpr')
    a.add_argument('-R', '--radius', type=float, default=10)
    a.add_argument('-o', '--outfile', default='fpt.txt')
    a.add_argument('-m', '--mode', choices='a w'.split(), default='w')

    opts = a.parse_args()
    return opts

def load_nwalkers(as_dir):
    nwalkers = sorted(glob.glob(os.path.join(as_dir, 'iteration', '*', 'nwalkers.txt')))

    def load(path):
        with open(path) as fd:
            n = int(fd.read().strip())
            return n

    nwalkers = map(load, nwalkers)
    return nwalkers

def load_folded_tpr(path):
    tpr = os.path.abspath(path)
    with pxul.os.TmpDir(), disable_gromacs_backups():
        pdb = 'topol.pdb'
        editconf(f=tpr, o=pdb)
        traj = mdtraj.load(pdb)
        phipsi = calc_phipsi(traj)
        return phipsi

def load_folded(path):
    _, ext = os.path.splitext(path)

    cases = {'.tpr': load_folded_tpr}

    load = cases[ext]
    return load(path)


def discovered_folded(folded, C, R):
    Cnew, _ = PC.online_poisson_cover(folded, R, Cprev=C, metric=dihedral_rmsd)
    return len(Cnew) == len(C)

def main(opts):
    cells = Cells.load_from_dir(os.path.join(opts.as_dir, 'cells'))
    folded = load_folded(opts.reference)
    nwalkers = load_nwalkers(opts.as_dir)

    found = False
    total_time = 0
    C = cells._cells[0]
    for i in xrange(len(nwalkers)):
        C = np.concatenate((C, cells._cells[i+1]))
        total_time = total_time + nwalkers[i] * opts.walker_time
        print i, len(C), nwalkers[i], total_time,

        found = discovered_folded(folded, C, opts.radius)
        print found
        if found: break

    if not found:
        total_time = float('inf')

    print total_time

    with open(opts.outfile, opts.mode) as fd:
        fd.write('{}\n'.format(total_time))


if __name__ == '__main__':
    opts = getopts()
    main(opts)
