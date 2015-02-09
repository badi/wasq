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
    a.add_argument('-r', '--reference', default='folded.tpr')
    a.add_argument('-R', '--radius', type=float, default=10)
    a.add_argument('-o', '--outfile', default='fpt.txt')
    a.add_argument('-m', '--mode', choices='a w'.split(), default='w')

    subparsers = a.add_subparsers()

    pAS = subparsers.add_parser('AS', help='Adaptive Sampling')
    pAS.add_argument('-d', '--dir', default='AS')
    pAS.add_argument('-t', '--time', type=float, default=1)
    pAS.set_defaults(func=main_AS)

    pMD = subparsers.add_parser('MD', help='Molecular Dynamics')
    pMD.add_argument('-x', '--trajs', required=True, nargs='+', help='Trajectory files')
    pMD.add_argument('-t', '--topology', default='topol.pdb', help='Topology file for loading trajectories')
    pMD.add_argument('-s', '--timestep', default=1, type=float, help='Time between frames in ps')
    pMD.set_defaults(func=main_MD)


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

def tpr_to_pdb(path, cont=lambda pdb:None):
    tpr = os.path.abspath(path)
    with pxul.os.TmpDir(), disable_gromacs_backups():
        pdb = 'topol.pdb'
        editconf(f=tpr, o=pdb)
        return cont(pdb)

def load_folded_pdb(path):
    traj = mdtraj.load(path)
    phipsi = calc_phipsi(traj)
    return phipsi

def load_folded_tpr(path):
    tpr = os.path.abspath(path)
    return tpr_to_pdb(tpr, cont=load_folded_pdb)

def load_folded(path):
    _, ext = os.path.splitext(path)

    cases = {'.tpr': load_folded_tpr,
             '.pdb': load_folded_pdb}

    load = cases[ext]
    return load(path)

def load_xtcs_tpr(xtcs, tpr):
    traj = tpr_to_pdb(tpr, cont=lambda pdb: load_xtcs_pdb(xtcs, pdb))
    return traj

def load_xtcs_pdb(xtcs, pdb):
    return mdtraj.load(xtcs, top=pdb)

def load_xtcs(xtcs, top):
    _, ext = os.path.splitext(top)

    cases = {'.tpr': load_xtcs_tpr,
             '.pdb': load_xtcs_pdb}

    load = cases[ext]
    return load(xtcs, top)

def discovered_folded(folded, C, R):
    Cnew, _ = PC.online_poisson_cover(folded, R, Cprev=C, metric=dihedral_rmsd)
    return len(Cnew) == len(C)

def main_AS(opts):
    print opts

    cells = Cells.load_from_dir(os.path.join(opts.dir, 'cells'))
    folded = load_folded(opts.reference)
    nwalkers = load_nwalkers(opts.dir)

    found = False
    total_time = 0
    C = cells._cells[0]
    for i in xrange(len(nwalkers)):
        C = np.concatenate((C, cells._cells[i+1]))
        total_time = total_time + nwalkers[i] * opts.time
        print i, len(C), nwalkers[i], total_time,

        found = discovered_folded(folded, C, opts.radius)
        print found
        if found: break

    if not found:
        total_time = float('inf')

    print total_time

    return total_time


def main_MD(opts):
    opts.trajs = map(os.path.abspath, opts.trajs)

    traj = load_xtcs(opts.trajs, opts.topology)
    phipsi = calc_phipsi(traj)
    folded = load_folded(opts.reference)

    for frame in xrange(len(traj)):
        time = frame * opts.timestep
        rmsd = dihedral_rmsd(folded[0].astype(np.float64), phipsi[frame].astype(np.float64))
        found = rmsd < opts.radius

        print frame, '/', len(traj), '\r',

        if found:
            print 'Folded', time, 'ps', 'distance =', rmsd
            break

    print
    return time

def main(opts):
    time = opts.func(opts)

    with open(opts.outfile, opts.mode) as fd:
        fd.write('{}\n'.format(time))

if __name__ == '__main__':
    opts = getopts()
    main(opts)
