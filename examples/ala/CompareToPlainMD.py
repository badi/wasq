#!/usr/bin/env python

from wasq.Cell import Cells
from wasq.AdaptiveSampling import SimulationState, calc_phipsi, dihedral_rmsd
import wasq.PoissonCover as PC
from wasq.VoronoiPlot import plot_step as plot_voronoi

import pxul
import trax
import mdtraj
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import argparse
import os




def getopts():
    p = argparse.ArgumentParser()
    p.add_argument('-r', '--radius', type=float, required=True)
    p.add_argument('-c', '--cells-dir', default='cells')
    p.add_argument('-x', '--trajs', nargs='+')
    p.add_argument('-t', '--top', default=None, help='Topology file for loading the trajectory (if needed)')
    p.add_argument('-d', '--timestep_ps', type=float, default=1)
    p.add_argument('-i', '--iterations', default='AS/iteration')
    p.add_argument('-w', '--walkerlen_ps', default=1, type=float)
    p.add_argument('-o', '--outpath', required=True)
    p.add_argument('-X', '--traxdir', default='.trax')

    opts = p.parse_args()

    return opts


def compare_C(C0, C1):
    assert C0.shape == C1.shape
    D = np.zeros(len(C0))
    for i in xrange(len(C0)):
        D[i] = dihedral_rmsd(C0[i].astype(np.double), C1[i].astype(np.double))

    return D


def contribution(R, new, reference):
    Cnew, _ = PC.online_poisson_cover(new, R, Cprev=reference, metric=dihedral_rmsd)
    contrib = len(Cnew) - len(reference)
    return contrib

def compute_contributions(R, cells, Cmd):
    niter = len(cells._cells)
    iterations = np.arange(niter)
    contrib_AS = np.zeros(niter)
    contrib_MD = np.zeros_like(contrib_AS)

    for i in iterations:
        Cas = np.vstack(cells._cells[:i+1])
        contrib_AS[i] = contribution(R, Cas, Cmd)
        contrib_MD[i] = contribution(R, Cmd, Cas)
        print i, contrib_AS[i], contrib_MD[i]

    return iterations, contrib_AS, contrib_MD


def plot_contributions(iterations, contrib_AS, contrib_MD):
    plt.scatter(iterations, contrib_AS, color='red',  label='AS')
    plt.scatter(iterations, contrib_MD, color='blue', label='MD')
    plt.legend(loc='best')
    plt.xlabel('AS Time Simulated (ns)')
    plt.ylabel('Contribution')


def trax_loggers(traxdir, names='Cmd contrib'):
    pxul.os.ensure_dir(traxdir)
    results = []
    for attr in names.split():
        cpt = os.path.join(traxdir, '{}.cpt'.format(attr))
        log = os.path.join(traxdir, '{}.log'.format(attr))
        logger = trax.SimpleTransactional(checkpoint=cpt, log=log)
        results.append(logger)
    return results

def main(opts):
    trax_Cmd, trax_contrib = trax_loggers(opts.traxdir)

    print 'Loading cells from', opts.cells_dir
    cells = Cells.load_from_dir(opts.cells_dir)

    print 'Loading trajectories', opts.trajs
    traj = mdtraj.load(opts.trajs, top=opts.top)

    print 'Calculating dihedrals'
    phipsi = calc_phipsi(traj)

    print 'Calculating MD covering with radius', opts.radius
    if os.path.exists(trax_Cmd.cpt_path):
        Cmd = trax_Cmd.recover(lambda obj, _: obj)
    else:
        _, Cmd = PC.simple_poisson_cover(phipsi, opts.radius, metric=dihedral_rmsd)
        trax_Cmd.checkpoint(Cmd)

    print 'Computing contributions'
    if os.path.exists(trax_contrib.cpt_path):
        iterations, contrib_AS, contrib_MD = trax_contrib.recover(lambda obj, _: obj)
    else:
        iterations, contrib_AS, contrib_MD = compute_contributions(opts.radius, cells, Cmd)
        trax_contrib.checkpoint((iterations, contrib_AS, contrib_MD))

    print 'Loading nwalkers/iteration from', opts.iterations
    times = np.zeros_like(iterations).astype(np.float)
    for i in xrange(len(times)-1):
        j = i+1
        path = os.path.join(opts.iterations, '{:05d}'.format(i), 'nwalkers.txt')
        nwalkers = open(path).read()
        times[j] = float(nwalkers) * opts.walkerlen_ps
    timescum = np.zeros_like(times)
    for i in xrange(len(times)):
        timescum[i] = times[:i].sum() / 10.**3

    print 'Plotting contributions'
    plot_contributions(timescum, contrib_AS, contrib_MD)
    plt.title('MD = {} $\mu s$'.format(len(traj)*opts.timestep_ps / 10.**6))
    plt.savefig(opts.outpath)


def test():
    class dummy: pass
    opts = dummy()
    opts.cells_dir = '/home/badi/sw/wasq.git/examples/ala/AS/cells'
    opts.trajs = ['/data/ALA_expl-solv/simulations/RUN0/traj-prot_nopbc.xtc']
    opts.top = '/data/ALA_expl-solv/simulations/poisson_cover/cell_0.pdb'
    opts.radius = 10

    return main(opts)

if __name__ == '__main__':
    opts = getopts()
    main(opts)
