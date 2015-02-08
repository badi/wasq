#!/usr/bin/env python

from wasq.AdaptiveSampling import *
import wasq.VoronoiPlot as V
import trax

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import argparse


def getopts():
    a = argparse.ArgumentParser()

    a.add_argument('-o', '--outfile', default='FE.png')
    a.add_argument('-T', '--traxdir', default='.trax')

    ps = a.add_subparsers()

    pAS = ps.add_parser('AS', help='Adaptive Sampling')
    pAS.add_argument('-c', '--cells',  default='AS/cells')
    pAS.add_argument('-n', '--numsims', type=int, help='Number of simulations to run from each cell')
    pAS.add_argument('-t', '--time', type=float, help='Time (in ps) for each simulation')
    pAS.add_argument('-r', '--')
    pAS.set_defaults(func=main_AS)

    pMD = ps.add_parser('MD', help='Molecular Dynamics')
    pMD.add_argument('-x', '--trajs', nargs='+',  help='Trajectory files')
    pMD.add_argument('-t', '--topology', help='Topology file for loading trajectory')
    pMD.add_argument('-r', '--radius', type=float, default=10)
    pMD.set_defaults(func=main_MD)

    opts = a.parse_args()
    return opts


def main_MD(opts):
    tx = trax.Trax(traxdir=os.path.join(opts.traxdir, 'MD'))

    print 'Loading trajectories', opts.trajs
    traj = tx.recover('traj', create=lambda:mdtraj.load(opts.trajs, top=opts.topology))

    print 'Calculating Phi/Psi angles'
    phipsi = tx.recover('phipsi', create=lambda:calc_phipsi(traj))

    print 'Finding cells (radius = {})'.format(opts.radius)
    cover = tx.recover('cover', create=lambda:PC.poisson_cover(phipsi, opts.radius, metric=dihedral_rmsd))
    _, C, _ = cover

    print 'Assigning {} conformations to {} cells'.format(len(traj), len(C))
    A = tx.recover('A', create=lambda:PC.assign(phipsi, C, metric=dihedral_rmsd))
    assert len(A) == len(phipsi), 'Mismatched length: A should be {} but is {}'.format(len(phipsi), len(A))

    print 'Calculating popupations'
    P = np.zeros(len(C))
    for a in A:
        P[a] += 1

    print 'Calculating Free Energy'
    E = -np.log(P)

    print 'Plotting to', opts.outfile

    # dummy figure to show colorbar
    alpha = 0.75
    sE = list(sorted(set(E)))
    energies = dict(zip(sE, range(len(sE))))
    colors = list(reversed('blue cyan green yellow orange red maroon'.split()))
    colormap = mpl.colors.LinearSegmentedColormap.from_list('mycolors', colors, N = len(sE))
    Z = [[0,0],[0,0]]
    levels = sE
    colorbar = plt.contourf(Z, levels, cmap=colormap, alpha=alpha)
    plt.clf()

    def get_colors(labels):
        def color(i):
            j = i % len(E)
            e = E[j]
            c = energies[e]
            return colormap(c)
        return color

    V.plot_cells(C, get_colors=get_colors, fill=True, alpha=alpha)
    cb = plt.colorbar(colorbar, format='%0.1f')
    cb.ax.set_ylabel('Free Energy')
    plt.savefig(opts.outfile)



def create_MD_Tasks(simstate):

def main_AS(opts):
    tx = trax.Trax(traxdir=os.path.join(opts.traxdir, 'AS'))

    print 'Loading cells', opts.cells
    cells = Cells.load_from_dir(opts.cells)

    

def main(opts):
    opts.func(opts)

if __name__ == '__main__':
    opts = getopts()
    main(opts)
