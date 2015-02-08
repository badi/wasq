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
    a.add_argument('-c', '--cells', default='AS/cells')
    a.add_argument('-x', '--trajectories', nargs='+')
    a.add_argument('-t', '--topology')
    a.add_argument('-r', '--radius', type=float)

    return a.parse_args()


def main(opts):
    tx = trax.Trax(traxdir=opts.traxdir)

    print 'Loading trajectories', opts.trajectories
    # traj = tx.recover('traj', create=lambda:mdtraj.load(opts.trajectories, top=opts.topology),
    #                   dumper=lambda obj, path: obj.save_hdf5(path),
    #                   loader=mdtraj.load_hdf5)
    traj = mdtraj.load(opts.trajectories, top=opts.topology)

    print 'Calculating Phi/Psi'
    phipsi = tx.recover('phipsi', create=lambda:calc_phipsi(traj))

    print 'Finding MD cells (radius = {})'.format(opts.radius)
    cover = tx.recover('cover', create=lambda:PC.poisson_cover(phipsi, opts.radius, metric=dihedral_rmsd))
    _, Cmd, _ = cover

    print 'Loading AS cells from', opts.cells
    cells = tx.recover('cells', create=lambda:Cells.load_from_dir(opts.cells))
    Cas   = cells.C

    print 'Finding contribution of MD'
    contrib_as = tx.recover('contrib_as', create=lambda:PC.online_poisson_cover(Cas, opts.radius, Cprev=Cmd, metric=dihedral_rmsd))
    Cnew, _ = contrib_as

    print 'Plotting to', opts.outfile
    origin = np.zeros(len(Cnew))
    origin[len(Cmd):] = 1

    perc_AS = 100 * len(origin[np.where(origin==1)]) / float(len(origin))
    perc_MD = 100 * len(origin[np.where(origin==0)]) / float(len(origin))

    print '{:0.2f}% AS'.format(perc_AS)
    print '{:0.2f}% MD'.format(perc_MD)
    print origin

    def get_colors(labels):
        def color(i):
            j = i % len(origin)
            c = origin[j]
            if c == 0:
                col = mpl.colors.colorConverter.to_rgba('white', alpha=0)
            elif c == 1:
                col = mpl.colors.colorConverter.to_rgba('blue', alpha=0.5)
            else: raise ValueError
            return col
        return color

    plt.figure(figsize=plt.figaspect(1/3.0))

    ax = plt.subplot(1,3,1)
    V.plot_cells(Cmd)
    plt.ylabel('$\Psi$')
    plt.xlabel('$\Phi$')
    plt.title('MD')

    plt.subplot(1,3,2, sharey=ax)
    V.plot_cells(Cas)
    plt.xlabel('$\Phi$')
    plt.title('AS')
    plt.gca().get_yaxis().set_visible(False)

    plt.subplot(1,3,3, sharex=ax, sharey=ax)
    V.plot_cells(Cnew, get_colors=get_colors, fill=True)
    plt.title('AS -> MD')
    plt.xlabel('$\Phi$')
    plt.gca().get_yaxis().set_visible(False)

    plt.savefig(opts.outfile, bbox_inches='tight')


if __name__ == '__main__':
    opts = getopts()
    main(opts)
