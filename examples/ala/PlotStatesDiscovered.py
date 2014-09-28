#!/usr/bin/env python

from wasq.AdaptiveSampling import *
import wasq.VoronoiPlot as V

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import argparse
import os
import glob


def getopts():
    a = argparse.ArgumentParser()
    a.add_argument('-a', '--as-dirs', nargs='+', required=True, default=['AS'])
    a.add_argument('-l', '--labels' , nargs='+')
    a.add_argument('-t', '--walker-time', type=float, required=True, help='Time (in ps) each walker simulates')
    a.add_argument('-o', '--outfile', default='cells.png')

    opts = a.parse_args()
    return opts

def plot_AS(as_dir, walker_time, label=None, color='blue'):
    cells = Cells.load_from_dir(os.path.join(as_dir, 'cells'))
    nwalkers = sorted(glob.glob(os.path.join(as_dir, 'iteration', '*', 'nwalkers.txt')))

    states  = np.zeros(cells.chunks)
    walkers = np.zeros(cells.chunks)
    times   = np.zeros(cells.chunks)
    for i, path in enumerate(nwalkers):
        j = i + 1
        w = int(open(path).read().strip())
        t = w * walker_time / 10.**3

        try:
            states [j] = sum(map(len, cells._cells[:j]))
            walkers[j] = w
            times  [j] = times[i] + t

        # for analysing a still-running AS
        except IndexError: break

    plt.plot(times, states, color=color, label=label)
    


def main(opts):
    color_map = mpl.colors.LinearSegmentedColormap.from_list('my_colors', 'blue cyan yellow orange red'.split(), N=len(opts.labels))
    for as_dir, label, color in zip(opts.as_dirs, opts.labels, range(len(opts.labels))):
        print as_dir, label
        plot_AS(as_dir, opts.walker_time, label=label, color=color_map(color))

    plt.xlabel('Time Simulated (ns)')
    plt.ylabel('States Discovered')
    plt.legend(title='Temp (K)', loc='lower right')
    plt.grid()
    plt.savefig(opts.outfile)


if __name__ == '__main__':
    opts = getopts()
    main(opts)
