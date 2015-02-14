#!/usr/bin/env python

from wasq.AdaptiveSampling import *
from wasq.util import save_cells_gps

import argparse

def getopts():
    a = argparse.ArgumentParser()

    a.add_argument('-c', '--cells', default='AS/cells', help='Path to AS cells directory')
    a.add_argument('-o', '--out',  default='cells.dat', help='Path to AWE cells file')
    a.add_argument('-n', '--ndx', default='atoms.ndx', help='Atom indices file to select atom coordinates')

    return a.parse_args()


def main(opts):

    print 'Loading cells from', opts.cells
    cells = Cells.load_from_dir(opts.cells)
    ndx = np.loadtxt(opts.ndx, dtype=np.int)
    with open(opts.out, 'w') as fd:
        save_cells_gps(fd, cells, ndx=ndx)


if __name__ == '__main__':
    main(getopts())
