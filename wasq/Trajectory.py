
import PoissonCover as PC
from metrics.dihedral_rmsd import dihedral_rmsd

import pxul

import mdtraj
import numpy as np


def calc_dihedral(traj, mdtraj_calculator):
    ind, radians = mdtraj_calculator(traj)
    degrees = radians * 180 / np.pi
    return degrees

def calc_phipsi(traj):
    phi = calc_dihedral(traj, mdtraj.compute_phi)
    psi = calc_dihedral(traj, mdtraj.compute_psi)
    return np.hstack([phi,psi])

def cover(traj, radius):
    phipsi = calc_phipsi(traj)
    I, C = PC.simple_poisson_cover(phipsi, radius, metric=dihedral_rmsd)
    return I, C



if __name__ == '__main__':
    import sys

    out = sys.argv[1]
    rad = float(sys.argv[2])
    pdb = sys.argv[3]
    tjs = sys.argv[4:]

    print 'Loading', tjs
    t = mdtraj.load(tjs, top=pdb)
    print t

    print 'Covering withh radius = ', rad
    I, C = cover(t, rad)

    print 'Saving results to', out
    with pxul.os.StackDir(out):
        np.savetxt('C.txt', C)
        np.savetxt('I.txt', I, fmt='%d')

        for k, i in enumerate(I):
            t[i].save('cell_{}.pdb'.format(k))
