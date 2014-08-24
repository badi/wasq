
import dihedral_rmsd_py as py
import dihedral_rmsd as cy

import numpy as np

import timeit

def get_arrays(shape=(2,), seed=42):
    np.random.seed(seed)
    X = np.random.random(shape)
    Y = np.random.random(shape)
    return X,Y

def test(module):
    X,Y = get_arrays()
    return module.dihedral_rmsd(X, Y)

def bench(module, number=1000):
    return timeit.timeit('test(%s)' % module,
                         number = number,
                         setup='from __main__ import test, %s' % module)

if __name__ == '__main__':

    X, Y = get_arrays()
    print py.dihedral_rmsd(X, Y)
    print cy.dihedral_rmsd(X, Y)

    p = bench('py')
    c = bench('cy')

    print
    print 'P', p
    print 'C', c
    print 'speedup', p / c
