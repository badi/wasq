
import numpy as np

def dihedral_rmsd(X, Y):
    d = np.mod(X, 360) - np.mod(Y, 360)
    return np.sqrt((np.abs(d)**2).mean())


def run_test(shape=(1,2)):
    X = np.random.random(shape)
    Y = np.random.random(shape)
    d = dihedral_rmsd(X, Y)

if __name__ == '__main__':

    import timeit
    print timeit.timeit('run_test()',
                        number = 100000,
                        setup='from __main__ import run_test')
