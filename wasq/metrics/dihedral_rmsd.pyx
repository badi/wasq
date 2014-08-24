import numpy as np
cimport numpy as np

cdef extern from "math.h":
     double fabs(double)
     double pow(double,double)
     double sqrt(double)


cdef double c_dihedral_rmsd(np.ndarray[double,ndim=1] X, np.ndarray[double,ndim=1] Y):
    cdef Py_ssize_t i, size = len(X)
    cdef double Z = 0, t = 0

    t = 0
    Z = 0

    for i in range(size):
        t = fabs(X[i] % 360) - (Y[i] % 360)
        t = pow(t,2)
        Z = Z + t

    return sqrt(Z / size)


def dihedral_rmsd(np.ndarray[double, ndim=1] X not None, np.ndarray[double, ndim=1] Y not None):
    return  c_dihedral_rmsd(X,Y)
