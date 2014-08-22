
import PoissonCover as PC

import pxul
import mdprep
from mdq.md.gmx import guamps_get, guamps_set, editconf, mdrun, tpr_set_scalar, tpr_get_scalar, disable_gromacs_backups

import mdtraj

import numpy as np
from scipy.spatial import distance

import array
import functools
import glob
import os
import shutil
import tempfile
import cPickle as pickle

__DEFAULT_METRIC = PC.DEFAULT_METRIC

def calc_dihedral(traj, mdtraj_calculator):
    ind, radians = mdtraj_calculator(traj)
    degrees = radians * 180 / np.pi
    return degrees

def calc_phipsi(traj):
    phi = calc_dihedral(traj, mdtraj.compute_phi)
    psi = calc_dihedral(traj, mdtraj.compute_psi)
    return np.hstack([phi,psi])

def dihedral_rmsd(X, Y):
    d = np.mod(X, 360) - np.mod(Y, 360)
    return np.sqrt((np.abs(d)**2).mean())

def load_guamps(path, dtype=str):
    data = open(path).read().strip().split('\n')

    if len(data) == 1: # scalar
        result = dtype(data[0])
        return result

    elif len(data) > 1: # vector
        assert len(data) > 4 # make sure the header is there
        cells, crds, dims = data[:3]
        cells, crds, dims = [
            int(line.strip().split()[-1])
            for line in [cells, crds, dims]
            ]
        values = data[4:]
        assert cells == 1
        assert len(values) == crds * dims
        values = np.array(values).astype(dtype)
        values = values.reshape((crds, dims))
        return values

    else:
        raise ValueError

def write_guamps(path, X):
    with open(path, 'w') as fd:
        try:
            cells = 1
            crds, dims = X.shape
            fd.write('ncells: %d\n' % cells)
            fd.write('ncoords: %d\n' % crds)
            fd.write('ndims: %d\n' % dims)
            fd.write('\n')
            np.savetxt(fd, X.flatten(), fmt='%f')
        except AttributeError:
            fd.write('%f\n' % X)


class SimulationState(object):
    def __init__(self, x=None, v=None, t=None):
        assert x is not None
        assert v is not None
        assert t is not None
        self.x = x
        self.v = v
        self.t = t

    @classmethod
    def from_tpr(cls, path):
        path = os.path.abspath(path)
        with pxul.os.TmpDir():
            x,v,t = 'x v t'.split()
            guamps_get(f=path, s='positions',  o=x)
            guamps_get(f=path, s='velocities', o=v)
            guamps_get(f=path, s='time',       o=t)
            return cls(x=load_guamps(x, dtype=float),
                       v=load_guamps(v, dtype=float),
                       t=load_guamps(t, dtype=float))

    @classmethod
    def from_trr(cls, path, frame):
        path = os.path.abspath(path)
        with pxul.os.TmpDir():
            x,v,t = 'x v t'.split()
            guamps_get(f=path, s='positions',  o=x, i=frame)
            guamps_get(f=path, s='velocities', o=v, i=frame)
            guamps_get(f=path, s='time',       o=t, i=frame)
            return cls(x=load_guamps(x, dtype=float),
                       v=load_guamps(v, dtype=float),
                       t=load_guamps(t, dtype=float))

    def writeX(self, path): write_guamps(path, self.x)
    def writeV(self, path): write_guamps(path, self.v)
    def writeT(self, path): write_guamps(path, self.t)


class GromacsWalker(object):
    def __init__(self, state, tpr=None, metric=dihedral_rmsd, threads=0):
        assert tpr is not None
        self.state = state
        self.tpr = tpr
        self._metric = metric
        self._threads = threads
        self._tmpdir = tempfile.mkdtemp()

    def __del__(self):
        import shutil
        shutil.rmtree(self._tmpdir)

    @classmethod
    def from_spec(cls, state, extra_files, extra_params):
        assert 'tpr' in extra_files
        return cls(state, **extra_params)

    @classmethod
    def from_tpr(cls, path, tpr=None, **kws):
        state = SimulationState.from_tpr(path)
        ref = tpr or path
        return cls(state, ref, **kws)

    @classmethod
    def from_trr(cls, path, frame=0):
        state = SimulationState.from_trr(path)
        return cls(state, self.tpr, metric=self._metric)

    def sample(self):
        x,v,t,tpr = 'x.gps v.gps t.gps topol.tpr'.split()
        # setup
        self.state.writeX(x)
        self.state.writeV(v)
        self.state.writeT(t)
        shutil.copy(self.tpr, tpr)

        # resume
        guamps_set(f=tpr, s='positions',  i=x)
        guamps_set(f=tpr, s='velocities', i=v)
        guamps_set(f=tpr, s='time',       i=t)

        # run
        trr = 'traj.trr'
        pdb = 'topol.pdb'
        mdrun(s=tpr, o=trr, c=pdb, nt=self._threads)

        # result
        traj = mdtraj.load(trr, top=pdb)
        state = np.array([SimulationState.from_trr(trr, i) for i in xrange(len(traj))])

        return traj, state

    def cover(self, traj, state, R, C, L, eps=0.000001):
        phipsi = calc_phipsi(traj)
        C, L = PC.labeled_online_poisson_cover(phipsi, R, L=state, C=C, CL=L, metric=self._metric)
        return C, L

    def run(self, R, C, L):
        with pxul.os.StackDir(self._tmpdir):
            traj, top = self.sample()
            return self.cover(traj, top, R, C, L)


class AdaptiveSampler(object):
    def __init__(self, R, C, S, iterations=float('inf'),
                 walker_class=GromacsWalker, extra_files=None, extra_params=None,
                 metric=dihedral_rmsd, checkpoint_dir='AS'):
        assert len(C) == len(S)
        self.R = R                           # :: float: the radius
        self.C = C                           # :: NxD array: N points in D dimensional space
        self.S = S                           # :: [SimulationState]: a label for each point of C
        self.current_iteration = 0           # :: int
        self.max_iterations = iterations     # :: int
        self.extra_files = extra_files or {} # :: dict(str -> filepath)
        self.extra_params= extra_files or {} # :: dict(str -> a)
        self.walker_class = walker_class     # :: has classmethod
                                             #        from_spec :: SimulationState -> dict(str->filepath) -> dict(str->a) -> obj
        self.metric = metric                 # :: a -> a -> a
        self.checkpoint_dir = checkpoint_dir # :: filepath

    @classmethod
    def from_tprs(cls, reference, initials, radius, iterations=float('inf'), **extra_params):
        C = []
        S = []
        inits = map(os.path.abspath, initials)
        with pxul.os.TmpDir(), disable_gromacs_backups():
            pdb = 'topol.pdb'
            for tpr in inits:
                editconf(f=tpr, o=pdb)
                traj = mdtraj.load(pdb)
                phipsi = calc_phipsi(traj)
                s = SimulationState.from_tpr(tpr)
                C.append(phipsi)
                S.append(s)
        C = np.vstack(C)
        S = np.hstack(S)
        extra_files = dict(tpr = reference)
        return cls(radius, C, S, iterations=iterations,
                   extra_files=extra_files, extra_params=extra_params)

    def select_by_kissing_number(self):
        """
        Select a subset of centroids to start new walkers from
        Accepts:
          C :: NxM array: the centroids
        Returns:
          I :: [Int]: indices into C to start new simulations from
        """
        ns = PC.neighborhood_size(self.C, self.R)
        _, dim = self.C.shape
        kissing_number = PC.KISSING_NUMBERS[dim]
        fringes = np.where(ns < kissing_number)[0]
        return fringes

    def select(self):
        return self.select_by_kissing_number()

    def _select(self):
        count = self.select()
        if len(count) < 1:
            raise StopIteration
        else:
            return count

    def iterate(self):
        "Run one iteration"
        idx = self._select()
        walkers = [self.walker_class.from_spec(st, self.extra_files, self.extra_params)
                   for st in self.S]

        # fan out
        print self.current_iteration, len(walkers)
        results = []
        for w in walkers:
            r =  w.run(self.R, self.C, self.S)
            results.append(r)

        # merge
        for Cw, Sw in results:
            self.C, self.S = PC.labeled_online_poisson_cover(Cw, self.R, L=Sw, C=self.C, CL=self.S, metric=self.metric)

    def write_log(self):
        iteration_dir = os.path.join(self.checkpoint_dir, 'iteration', '%05d' % self.current_iteration)
        with pxul.os.StackDir(iteration_dir):
            np.savetxt('C.txt', self.C)
            with open('S.pkl', 'wb') as fd:
                pickle.dump(self.S, fd, pickle.HIGHEST_PROTOCOL)

    def run(self, eps=0.000001):
        while self.current_iteration < self.max_iterations:
            self.write_log()
            try: self.iterate()
            except StopIteration: break
            self.current_iteration += 1


def getopts():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    p = ArgumentParser()
    p.add_argument('-r', '--radius', default=20)
    p.add_argument('-i', '--iterations', type=float, default=float('inf'))
    p.add_argument('tprs', metavar='TPR', nargs='+',
                   help='Coordinates for the initial states. The first one will be used to propagate simulation parameters.')
    opts =  p.parse_args()
    opts.tprs = map(os.path.abspath, opts.tprs)
    return opts

def main(opts):
    ref = opts.tprs[0]
    sampler = AdaptiveSampler.from_tprs(ref, opts.tprs, opts.radius, iterations=opts.iterations)
    sampler.run()

if __name__ == '__main__':
    opts = getopts()
    main(opts)
