
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


def set_tpr_values(keys, values, tpr):
    assert len(keys) == len(values)
    for k,v in zip(keys,values):
        tpr_set_scalar(tpr, k, v)

def set_tpr_outputfreq(freq, tpr):
    for attr in 'nstxout nstxtcout nstfout nstvout nstlog'.split():
        tpr_set_scalar(tpr, attr, freq)

def prepare_tpr(tpr, outdir='.', nsteps=10000, outfreq=1000):
    pxul.os.ensure_dir(outdir)
    set_tpr_outputfreq(outfreq, tpr)
    tpr_set_scalar(tpr, 'nsteps', nsteps)
    x = os.path.join(outdir,'x.gps')
    v = os.path.join(outdir,'v.gps')
    t = os.path.join(outdir,'t.gps')
    x,v,t = map(os.path.normpath, [x,v,t])
    guamps_get(f=tpr, s='positions' , o=x)
    guamps_get(f=tpr, s='velocities', o=v)
    guamps_get(f=tpr, s='time'      , o=t)
    return x,v,t

def sample(x,v,t,tpr, time_ps=100, outputfreq_ps=10, workarea='.'):
    x,v,t = map(lambda p: os.path.relpath(p, workarea),
                [x,v,t])

    with pxul.os.StackDir(workarea), disable_gromacs_backups():
        # resume the simulation
        guamps_set(f=tpr, s='positions', i=x)
        guamps_set(f=tpr, s='velocities',i=v)
        guamps_set(f=tpr, s='time',      i=t)

        # make sure nsteps and output frequency are correct
        dt = tpr_get_scalar(tpr, 'deltat', float)
        nsteps = int(time_ps / dt)
        tpr_set_scalar(tpr, 'nsteps', nsteps)
        freq = int(outputfreq_ps / dt)
        set_tpr_outputfreq(freq, tpr)

        # run the simulation
        trr = 'traj.trr'
        pdb = 'topol.pdb'
        mdrun(s=tpr, o=trr, c=pdb, nt=0)


def prepare(pdbs, workarea='AS', ff='amber03', nsteps=10000, outfreq=1000, preparer=None):
    preparer = preparer or mdprep.gmx_prep.PrepareImplicitSolventSystem
    prefix = os.path.join(workarea, 'previous')

    for i, pdb in enumerate(pdbs):

        pdb = os.path.abspath(pdb)
        walker = 'walker.%d' % i
        outdir = os.path.join(prefix, walker)
        pxul.os.ensure_dir(outdir)

        with pxul.os.StackDir(outdir):
            prep = preparer()
            res  = prep.prepare(pdb, ff=ff, seed=np.random.randint(999999))
            tpr  = res['tpr']
            xvt  = prepare_tpr(tpr, nsteps=nsteps, outfreq=outfreq)
            tpr2 = 'topol.tpr'
            os.rename(tpr, tpr2)

        x,v,t,tpr = map(functools.partial(os.path.join, outdir),
                        xvt + (tpr2,))
        nxt = os.path.join(workarea, 'next', walker)
        pxul.os.ensure_dir(nxt)

        result = dict()
        for key, path in zip('x v t tpr'.split(), [x,v,t,tpr]):
            base = os.path.basename(path)
            dst = os.path.join(nxt, base)
            shutil.copy(path, dst)
            pxul.logging.logger.info(path, '->', dst)
            result[key] = base

        tpr_set_scalar(tpr, 'nsteps', nsteps)
        set_tpr_outputfreq(outfreq, tpr)

        yield nxt, result



def sample_walker(walker_dir, time_ps, outfreq_ps=10):
    with pxul.os.StackDir(walker_dir):
        sample('x.gps', 'v.gps', 't.gps', 'topol.tpr', time_ps=time_ps, outputfreq_ps=outfreq_ps)


def calc_dihedral(traj, mdtraj_calculator):
    ind, radians = mdtraj_calculator(traj)
    degrees = radians * 180 / np.pi
    return degrees

def calc_phipsi(traj):
    phi = calc_dihedral(traj, mdtraj.compute_phi)
    psi = calc_dihedral(traj, mdtraj.compute_psi)
    return np.hstack([phi,psi])

def cover_walker(walker_dir, radius, C=None, metric=__DEFAULT_METRIC):
    with pxul.os.StackDir(walker_dir):
        trj = mdtraj.load('traj.xtc', top='topol.pdb')
        S = calc_phipsi(trj)
        C = PC.online_poisson_cover(S, radius, C=C, metric=metric)
        A = PC.assign(S, C, metric=metric)
        return C, A, trj


def step_walker(i, radius, workarea='AS', C=None, metric=__DEFAULT_METRIC, time_ps=100, outfreq_ps=10):
    walker_dir = os.path.join(workarea, 'next', 'walker.%d' % i)
    sample_walker(walker_dir, time_ps=time_ps, outfreq_ps=outfreq_ps)
    C, A, trj = cover_walker(walker_dir, radius, C=C, metric=metric)
    return C, A, trj


def initialize(starting_pdbs, radius=10):
    angles = []
    for pdb in starting_pdbs:
        trj = mdtraj.load(pdb)
        ang = calc_phipsi(trj)
        angles.append(ang)
    angles = np.vstack(angles)
    return PC.poisson_cover(angles, radius)

def main_loop(walker_ids, radius, workarea='AS', C=None, metric=__DEFAULT_METRIC, time_ps=100, outfreq_ps=10, eps=0.000001):
    cdist = PC.make_cdist(metric=metric)

    C = C
    trajs = []
    nexts = [] # :: [(walker_id, traj idx)]

    for i in walker_ids:
        C, _, trj = step_walker(i, radius, workarea=workarea, C=C, metric=metric, time_ps=time_ps, outfreq_ps=outfreq_ps)
        trajs.append(trj)

    fringes = PC.get_fringes(C, radius, metric=metric)

    # reassign to find all centroids
    for fringe in fringes:
        for i in walker_ids:
            phipsi = calc_phipsi(trajs[i])
            d = cdist([fringe], phipsi)
            if np.abs(d.min()) <= eps:
                nexts.append((i, d.argmin()))
                break

    # prepare for next iteration
    new_walkers = []
    with pxul.os.StackDir(workarea):
        shutil.rmtree('previous')
        os.rename('next', 'previous')
        new_walker = 0
        for i, k in nexts:
            prv = os.path.join('previous', 'walker.%d' % i)
            nxt = os.path.join('next',     'walker.%d' % new_walker)
            pxul.os.ensure_dir(nxt)
            shutil.copy(os.path.join(prv, 'topol.tpr'),
                        os.path.join(nxt, 'topol.tpr'))

            # extract x, v, t of centroids from trr
            trr = os.path.join(prv, 'traj.trr')
            x,v,t = map(functools.partial(os.path.join, nxt),
                        'x.gps v.gps t.gps'.split())
            guamps_get(f=trr, s='positions', o=x, i=k)
            guamps_get(f=trr, s='velocities',o=v, i=k)
            guamps_get(f=trr, s='time',      o=t, i=k)

            new_walkers.append(new_walker)
            new_walker += 1

    return C, new_walkers


def dihedral_rmsd(X, Y):
    d = np.mod(X, 360) - np.mod(Y, 360)
    return np.sqrt((np.abs(d)**2).mean())

def main(opts):
    if opts.water is 'none':
        preparer = mdprep.gmx_prep.PrepareImplicitSolventSystem
    else:
        preparer = mdprep.gmx_prep.PrepareSolvatedSystem

    workarea = 'AS'
    if os.path.exists(workarea):
        shutil.rmtree(workarea)

    pxul.logging.set_warning()

    W = range(len(opts.pdbs))
    # prepare is lazy
    for _ in prepare(opts.pdbs, ff=opts.ff, nsteps=opts.steptime*1000, outfreq=opts.outputfreq*1000, preparer=preparer): continue
    C = initialize(opts.pdbs, opts.radius)
    s0, s1 = len(C), len(C)
    statedir = os.path.join(workarea, 'state')
    statetmpl = os.path.join(statedir, '%05d.dat')
    pxul.os.ensure_dir(statedir)
    np.savetxt(statetmpl % 0, C)

    for i in xrange(opts.iterations):
        print i, len(W)
        C1, W = main_loop(W, opts.radius, C=C, metric=dihedral_rmsd, time_ps=opts.steptime, outfreq_ps=opts.outputfreq)
        np.savetxt(statetmpl % (i+1), C1)
        if C.shape == C1.shape and np.abs(C-C1).max() == 0:
            print 'Converged'
            break
        C = C1


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


def initialize_walkers(reference_tpr, initial_tprs):
    return [
        GromacsWalker.from_tpr(tpr, reference=reference_tpr)
        for tpr in initial_tprs
        ]

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


def test():
    tpr = '/tmp/ala.tpr'
    pdb = '/tmp/ala.pdb'
    s = AdaptiveSampler.from_tprs(tpr, [tpr], 10)
    return s


def getopts():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    p = ArgumentParser()#formatter_class=ArgumentDefaultsHelpFormatter)
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
