
import PoissonCover as PC
from metrics.dihedral_rmsd import dihedral_rmsd
from Trajectory import calc_dihedral, calc_phipsi

import pxul
import mdprep
from mdq.md.gmx import guamps_get, guamps_set, editconf, mdrun, tpr_set_scalar, tpr_get_scalar, disable_gromacs_backups
import pwq

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
import copy

__DEFAULT_METRIC = PC.DEFAULT_METRIC


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


class AbstractWalker(object):
    def __init__(self, state, reference=None, reference_id=None, metric=dihedral_rmsd, **kws):
        assert reference is not None
        assert reference_id is not None

        self.state = state
        self.reference = reference
        self.reference_id = reference_id
        self._metric = metric

        self._kws = kws

    @property
    def local_reference(self):
        base, ext = os.path.splitext(os.path.basename(self.reference))
        name = 'ref{i}_{base}{ext}'.format(i=self.reference_id, base=base, ext=ext)
        return name

    def sample(self):
        """
        return :: tuple =
          trajectory :: mdtraj.Trajectory
          state      :: [SimulationState]
        """
        raise NotImplementedError

    def cover(self, traj, state, R, C, L, eps=0.000001):
        phipsi = calc_phipsi(traj)
        C, L = PC.online_poisson_cover(phipsi, R, L=state, Cprev=C, Lprev=L, metric=self._metric)
        return C, L

    def run(self, R, C, L, workarea=None):
        """
        Returns tuple ::
          C :: NxD array: centroids
          L :: N   array: labels
        """
        if workarea is None:
            dir_ctx = pxul.os.TmpDir
        elif type(workarea) is str:
            dir_ctx = lambda: pxul.os.StackDir(workarea)
        else:
            dir_ctx = lambda: pxul.os.StackDir(workarea())

        with dir_ctx():
            traj, state = self.sample()
            return self.cover(traj, state, R, C, L)


class GromacsWalker(AbstractWalker):
    def __init__(self, state, threads=0, **kws):
        super(GromacsWalker, self).__init__(state, **kws)
        self._threads = threads

    @classmethod
    def from_tpr(cls, path, **kws):
        state = SimulationState.from_tpr(path)
        return cls(state, reference=os.path.abspath(path), **kws)

    @classmethod
    def from_trr(cls, path, frame=0, **kws):
        state = SimulationState.from_trr(path, frame=frame)
        return cls(state, self.tpr, metric=self._metric)

    def sample(self):
        with disable_gromacs_backups():
            x,v,t = 'x.gps v.gps t.gps'.split()
            tpr = self.local_reference
            # setup
            self.state.writeX(x)
            self.state.writeV(v)
            self.state.writeT(t)
            if not os.path.exists(tpr):
                shutil.copy(self.reference, tpr)

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
            walkers = np.zeros(len(traj), dtype=np.object)
            for i in xrange(len(traj)):
                state = SimulationState.from_trr(trr, i)
                new_walker = copy.copy(self)
                new_walker.state = state
                walkers[i] = new_walker

            return traj, walkers

class AbstractAdaptiveSampler(object):
    def __init__(self, R, C, S, iterations=float('inf'),
                 walker_class=GromacsWalker, extra_files=None, extra_params=None,
                 metric=dihedral_rmsd, workarea='AS'):
        assert len(C) == len(S)
        self.R = R                           # :: float: the radius
        self.C = C                           # :: NxD array: N points in D dimensional space
        self.S = S                           # :: [SimulationState]: a label for each point of C
        self.current_iteration = 0           # :: int
        self.max_iterations = iterations     # :: int
        self.metric = metric                 # :: a -> a -> a
        self.workarea = workarea             # :: filepath

    @classmethod
    def from_tprs(cls, initials, radius, iterations=float('inf'),
                  walker_class=AbstractWalker, walker_kws=None, **init_kws):
        walker_kws = walker_kws or {}
        C = []
        W = np.zeros(len(initials), dtype=np.object)
        inits = map(os.path.abspath, initials)
        with pxul.os.TmpDir(), disable_gromacs_backups():
            pdb = 'topol.pdb'
            for tprid, tpr in enumerate(inits):
                editconf(f=tpr, o=pdb)
                traj = mdtraj.load(pdb)
                phipsi = calc_phipsi(traj)
                W[tprid] = walker_class.from_tpr(tpr, reference_id = tprid, **walker_kws)
                C.append(phipsi)
        C = np.vstack(C)

        return cls(radius, C, W, **init_kws)

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

    def current_walkers(self):
        return self.S

    def run_walker(self, walker):
        raise NotImplementedError

    def collect_results(self):
        raise NotImplementedError

    def iterate(self):
        "Run one iteration"
        selected = set(self._select())
        walkers = self.current_walkers()

        count_submitted = 0
        for i,w in enumerate(walkers):
            if i not in selected: continue
            self.run_walker(w)
            selected.remove(i)
            count_submitted += 1

        print self.current_iteration, count_submitted

        for Cw, Sw in self.collect_results():
            self.C, self.S = PC.online_poisson_cover(Cw, self.R, L=Sw, Cprev=self.C, Lprev=self.S, metric=self.metric)

    def write_log(self):
        iteration_dir = os.path.join(self.workarea, 'iteration', '%05d' % self.current_iteration)
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


class LocalAdaptiveSampler(AbstractAdaptiveSampler):
    "Run adaptive sampling on the local  machine"
    def __init__(self, *args, **kws):
        super(LocalAdaptiveSampler, self).__init__(*args, **kws)
        self._iteration_results = list()

    def run_walker(self, walker):
        r = walker.run(self.R, self.C, self.S)
        self._iteration_results.append(r)

    def collect_results(self):
        for r in self._iteration_results:
            yield r
        self._iteration_results = list()



class WorkQueueTaskWrapper(object):
    def __init__(self, obj, *args, **kws):
        self.obj = obj
        self.args = args
        self.kws = kws

    def run(self):
        return self.obj.run(*self.args, **self.kws)


class PythonTaskWorkQueueAdaptiveSampler(AbstractAdaptiveSampler):
    "Run using WorkQueue, but the Tasks are just pickled Walkers"

    def __init__(self, *args, **kws):
        super(PythonTaskWorkQueueAdaptiveSampler, self).__init__(*args, **kws)
        self._wq = None
        self.task_files_dir = os.path.join(self.workarea, 'task_files')
        pxul.os.ensure_dir(self.task_files_dir)

        self._walker_tmpl = 'walker.pkl'
        self._result_tmpl = 'result.pkl'

    def set_workqueue(self, wq):
        self._wq = wq

    def run_walker(self, walker):
        wrapped_walker = WorkQueueTaskWrapper(walker, self.R, self.C, self.S, workarea=os.getcwd)

        wasq_root = os.environ['WASQ_ROOT']
        runtask = os.path.join(wasq_root, 'wasq', 'runtask.py')

        t = pwq.Task('python runtask.py')
        t.specify_input_file(runtask, 'runtask.py', cache=True)

        walker_pkl = self.walker_path(t)
        result_pkl = self.result_path(t)

        pickle.dump(wrapped_walker, open(walker_pkl, 'wb'), pickle.HIGHEST_PROTOCOL)

        t.specify_input_file(walker_pkl, 'walker.pkl', cache=False)
        t.specify_output_file(result_pkl, 'result.pkl', cache=False)

        t.specify_input_file(walker.reference, walker.local_reference, cache=True)

        self._wq.submit(t)

    def walker_path(self, task): return os.path.join(self.task_files_dir, '%s.%s' % (self._walker_tmpl, task.uuid))
    def result_path(self, task): return os.path.join(self.task_files_dir, '%s.%s' % (self._result_tmpl, task.uuid))

    def collect_results(self):

        while not self._wq.empty():
            t = self._wq.wait(5)

            # success
            if t and t.result == 0:
                walker_pkl = self.walker_path(t)
                os.unlink(walker_pkl)

            # failure
            elif t and t.result != 0:
                msg = 'task %s failed with code %s\n' % (t.command, t.result)
                msg += t.output
                raise Exception, msg

        for result_pkl in glob.iglob(os.path.join(self.task_files_dir, '{}.*'.format(self._result_tmpl))):
            result = pickle.load(open(result_pkl, 'rb'))
            yield result
            os.unlink(result_pkl)


def test(opts):
    sampler = PythonTaskWorkQueueAdaptiveSampler.from_tprs(opts.ref, opts.tprs, opts.radius, iterations=opts.iterations)
    walkers = sampler.current_walkers()
    w = walkers.next()
    sampler.run_walker(w)
    

def getopts():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

    walker_choices = dict(gromacs = GromacsWalker)

    p = ArgumentParser()
    p.add_argument('-d', '--debug', default=False, help='Turn on debugging')
    p.add_argument('-p', '--port', default=9123, help='Start WorkQueue on this port')
    p.add_argument('-r', '--radius', default=20.0, type=float, help='Radius to use when covering data points')
    p.add_argument('-i', '--iterations', type=float, default=float('inf'), help='Number of AS iterations to run')
    p.add_argument('-t', '--type', default='gromacs', choices=walker_choices.keys())
    p.add_argument('tprs', metavar='TPR', nargs='+', help='Coordinates for the initial states.')

    opts = p.parse_args()
    opts.walker_class = walker_choices[opts.type]
    opts.tprs = map(os.path.abspath, opts.tprs)
    return opts

def main(opts):
    sampler = PythonTaskWorkQueueAdaptiveSampler.from_tprs(opts.tprs, opts.radius, iterations=opts.iterations,
                                                           walker_class=GromacsWalker, walker_kws=dict(threads=1))
    mkq = pwq.MkWorkQueue().replicate().port(opts.port)
    if opts.debug:
        mkq.debug()
    q = mkq()
    print 'WorkQueue running on', q.port
    sampler.set_workqueue(q)
    sampler.run()

if __name__ == '__main__':
    opts = getopts()
    main(opts)
