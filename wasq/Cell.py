
import pxul

import numpy as np

import itertools
import cPickle as pickle
import os


class Cells(object):
    def __init__(self, initial, labels):
        assert len(initial) == len(labels), '|cells| = {} but |labels| = {}'.format(len(initial),
                                                                                    len(labels))
        self._cells  = [initial]
        self._labels = [labels]
        self._count    = len(initial)
        self._dim    = initial.shape[-1]

        self._view_index = 1
        self._cells_view = initial
        self._labels_view = labels
        self._view_idx_to_chunk = [0 for _ in xrange(self._count)]

    def __len__(self):
        return self._count

    @property
    def _view_updated(self):
        "Predicate indicating that the views represent all the data"
        return self._view_index == len(self._cells)

    @property
    def shape(self):
        return self._count, self._dim

    @property
    def dim(self):
        "The dimension of the cells"
        return self._dim

    def _view_attr(self, attr):
        "Return the given attributed after updated the view if necessary"
        if not self._view_updated:
            self.update_view()
        return getattr(self, attr)

    def update_view(self):
        self._cells_view  = np.vstack([self._cells_view ] + self._cells[self._view_index:])
        self._labels_view = np.hstack([self._labels_view] + self._labels[self._view_index:])
        self._view_index = len(self._cells)

    @property
    def C(self):
        "View the cells into a single NxD array"
        return self._view_attr('_cells_view')

    @property
    def L(self):
        "View the labels into a single N array"
        return self._view_attr('_labels_view')

    @property
    def chunks(self):
        "Returns the number of chunks"
        return len(self._cells)

    def chunk_of(self, i):
        "Returns the chunk that the index into the view belongs to"
        return self._view_idx_to_chunk[i]

    def learn(self, cells, labels):
        assert len(cells) == len(labels), '|cells| = {} but |labels| = {}'.format(len(cells),
                                                                                  len(labels))

        assert cells.shape[1] == self._dim, 'Expect {dim}-D data but got {shape[1]}'.format(dim=self._dim, shape=cells.shape)

        for _ in cells: self._view_idx_to_chunk.append(len(self._cells))
        self._cells.append(cells)
        self._labels.append(labels)
        self._count += len(cells)

    def write_to_dir(self, prefix, cellname='C', labelname='L', suffix='pkl', padding=6, overwrite=False):
        maxfiles = int(padding*'9')
        assert len(self._cells) <= maxfiles, \
          'Padding {} limits number of files to {}, but {} are needed'.format(padding, maxfiles, len(self._cells))

        template = '{name}-{i:0%d}.{suffix}' % padding

        with pxul.os.StackDir(prefix):
            for i, c, l in itertools.izip(xrange(len(self._cells)), self._cells, self._labels):
                cn = template.format(name = cellname, i = i, suffix = suffix)
                ln = template.format(name = labelname,i = i, suffix = suffix)

                for p,a in [(cn, c), (ln, l)]:
                    if os.path.exists(p) and not overwrite: continue
                    with open(p, 'wb') as fd: pickle.dump(a, fd, protocol=pickle.HIGHEST_PROTOCOL)

    @classmethod
    def load_from_dir(cls, path, cellname='C', labelname='L', suffix='pkl'):
        import glob
        import cPickle as pickle

        def paths(name):
            paths = glob.glob(os.path.join(path, '{name}-*.{suffix}'.format(name=name, suffix=suffix)))
            paths.sort()
            return paths
        cellpaths  = paths(cellname)
        labelpaths = paths(labelname)
        assert len(cellpaths) > 0, 'No files found in {}'.format(path)
        assert len(cellpaths) == len(labelpaths), 'Number of labels files {} should equal number of cell files {}'.format(
            len(labelpaths), len(cellpaths))

        def load(path):
            if suffix == 'pkl':
                with open(path, 'rb') as fd: return pickle.load(fd)

        cells = None
        for cpath, lpath in itertools.izip(cellpaths, labelpaths):
            C = load(cpath)
            L = load(lpath)

            if cells is None:
                cells = cls(C, L)
            else:
                cells.learn(C, L)

        return cells
