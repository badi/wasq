import numpy as np
import textwrap


def load_cells_gps(path):
    def load_header(fd):
        m_ncells = re.match(r'ncells: (\d+)', fd.readline())
        m_ncrd   = re.match(r'ncoords: (\d+)', fd.readline())
        m_dim    = re.match(r'ndims: (\d+)', fd.readline())

        if not m_ncells: raise Exception, 'Failed to parse ncells from %s' % path
        if not m_ncrd  : raise Exception, 'Failed to parse ncoords from %s' % path
        if not m_dim  : raise Exception, 'Failed to parse ndims from %s' % path

        return int(m_ncells.group(1)), int(m_ncrd.group(1)), int(m_dim.group(1))

    with open(path) as fd:
        ncells, ncoords, dim = load_header(fd)
        cells = np.loadtxt(fd).reshape((ncells, ncoords, dim))

    return cells


def save_cells_gps(fd, cells, ndx=None, dim=3):
    if ndx is None and len(cells) > 0:
        ndx = len(cells[0])

    centroids = -np.ones((len(cells), len(ndx), dim))
    for i, lbl in enumerate(cells.L):
        centroids[i] = lbl.x[ndx]

    ncells = len(cells)
    ncoords = len(ndx)

    header = textwrap.dedent("""\
       ncells: {}
       ncoords: {}
       ndims: {}

       """.format(ncells, ncoords, dim))
    fd.write(header)
    np.savetxt(fd, centroids.flatten(), fmt='%f')
