
import numpy as np
from scipy.spatial import distance

import functools



DEFAULT_METRIC = 'euclidean'
KISSING_NUMBERS = dict(zip(range(1,12), # dims
                           [2,6,12,24,40,71,126,240,272,336,438,756]))



class NeighborSearch(object):
    def __init__(self, data, metric='euclidean'):
        self._data = data
        self._metric = metric

    def query_ball_point(self, x, r):
        "Return the indices I of the dataset where d(x, data_i) <= r s.t. i \in I"
        d = distance.cdist([x], self._data, metric=self._metric)
        d = d[0] # since cdist([x],...)
        I = np.where(d <= r)[0]
        return I


def make_cdist(metric=DEFAULT_METRIC):
    @functools.wraps(distance.cdist)
    def wrapper(*args, **kws):
        kws['metric'] = metric
        return distance.cdist(*args, **kws)
    return wrapper


def _poisson_cover(S, R, metric=DEFAULT_METRIC):
    """
    Parameters:
      S :: NxD matrix -- data points
      R :: float      -- radius

    Returns:
      I :: [int]      -- indices into S, the selected 'centroids'
    """

    cdist = make_cdist(metric=metric)

    active = [0]
    output = set()
    avail  = set(range(1,len(S)))
    nbrs   = NeighborSearch(S, metric=metric)
    _,dim  = S.shape
    knums  = KISSING_NUMBERS[dim]

    while avail or active:
        # current target
        i = (avail or active).pop()
        avail.discard(i)
        x = S[i]

        # skip if d(x, O) < R
        if output and cdist([x], S[list(output)]).min() < R: continue
        output.add(i)

        # get all neighbors n of x where R < dist(n,x) <= 2R
        ns = np.array(nbrs.query_ball_point(x, 2*R))
        ds = cdist([x], S[ns])
        js = np.where((R<ds)&(ds<=2*R))[0]

        # if no neighbors, ignore these coordinates in the future
        if (len(js) < 1):
            for n in ns: avail.discard(n)
            continue

        ns = ns[js]
        ds = ds[js]

        # consider up to K neighbors (R < d <= 2R) of x
        kf = 0
        for n in ns:
            d = cdist([S[n]], S[list(output)]).min()
            if d > R:
                kf += 1
                active.append(n)
            avail.discard(n)
            if kf >= knums: break

    I = np.array(list(output))
    return I

def simple_poisson_cover(S, R, metric=DEFAULT_METRIC):
    I = _poisson_cover(S, R, metric=metric)
    return I, S[I]

def labeled_poisson_cover(S, R, L, metric=DEFAULT_METRIC):
    I, S2 = simple_poisson_cover(S, R, metric=metric)
    L2    = L[I]
    return I, S2, L2


def poisson_cover(S, R, L=None, metric=DEFAULT_METRIC):
    """
    Parameters:
      S :: NxD matrix of float -- N data points in D dimensions
      R :: float               -- the radius
      L :: Option [a]          -- optional labels for each element in S. Required: |S| == |L|

    Returns:
      (I, Sc, Lc) :: ([int], CxD matrix, Option [a])
        where
          I  -- the indices into S
          Sc -- the values of S at I
          Lc -- the values of L at I
    """

    I, Sc = simple_poisson_cover(S, R, metric=metric)
    Lc    = None

    if L is not None:
        assert len(L) == len(S), '|S| = {} but |L| = {}'.format(len(S), len(L))
        Lc = L[I]

    return I, Sc, Lc


def assign(S, C, metric=DEFAULT_METRIC):
    cdist = make_cdist(metric=metric)
    A = -np.ones(len(S), dtype=np.int)
    for i, crd in enumerate(S):
        d = cdist([crd], C)
        A[i] = d.argmin()
    return A

def neighborhood_size(C, R, metric=DEFAULT_METRIC):
    t = NeighborSearch(C, metric=metric)
    ns = np.zeros(len(C), dtype=np.int32)
    for i, x  in enumerate(C):
        ixs = t.query_ball_point(x, 2*R)
        ns[i] = len(ixs) - 1
    return ns


def get_fringes(C, R, metric=DEFAULT_METRIC):
    ns = neighborhood_size(C, R, metric=metric)
    _, dim = C.shape
    kissing_number = KISSING_NUMBERS[dim]
    fringes = C[ns < kissing_number]
    return fringes


def online_poisson_cover(S, R, L=None, Cprev=None, Lprev=None, metric=DEFAULT_METRIC):
    """
    Parameters:
      S
      R
      L
      Cprev
      Lprev

    Returns:
      Cnew, Lnew
    """

    def subsample_labels(L, I):
        if L is not None: return L[I]
        else: pass

    cdist = make_cdist(metric=metric)

    if Cprev is None or len(Cprev) == 0:
        _, Cnew, Lnew = poisson_cover(S, R, L=L, metric=metric)
        return Cnew, Lnew
    else:
        D = cdist(S, Cprev)
        new_points = [] # find all new points (further than R away from Cprev)
        for i in xrange(len(S)):
            a = D[i].argmin()
            if D[i,a] > R:
                new_points.append(i)

        if len(new_points) <= 0: return Cprev, Lprev

        _, Cn, Ln = poisson_cover(S[new_points], R, L=L[new_points], metric=metric)
        CS = np.vstack([Cprev, Cn])
        LS = np.hstack([Lprev, Ln])
        _, Cnew, Lnew = poisson_cover(CS, R, L=LS, metric=metric)
        return Cnew, Lnew
