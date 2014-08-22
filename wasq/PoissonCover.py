
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


def poisson_cover(S, R, metric=DEFAULT_METRIC):

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

    O = np.array(list(output))
    return S[O]

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


def online_poisson_cover(S, R, C=None, metric=DEFAULT_METRIC, **kws):
    cdist = make_cdist(metric=metric)
    if C is None or len(C) == 0:
        return poisson_cover(S, R, metric=metric, **kws)
    else:
        D = cdist(S, C)
        to_run_on = [] # find all points further than R from C
        for i in xrange(len(S)):
            a = D[i].argmin()
            if D[i,a] > R:
                to_run_on.append(i)
        if len(to_run_on) > 0:
            C1 = poisson_cover(S[to_run_on], R, metric=metric, **kws)
            CS = np.vstack([C, C1])
            C2 = poisson_cover(CS, R, metric=metric)
            return C2
        else:
            return C



################################################################################

def labeled_poisson_cover(S, R, L=None, metric=DEFAULT_METRIC):

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

    O = np.array(list(output))
    return S[O], L[O]




def labeled_online_poisson_cover(S, R, L=None, C=None, CL=None, metric=DEFAULT_METRIC, **kws):
    cdist = make_cdist(metric=metric)
    if C is None or len(C) == 0:
        return poisson_cover(S, R, L=L, metric=metric, **kws)
    else:
        D = cdist(S, C)
        new_points = [] # find all new points (further than R from C)
        for i in xrange(len(S)):
            a = D[i].argmin()
            if D[i,a] > R:
                new_points.append(i)
        if len(new_points) > 0:
            Cn, Ln = labeled_poisson_cover(S[new_points], R, L=L[new_points], metric=metric, **kws)
            CS = np.vstack([C, Cn])
            LS = np.hstack([CL, Ln])
            C2, L2 = labeled_poisson_cover(CS, R, L=LS, metric=metric)
            return C2, L2
        else:
            return C, CL
