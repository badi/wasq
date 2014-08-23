
from wasq.AdaptiveSampling import *

import mdtraj

import cPickle as pickle
import copy




if __name__ == '__main__':
    import sys
    spath,toppath,opath = sys.argv[1:4]

    S = pickle.load(open(spath,'rb'))
    top = mdtraj.load(toppath)

    traj = copy.copy(top)
    for st in S:
        t = copy.copy(top)
        t.xyz = st.x
        traj = traj.join(t)

    print traj
    traj.save(opath)
