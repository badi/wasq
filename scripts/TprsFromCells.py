
from wasq.AdaptiveSampling import *

from mdq.md.gmx import guamps_set, disable_gromacs_backups
import pxul

import cPickle as pickle
import copy
import shutil

if __name__ == '__main__':
    import sys
    spath,tprpath,odir = sys.argv[1:4]

    tprpath = os.path.abspath(tprpath)
    odir = os.path.abspath(odir)

    S = pickle.load(open(spath,'rb'))

    with pxul.os.StackDir(odir), disable_gromacs_backups():
        for i, s in enumerate(S):
            tpr = 'centroid_%05d.tpr' % i
            print tpr

            shutil.copy(tprpath, tpr)

            x,v,t = 'x.gps v.gps t.gps'.split()
            s.writeX(x)
            s.writeV(v)
            s.writeT(t)

            guamps_set(f=tpr, s='positions',  i=x)
            guamps_set(f=tpr, s='velocities', i=v)
            guamps_set(f=tpr, s='time',       i=t)

            map(os.unlink, [x,v,t])
