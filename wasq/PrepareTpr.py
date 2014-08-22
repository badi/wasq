
import pxul
import mdprep
from mdq.md.gmx import guamps_get, guamps_set, tpr_set_scalar, tpr_get_scalar, disable_gromacs_backups

import numpy as np

import argparse
import os
import shutil



def set_tpr_values(keys, values, tpr):
    assert len(keys) == len(values)
    for k,v in zip(keys,values):
        tpr_set_scalar(tpr, k, v)

def set_tpr_outputfreq(freq, tpr):
    for attr in 'nstxout nstxtcout nstfout nstvout nstlog'.split():
        tpr_set_scalar(tpr, attr, freq)


def prepare(input_pdb, target_tpr, ff='amber03', water='tip3p', time_ps=10, outputfreq_ps=1):
    print input_pdb
    pdb = os.path.abspath(input_pdb)
    dst = os.path.abspath(target_tpr)
    odir = os.path.dirname(os.path.abspath(target_tpr))

    if water == 'none':
        preparer = mdprep.gmx_prep.PrepareImplicitSolventSystem
    else:
        preparer = mdprep.gmx_prep.PrepareSolvatedSystem

    with pxul.os.TmpDir(), disable_gromacs_backups():
        prep = preparer()
        res = prep.prepare(pdb, ff=ff, seed=np.random.randint(9999999))
        tpr = res['tpr']

        dt = tpr_get_scalar(tpr, 'deltat', float)
        nsteps = int(time_ps / dt)
        freq   = int(outputfreq_ps / dt)
        tpr_set_scalar(tpr, 'nsteps', nsteps)
        set_tpr_outputfreq(freq, tpr)

        shutil.copy(tpr, dst)


def getopts():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-f', '--ff', default='amber03', help='Forcefield')
    p.add_argument('-w', '--water', default='tip3p', help='Water model')
    p.add_argument('-T', '--simulation-time', default=10, type=float, help='Simulation time of final tpr in ps')
    p.add_argument('-O', '--output-freq', default=1, type=float, help='Trajectory frequency in ps')
    p.add_argument('-p', '--pdb', required=True)
    p.add_argument('-o', '--tpr', required=True)
    p.add_argument('-d', '--debug', default='critical', choices='debug info2 info1 info warning error critical'.split())
    return p.parse_args()


def main(opts):
    loglevel = 'set_' + opts.debug
    getattr(pxul.logging, loglevel)()
    prepare(opts.pdb, opts.tpr, ff=opts.ff, water=opts.water, time_ps=opts.simulation_time, outputfreq_ps=opts.output_freq)



if __name__ == '__main__':
    opts = getopts()
    main(opts)
