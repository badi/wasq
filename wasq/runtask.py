
from wasq.AdaptiveSampling import *

import pxul

import cPickle as pickle
import glob
import sys
import tarfile



def getopts():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    p = ArgumentParser(formatter_class = ArgumentDefaultsHelpFormatter)
    p.add_argument('-i', '--ipath', metavar='FILE', default='walker.pkl', help='Input state path')
    p.add_argument('-o', '--opath', metavar='FILE', default='result.pkl', help='Output state path')
    return p.parse_args()

def main(opts):
    walker = pickle.load(open(opts.ipath, 'rb'))

    try:
        result = walker.run()
        retcode = 0

    except:
        files = glob.glob('*')
        files = filter(lambda f: os.access(f, os.R_OK), files)
        tarball = 'debug.tar.bz2'
        with tarfile.open(tarball, 'w:bz2') as tf:
            for path in files:
                tf.add(path, recursive=True)
        with open(tarball, 'rb') as fd:
            result = dict(debug = fd.read())
            retcode = 1

    pickle.dump(result, open(opts.opath, 'wb'), pickle.HIGHEST_PROTOCOL)
    sys.exit(retcode)


if __name__ == '__main__':
    opts = getopts()
    main(opts)
