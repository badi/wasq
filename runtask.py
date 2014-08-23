
from wasq.AdaptiveSampling import *

import cPickle as pickle



def getopts():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    p = ArgumentParser(formatter_class = ArgumentDefaultsHelpFormatter)
    p.add_argument('-i', '--ipath', metavar='FILE', default='walker.pkl', help='Input state path')
    p.add_argument('-o', '--opath', metavar='FILE', default='result.pkl', help='Output state path')
    return p.parse_args()

def main(opts):
    walker = pickle.load(open(opts.ipath, 'rb'))
    result = walker.run()
    pickle.dump(result, open(opts.opath, 'wb'), pickle.HIGHEST_PROTOCOL)
    print result[0]


if __name__ == '__main__':
    opts = getopts()
    main(opts)
