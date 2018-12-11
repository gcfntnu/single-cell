#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import sys
import os
import argparse

import scanpy.api as sc

GENOME = {'homo_sapiens': 'GRCh38',
          'human': 'GRCh38',
          'hg38': 'GRCh38',
          'mus_musculus': 'mm10',
          'mouse': 'mm10',
          'mm10': 'mm10',
          '': ''
}

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--infile', help='input filename')
parser.add_argument('-o', '--outfile', help="output filefilename",)
parser.add_argument('--genome', help="Reference genome", default=None)



if __name__ == '__main__':
    args = parser.parse_args()
    print(args)
    
    if not os.path.exists(args.infile):
        raise IOError

    if args.infile.endswith('.h5'):
        if not args.genome:
            raise ValueError('loading a .h5 file in scanpy requires the `genome` parameter')
        data = sc.read_10x_h5(args.infile, genome=args.genome)
    else:
        data = sc.read_10x_mtx(args.infile, gex_only=False)
    
    data.var_names_make_unique()
    data.write(args.outfile)
