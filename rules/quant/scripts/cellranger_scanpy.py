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
          'GRCh38': 'GRCh38',
          'mus_musculus': 'mm10',
          'mouse': 'mm10',
          'mm10': 'mm10',
          'GRCm38': 'mm10'
}

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--infile', help='input filename')
parser.add_argument('-o', '--outfile', help="output filename",)
parser.add_argument('-f', '--format', choices=['anndata', 'loom', 'csvs'], default='anndata', help="output file format")
parser.add_argument('--sample-sheet', help="samplesheet filename")
parser.add_argument('--genome', help="Reference genome", default=None)
parser.add_argument('--filter-org', help="Filter data (genes) by organism", default=None)



if __name__ == '__main__':
    args = parser.parse_args()
    
    if not os.path.exists(args.infile):
        raise IOError

    if args.infile.endswith('.h5'):
        if not args.genome:
            raise ValueError('loading a .h5 file in scanpy requires the `genome` parameter')
        data = sc.read_10x_h5(args.infile, genome=args.genome)
    else:
        data = sc.read_10x_mtx(args.infile, gex_only=True, var_names='gene_ids')


    if args.filter_org is not None:
        org = args.filter_org
        if org in GENOME:
            org = GENOME[org]
        PREFIX = {'GRCh38': 'hg', 'GRCm38':'', 'mm10':''}
    
    if args.format == 'anndata':
        data.write(args.outfile)
    elif args.format == 'loom':
        data.write_loom(args.outfile)
    elif args.format == 'csvs':
        data.write_csvs(args.outpfile)
    else:
        raise ValueError
