#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import sys
import os
import argparse

import scanpy  as sc

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

parser.add_argument('input', help='input file(s)', nargs='+')
parser.add_argument('-o', '--outfile', help='output filename', required=True)
parser.add_argument('-f', '--format', choices=['anndata', 'loom', 'csvs'], default='anndata', help='output file format')
parser.add_argument('--samples', help='comma separated list of sample names. Used as `batch_categories` in anndata. Also used to match `Sample_ID` in optional samplesheet', default=None)
parser.add_argument('--sample-sheet', help='samplesheet filename', default=None)
parser.add_argument('--genome', help='reference genome', default=None)
parser.add_argument('--filter-org', help='filter data (genes) by organism', default=None)
parser.add_argument('--var-names', help='filter data (genes) by organism', default='gene_ids', choices=['gene_ids', 'gene_symbols'])
parser.add_argument('--gex-only', help='only keep 'Gene Expression' data and ignore other feature types.', default=True)
parser.add_argument('--verbose', help='verbose output.', action='store_true')



if __name__ == '__main__':
    args = parser.parse_args()

    for fn in args.input:
        if not os.path.exists(fn):
            raise IOError('file does not exist! {}'.format(fn))
    if args.samples is not None:
        samples = args.samples.split(',')
        if len(samples) != len(args.input):
            msg = 'number of sampleslisted in --samples : {} does not match number of of inputs: {} '
            raise ValueError(msg.format(len(samples), len(args.input)))
    if args.sample_sheet is not None:
        sample_info = pd.read_csv(args.sample_sheet, sep='\t')
        if not 'Sample_ID' in sample_info.columns:
            raise ValueError('sample_sheet needs a column called `Sample_ID`')
        sample_info.index = sample_info['Sample_ID']
        if args.samples is not None:
            if all(i in sample_info.index for i in samples):
                sample_info = sample_info.loc[samples,:]
            else:
                msg = 'samples not present in samplesheet: {}'
                raise ValueError(msg.format(set(samples).difference(sample_info.index)))
    data_list = []
    for fn in args.input:
        if fn.endswith('.h5'):
            if not args.genome:
                raise ValueError('loading a .h5 file in scanpy requires the `genome` parameter')
            data = sc.read_10x_h5(fn, genome=args.genome)
        else:
            dirname = os.path.dirname(fn)
            data = sc.read_10x_mtx(dirname, gex_only=args.gex_only, var_names=args.var_names, make_unique=True)

        # remove the gem channel numbering if cells originate from one gem pool 
        cells = data.obs_names.copy()
        barcodes = [i.split('-')[0] for i in cells]
        if len(set(barcodes)) == len(barcodes):
            data.obs_names = barcodes
            
        if args.filter_org is not None:
            org = args.filter_org
            if org in GENOME:
                org = GENOME[org]
            PREFIX = {'GRCh38': 'hg', 'GRCm38':'', 'mm10':''}
            raise NotImplementedError
        
        data_list.append(data)

    data = data_list.pop(0)
    if len(data_list) > 0:
        if args.samples is not None:
            batch_key = 'Sample_ID'
            batch_categories = samples
        else:
            batch_key = 'batch'
            batch_categories = None    
        data = data.concatenate(*data_list, batch_key=batch_key, batch_categories=batch_categories)
    
    if args.format == 'anndata':
        data.write(args.outfile)
    elif args.format == 'loom':
        data.write_loom(args.outfile)
    elif args.format == 'csvs':
        data.write_csvs(args.outpfile)
    else:
        raise ValueError
