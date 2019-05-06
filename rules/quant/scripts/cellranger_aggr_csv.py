#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import sys
import os
import argparse

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('input', help='sample molecule info files (hd5)', nargs='+')
parser.add_argument('-o', '--outdir', help='output directory', required=True)
parser.add_argument('--sample-info', help='samplesheet info, tab seprated file assumes `Sample_ID` in header', required=True)
parser.add_argument('--batch', default=None, help='Column name in sample-info describing batch (optional)')
parser.add_argument('--groupby', default='all', help='Column name in sample-info describing groups of sample to be aggregated')
parser.add_argument('-v ', '--verbose', help='verbose output.', action='store_true')


def read_sample_info(args):
    import pandas as pd
    df = pd.read_csv(args.sample_info_info, sep='\t')
    assert('Sample_ID' in df.columns)
    df.index = df['Sample_ID']
    if args.batch is not None:
        assert(args.batch in df.columns)
        if args.verbose:
            _batches = list(set(df[args.batch]))
            print('identified batch column: {} with {} unique batches'.format(args.batch, len(_batches)))
    if args.groupby is not None or args.groupby != 'all':
        assert(args.groupby in df.columns)
        _groups = list(set(df[args.groupby]))
        print('identified groupby column: {} with {} unique groups'.format(args.batch, len(_groups)))
    
    if args.verbose:
        print(df.head())
    return df

def write_csv(input_files, valid_samples, fh, batch=None):
    if batch:
        fh.write('library_id,molecule_h5,batch\n')
    else:
        fh.write('library_id,molecule_h5\n')
        
    for fn in input_files:
        sample = fn.split(os.path.sep)[-3]
        if sample in valid_samples:
            if batch is not None:
                fh.write('{},{},{}\n'.format(sample, mol_h5, df.loc[sample][batch]))
            else:
                fh.write('{},{}\n'.format(sample, mol_h5)) 
        else:
            print('{} not in {}'.format(sample, str(valid_samples)))
        if args.verbose:
            print('wrote file: {}'.format(fh.name))
    
if __name__ == '__main__':
    args = parser.parse_args()
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
        
    df = read_sample_info(args)

    if args.groupby is not None or args.groupby != 'all':
        groups = df.groupby(arg.groupby).groups
    else:
        groups = {'all': list(df.index)}

    for name, samples in groups.items():
        output_fn = os.path.join(outdir, '{}_cellranger_aggr.csv'.format(name))
        with open(output_fn) as fh:
            write_csv(args.input, samples, fh, args.batch)
