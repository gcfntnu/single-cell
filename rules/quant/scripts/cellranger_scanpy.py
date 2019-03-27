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
parser.add_argument('--sample-sheet', help='samplesheet filename, tab seprated file assumes `Sample_ID` in header', default=None)
parser.add_argument('--feature-info', help='extra feature info filename, tab seprated file assumes `gene_id` in header', default=None)
parser.add_argument('--genome', help='reference genome', default=None)
parser.add_argument('--filter-org', help='filter data (genes) by organism', default=None)
parser.add_argument('--var-names', help='filter data (genes) by organism', default='gene_ids', choices=['gene_ids', 'gene_symbols'])
parser.add_argument('--gex-only', help='only keep `Gene Expression` data and ignore other feature types.', default=True)
parser.add_argument('--normalize', help='normalize depth across the input libraries', default='none', choices=['none', 'mapped'])
parser.add_argument('--verbose', help='verbose output.', action='store_true')


def downsample_gemgroup(data_list):
    """downsample data total read count to gem group with lowest total count
    """
    min_count = 1E99
    for i, data in enumerate(data_list):
        isum = data.X.sum()
        if isum < min_count:
            min_count = isum
            idx = i
    for j, data in enumerate(data_list):
        if j != idx:
            data = sc.pp.downsample(data, total_counts = min_count)
        data_list[j] = data
    return data_list

if __name__ == '__main__':
    aggr_input = False
    samples = []
    sample_info = None
    args = parser.parse_args()
    for fn in args.input:
        if not os.path.exists(fn):
            raise IOError('file does not exist! {}'.format(fn))
        
    if args.samples is not None:
        samples = args.samples.split(',')
        if len(samples) != len(args.input):
            msg = 'number of sampleslisted in --samples : {} does not match number of of inputs: {} '
            raise ValueError(msg.format(len(samples), len(args.input)))
        if len(samples) == 1 and samples[0] == 'aggr':
            aggr_input = True
            
    if args.sample_sheet is not None:
        sample_info = pd.read_csv(args.sample_sheet, sep='\t')
        if not 'Sample_ID' in sample_info.columns:
            raise ValueError('sample_sheet needs a column called `Sample_ID`')
        sample_info.index = sample_info['Sample_ID']
        
        if args.samples is not None:
            if all(i in sample_info.index for i in samples):
                sample_info = sample_info.loc[samples,:]
            else:
                if not aggr_input:
                    msg = 'samples not present in samplesheet: {}'
                    raise ValueError(msg.format(set(samples).difference(sample_info.index)))
    data_list = []
    for i, fn in enumerate(args.input):
        if fn.endswith('.h5'):
            if not args.genome:
                raise ValueError('loading a .h5 file in scanpy requires the `genome` parameter')
            data = sc.read_10x_h5(fn, genome=args.genome)
        else:
            #dirname = os.path.dirname(fn)
            data = sc.read_10x_mtx(fn, gex_only=args.gex_only, var_names=args.var_names, make_unique=True)

        # remove the gem group numbering if cells originate from one gem pool 
        cells = data.obs_names.copy()
        barcodes = [i.split('-')[0] for i in cells]
        if len(set(barcodes)) == len(barcodes):
            data.obs_names = barcodes
            
        if args.filter_org is not None:
            org = args.filter_org
            if org in GENOME:
                org = GENOME[org]
            PREFIX = {'GRCh38': 'hg', 'GRCm38':'mm10'}
            
            raise NotImplementedError
        
        if samples:
            gem_groups = [int(c.split('-')[-1])-1 for c in cells]
            gem_names = [samples[g] for g in gem_groups]
            data.obs['sample'] = gem_names
            
        if sample_info:
            df = sample_info.loc[gem_names,:]
            df.index = data.obs_names
            obs = pd.concat(data.obs, df, axis=1)
        data_list.append(data)

    if len(data_list) > 1:
        if args.normalize == 'mapped':
            data_list = downsample_gemgroup(data_list)

    data = data_list.pop(0)
    if len(data_list) > 0:
        if args.samples is not None:
            batch_key = 'Sample_ID'
            batch_categories = samples
        else:
            batch_key = 'batch'
            batch_categories = None    
        data = data.concatenate(*data_list, batch_key=batch_key, batch_categories=batch_categories)
        
    if 'gene_symbols' in data.var.columns:
        mito_genes = adata.var.gene_symbols.str.startswith('MT-')
        adata.obs['fraction_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    
    if args.format == 'anndata':
        data.write(args.outfile)
    elif args.format == 'loom':
        data.write_loom(args.outfile)
    elif args.format == 'csvs':
        data.write_csvs(args.outpfile)
    else:
        raise ValueError
