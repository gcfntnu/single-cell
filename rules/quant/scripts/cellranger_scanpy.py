#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import sys
import os
import argparse

import scanpy  as sc
import pandas as pd
import numpy as np
from vpolo.alevin import parser as alevin_parser

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
parser.add_argument('-f', '--input-format', choices=['cellranger_aggr', 'cellranger', 'star', 'alevin', 'umitools'], default='cellranger_aggr', help='input file format')
parser.add_argument('-F', '--output-format', choices=['anndata', 'loom', 'csvs'], default='anndata', help='output file format')
parser.add_argument('--sample-info', help='samplesheet info, tab seprated file assumes `Sample_ID` in header', default=None)
parser.add_argument('--feature-info', help='extra feature info filename, tab seprated file assumes `gene_id` in header', default=None)
parser.add_argument('--genome', help='reference genome', default=None)
parser.add_argument('--filter-org', help='filter data (genes) by organism', default=None)
parser.add_argument('--gex-only', help='only keep `Gene Expression` data and ignore other feature types.', default=True)
parser.add_argument('--normalize', help='normalize depth across the input libraries', default='none', choices=['none', 'mapped'])
parser.add_argument('--batch', help='column name in `sample-info` with batch covariate', default=None)
parser.add_argument('--no-zero-cell-rm', help='do not remove cells with zero counts', action='store_true')
parser.add_argument('-v ', '--verbose', help='verbose output.', action='store_true')


def downsample_gemgroup(data_list):
    """downsample data total read count to gem group with lowest total count
    """
    min_count = 1E99
    sampled_list = []
    for i, data in enumerate(data_list):
        isum = data.X.sum()
        if isum < min_count:
            min_count = isum
            idx = i
    for j, data in enumerate(data_list):
        if j != idx:
            sc.pp.downsample_counts(data, target_counts = min_count)
        sampled_list.append(data)
    return sampled_list

def read_cellranger(fn, args, rm_zero_cells=True, **kw):
    """read cellranger results

    Assumes the Sample_ID may be extracted from cellranger output dirname, 
    e.g ` ... /Sample_ID/outs/filtered_feature_bc_matrix.h5 `
    """
    if fn.endswith('.h5'):
        if not args.genome:
            raise ValueError('loading a .h5 file in scanpy requires the `genome` parameter')
        dirname = os.path.dirname(fn)
        data = sc.read_10x_h5(fn, genome=GENOME[args.genome])
        data.var['gene_symbols'] = list(data.var_names)
        data.var_names = list(data.var['gene_ids'])
    else:
        mtx_dir = os.path.dirname(fn)
        dirname = os.path.dirname(mtx_dir)
        data = sc.read_10x_mtx(mtx_dir, gex_only=args.gex_only, var_names='gene_ids', make_unique=True)
        data.var['gene_ids'] = list(data.var_names)
    sample_id = os.path.basename(os.path.dirname(dirname))
    data.obs['library_id'] = [sample_id] * data.obs.shape[0]
    barcodes = [b.split('-')[0] for b in data.obs.index]
    if len(barcodes) == len(set(barcodes)):
        data.obs_names = barcodes
    data.obs.index.name = 'barcodes'
    data.var.index.name = 'gene_id'
    return data
        
def read_cellranger_aggr(fn, args, **kw):
    data = read_cellranger(fn, args)
    if 'library_id' in data.obs:
        data.obs.rename(index=str, columns={'library_id': 'group'}, inplace=True)
    dirname = os.path.dirname(fn)
    if not fn.endswith('.h5'):
        dirname = os.path.dirname(dirname)

    # fix cellranger aggr enumeration to start at 0 (matches scanpy enum)
    barcodes =  [i[0] for i in data.obs.index.str.split('-')]
    if any(data.obs.index.str.contains('-')):
        barcodes_enum = [str(int(i[1])-1) for i in data.obs.index.str.split('-')]
    else:
        barcodes_enum = ['0'] * len(barcodes)
    data.obs_names = ['-'.join(e) for e in zip(barcodes, barcodes_enum)]

    aggr_csv = os.path.join(dirname, 'aggregation.csv')
    if os.path.exists(aggr_csv):
        aggr_csv = pd.read_csv(aggr_csv)
        sample_map = dict((str(i), n) for i,n in enumerate(aggr_csv['library_id']))
        samples = [sample_map[i] for i in barcodes_enum]
        data.obs['library_id'] = samples
  
    return data

def read_star(fn, args, **kw):
    mtx_dir = os.path.dirname(fn)
    dirname = os.path.dirname(mtx_dir)
    data = sc.readwrite._read_legacy_10x_mtx(mtx_dir, var_names='gene_ids', make_unique=True)
    data.var['gene_ids'] = list(data.var_names)
    sample_id = os.path.basename(dirname)
    data.obs['library_id'] = [sample_id] * data.obs.shape[0]
    barcodes = [b.split('-')[0] for b in data.obs.index]
    if len(barcodes) == len(set(barcodes)):
        data.obs_names = barcodes
    
    return data


def read_alevin(fn, args, **kw):
    avn_dir = os.path.dirname(fn)
    dirname = os.path.dirname(avn_dir)
    if fn.endswith('.gz'):
        df = alevin_parser.read_quants_bin(input_dir)
    else:
        df = alevin_parser.read_quants_csv(input_dir)
    row = {'row_names': df.index.values.astype(str)}
    col = {'col_names': np.array(df.columns, dtype=str)}
    data = AnnData(df.values, row, col, dtype=np.float32)
    data.var['gene_ids'] = list(data.var_names)
    sample_id = os.path.basename(dirname)
    data.obs['library_id'] = [sample_id] * data.obs.shape[0]
    
    return data
    
def read_umitools(fn, **kw):
    raise NotImplementedError

READERS = {'cellranger_aggr': read_cellranger_aggr, 'cellranger': read_cellranger, 'star': read_star, 'umitools': read_umitools, 'alevin': read_alevin}
        
if __name__ == '__main__':
    args = parser.parse_args()
    
    reader = READERS.get(args.input_format.lower())
    if reader is None:
        raise ValueError('{} is not a supported input format'.format(args.input_format))
    for fn in args.input:
        if not os.path.exists(fn):
            raise IOError('file does not exist! {}'.format(fn))
    n_input = len(args.input)
    if n_input > 1:
        assert(args.input_format != 'cellranger_aggr')
            
    if args.sample_info is not None:
        sample_info = pd.read_csv(args.sample_info, sep='\t')
        if not 'Sample_ID' in sample_info.columns:
            raise ValueError('sample_sheet needs a column called `Sample_ID`')
        sample_info.index = sample_info['Sample_ID']
        if args.batch is not None:
            batch_categories = sample_info[batch].astype('category')
    else:
        sample_info = None
        if args.batch is not None:
            raise ValueError('cannot use option `batch` when option `--sample-info` not used')
        batch_categories = None
        
    if args.feature_info is not None:
        feature_info = pd.read_csv(args.feature_info, sep='\t')
        if not 'gene_id' in feature_info.columns:
            raise ValueError('feature_info needs a column called `gene_id`')
        
        feature_info.index = feature_info['gene_id']
    else:
        feature_info = None
    
    data_list = []
    for i, fn in enumerate(args.input):
        fn = os.path.abspath(fn)
        data = reader(fn, args)
        data_list.append(data)

    if len(data_list) > 1:
        if args.normalize == 'mapped':
            data_list = downsample_gemgroup(data_list)

    data = data_list.pop(0)
    if len(data_list) > 0:
        data = data.concatenate(*data_list, batch_categories=batch_categories)
        # clean up feature info (assumes inner join)
        # fixme: maybe outer join support is what we want ?
        if any(i.endswith('-0') for i in data.var.columns):
            keep = [i for i in data.var.columns if i.endswith('-0')]
            data.var  = data.var.loc[:,keep]
            data.var.columns = [i.split('-')[0] for i in data.var.columns]
    if sample_info:
        lib_ids = set(data.obs['library_id'])
        for l in lib_ids:
            if l not in sample_info.index:
                raise ValueError('Library `{}` not present in sample_info'.format(l))
        obs = sample_info.loc[data.obs['library_id'],:]
        obs.index = data.obs.index.copy()
        data.obs = data.obs.merge(obs, how='left', on='barcode', copy=False)

    if not args.no_zero_cell_rm:
        keep = data.X.sum(1).A.squeeze() > 0
        data = data[keep,:]
        keep = data.X.sum(0).A.squeeze() > 0 
        #keep = (data.X != 0).any(axis=0)
        data = data[:,keep]
        
    if feature_info:
        data.var = data.var.merge(feature_info, how='left', on='gene_id', copy=False)
        
    if 'gene_symbols' in data.var.columns:
        mito_genes = data.var.gene_symbols.str.lower().str.startswith('mt-')
        data.obs['fraction_mito'] = np.sum(data[:, mito_genes].X, axis=1).A1 / np.sum(data.X, axis=1).A1
    data.obs['n_counts'] = data.X.sum(axis=1).A1

    if args.verbose:
        print(data)
        print(data.X.A.sum())
        
    if args.output_format == 'anndata':
        data.write(args.outfile)
    elif args.output_format == 'loom':
        data.write_loom(args.outfile)
    elif args.output_format == 'csvs':
        data.write_csvs(args.outpfile)
    else:
        raise ValueError("Unknown output format: {}".format(args.output_format))
