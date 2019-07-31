import sys
import os
import argparse

import scvelo as sv
import scanpy as sc

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('infiles', help='input filenames (velocyto .loom files)', nargs='+')
parser.add_argument('-o', '--outfile', help="output filename")
parser.add_argument('-f', '--format', choices=['anndata', 'loom', 'csvs'], default='anndata', help="output file format")
parser.add_argument('--sample-sheet', help="samplesheet filename")


def aggr_velocyto_loom(infiles):
    """read and concatenate velocyto loom files
    """
    infiles = sorted(infiles)
    fn = infiles.pop(0)
    adata = sc.read_loom(fn, var_names='Accession')
    adata.var.index.name = 'gene_ids'
    adata.var.rename(columns={'Gene': 'gene_symbols'})
    sv.utils.clean_obs_names(adata)
    sample_id = os.path.splitext(os.path.basename(fn))[0]
    adata.obs['library_id'] = sample_id
    adata.obs.index.name = 'barcodes'
    adata.var.index.name = 'gene_ids'

    if len(infiles) > 0:
        for fn in infiles:
            _adata = sc.read_loom(fn, var_names='Accession')
            sv.utils.clean_obs_names(_adata)
            sample_id = os.path.splitext(os.path.basename(fn))[0]
            _adata.obs['library_id'] = [sample_id] * _adata.obs.shape[0]
            _adata.obs.index.name = 'barcodes'
            _adata.var.index.name = 'gene_id'
            _adata.var.rename(columns={'Gene': 'gene_id'})
            adata = adata.concatenate(_adata, batch_key='library_id')
            
    return adata

if __name__ == '__main__':
    args = parser.parse_args()
    
