import os
import argparse
import subprocess

parser = argparse.ArgumentParser(description='Create testdata from subsampling')
parser.add_argument('input', nargs='+', help='fastq input files')
parser.add_argument('-o', help='output dir', dest='output')
parser.add_argument('-n', dest='number', help='number of reads in output fastq file')

args = parser.parse_args()


if not os.path.exists(args.output):
    subprocess.call('mkdir','-p ', args.output)

for fn in args.input:
    out_fn = os.path.join(args.output, fn)
    subprocess.call('zcat', fn, '|', 'seqtk sample', '-s123', '-', args.number, '|', 'gzip', '>', out_fn)
    
    
