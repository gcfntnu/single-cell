#!/usr/bin/env python

import sys
import argparse

import scanpy.api as sc

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('infile', help='input file', type=argparse.FileType('r'))
parser.add_argument('-o', '--outfile', help="output file, if empty stdout is used", default=sys.stdout, type=argparse.FileType('w'))
    

args = parser.parse_args(sys.argv)

sc.read_
