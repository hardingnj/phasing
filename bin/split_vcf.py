#! /usr/bin/python

__author__ = 'Nicholas Harding'

import argparse
import gzip
import re
import os

parser = argparse.ArgumentParser(description='Tool to filter a vcf file')

parser.add_argument('input', help='input compressed vcf file')
parser.add_argument('--out', action='store', dest='out_dir',
                    default=os.getcwd())
parser.add_argument('--start', action='store', dest='start', type=int)
parser.add_argument('--end', action='store', dest='end', type=int)

args = parser.parse_args()
comment = re.compile('^#')

stem = os.path.split(args.input)[-1].split('.')[0]
output_fn = args.out_dir + "/{0}_{1}_{2}.vcf.gz".format(stem,
                                                        args.start, args.end)

input_fh = gzip.open(args.input, 'r')
output_fh = gzip.open(output_fn, 'w-')

print(args.start, args.end)
for l in input_fh:
    if comment.match(l):
        output_fh.write(l)
    else:
        arr = l[:20].split("\t")
        pos = int(arr[1])
        if pos < args.start:
            continue
        elif pos > args.end:
            break
        else:
            output_fh.write(l)

output_fh.close()
input_fh.close()
