#!/usr/bin/env python3
import pandas as pd
import argparse as ap
import gzip

parser=ap.ArgumentParser()
parser.add_argument('-i', help='input file with variants, gz format', type=str, required=True)
parser.add_argument('-o', help='filtered variants, in which * is substituted with N', type=str, required=True)
args=parser.parse_args()

df=pd.read_csv(args.i, comment='#', compression='gzip', sep='\t', header=None)
df[4]=df[4].apply(lambda x: x.replace('*','N'))
comment=''
with gzip.open(args.i,'r') as infile:
    lines=[x.decode() for x in infile.readlines()]
    i=0
    while lines[i].startswith('#'):
        comment+=lines[i]
        i+=1
mat=df.values
with open(args.o, 'w') as out:
    out.write(comment)
    for i in range(len(mat)):
        for j in range(len(mat[i])):
            out.write(str(mat[i][j]))
            out.write('\t')
        out.write('\n')