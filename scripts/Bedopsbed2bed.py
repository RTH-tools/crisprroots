#!/usr/bin/env python3

'''
This script takes as input a file in bed format from Bedops, which is the result of transforming a VCF into BED format.
The script saved the first 3 columns (chr, start, end) and adds two columns.
4th column (score) info about the event from mutect2 (ref, alt, pos on reference, tlod)
5th column (name) the length of the event: negative for deletions, positive for insertions, 0 for SNVs.
'''

# **********************************************************************************************************************

import pandas as pd
import argparse as ap
import os
from io import StringIO

# **********************************************************************************************************************

def get_eventlength(ser: pd.Series):
    ref = ser[5]
    lst_alt = ser[6].split(',')
    length = 0
    for x in lst_alt:
        if len(x) - len(ref) < length:
            length = len(x) - len(ref)  # deletions
        if length == 0 and len(x) - len(ref) > 0:
            length = len(x) - len(ref)  # insertions (if already deletion is present, keep deletion)
    return length


def extract_tlod(x: str):
    lst_x = x.split(';')
    for x in lst_x:
        if 'TLOD=' in x:
            return float(x.rstrip().split('=')[1].split(',')[0])


# **********************************************************************************************************************

parser = ap.ArgumentParser()
parser.add_argument('-i', help='Path to input file in Bedops bed format', type=str, required=True)
parser.add_argument('-o', help='Path to output file in bed format', type=str, required=True)
args = parser.parse_args()


N_INFO_ROWS = 10

with open(args.i, 'r') as infile:
    lines = infile.readlines()
    header_rows = [line for line in lines if line.startswith('_header')]
    df_rows = [line for line in lines if not line.startswith('_header')]
sorted_samples = [x.rstrip() for x in header_rows[-1].split('\t')][N_INFO_ROWS+2:]
if len(df_rows)==0:
    pd.DataFrame(columns=['#CHR', 'START', 'END', 'EVENT', 'EVENTLENGTH']).to_csv(args.o)
else:
    df = pd.read_csv(StringIO('\n'.join(df_rows)), sep='\t', header=None)
    df['EVENT'] = df[[0, 1, 2, 5, 6, 8]].apply(lambda x: 'pos=%s:%i-%i|ref=%s|alt=%s|tlod=%.2f|' % (
        x[0], x[1], x[2], x[5], ';'.join(x[6].split(',')), extract_tlod(x=x[8])), axis=1)
    df['EVENTLENGTH'] = df[[5, 6]].apply(lambda x: get_eventlength(ser=x), axis=1)
    for i,s in enumerate(sorted_samples):
        df[s] = df[N_INFO_ROWS+i]
        df[s] = df[s].apply(lambda x: ':'.join(x.split(':')[:2]))
    for s in sorted_samples:
        df['EVENT'] = df['EVENT'] + s + ':' + df[s] + ';'
    df = df[[0,1,2,'EVENT', 'EVENTLENGTH']]
    df.columns = ['#CHR', 'START', 'END', 'EVENT', 'EVENTLENGTH']
    df.to_csv(args.o, sep='\t', index=None)

