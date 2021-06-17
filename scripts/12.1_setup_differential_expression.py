#!/usr/bin/env python3

'''
Preparing the featurecounts matrix for DESeq2
'''

import os
import pandas as pd
import argparse as ap

# **********************************************************************************************************************

parser = ap.ArgumentParser()
parser.add_argument('-o', help='Path to output diffexp folder', required=True, type=str)
parser.add_argument('-i', help='Path to featurecounts input table', required=True, type=str)
args = parser.parse_args()

df_fc = pd.read_csv(args.i, sep='\t', skiprows=1, index_col=0)
df_fc.index = ['|'.join([x, df_fc['gene_name'].loc[x]]) if isinstance(df_fc['gene_name'].loc[x], str) else
               '|'.join([x, 'Unknown_%s' % x]) for x in df_fc.index]  # Eg. Phantom genes do not have a name
df_fc.drop(['Chr', 'Start', 'End', 'Strand', 'Length', 'gene_name'], axis=1, inplace=True)
df_fc.columns = [os.path.basename(sample).split('.')[0] for sample in df_fc.columns]
if not os.path.exists(args.o):
    os.mkdir(args.o)
featurecounts_out = os.path.join(args.o, 'feature_count.tsv')
df_fc.to_csv(featurecounts_out, sep='\t', index_label='GeneID')
