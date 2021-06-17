#!/usr/bin/env python3
'''
Candidate off-targets located at the exact same position and strand are collapsed in one line, in which the information regarding all overlapping genes is included.
'''

# **********************************************************************************************************************

import pandas as pd
import argparse as ap
import os

# **********************************************************************************************************************

parser = ap.ArgumentParser()
parser.add_argument('-b', '--bed_intersect',
                    help='File containing the genomic intersection between off-target candidates and differentially expressed genes in BED format',
                    type=str, required=True)
parser.add_argument('-o', '--out_file',
                    help='Collapsed output file, containing the collapsed genomic intersection between off-target candidates and differentially expressed genes in BED format. ',
                    type=str, required=True)
args = parser.parse_args()

if os.stat(args.bed_intersect).st_size == 0:
    pd.DataFrame().to_csv(args.out_file, sep='\t', header=None, index=None)
else:
    df = pd.read_table(args.bed_intersect, sep='\t', header=None)
    df[1] = df[1].astype(str)
    df[2] = df[2].astype(str)
    df['coordinates'] = df[[0, 1, 2, 5]].apply(lambda x: '|'.join(x), axis=1)
    collapsed_df = pd.DataFrame()
    for id, df_tmp in df.groupby('coordinates'):
        df_tmp = df_tmp.copy()
        df_tmp[10] = ' | '.join([x if x != '.' else 'Intergenic' for x in df_tmp[10]])
        collapsed_df = collapsed_df.append(df_tmp.iloc[0].copy())
    collapsed_df.drop('coordinates', axis=1, inplace=True)
    collapsed_df[[0, 1, 2, 3, 4, 5, 6, 10]].to_csv(args.out_file, sep='\t', header=None, index=None)
