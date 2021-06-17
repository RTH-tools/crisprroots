#!/usr/bin/env python3

'''
Remove off-targets that are on a hemizygous chromosome and overlap an expressed gene.
If any variation occurred, it should be captured in the variant-based off-target screening.
'''

# **********************************************************************************************************************

import pandas as pd
import argparse as ap
import os

EXPRESSED = 10


# **********************************************************************************************************************

def get_exp(info: str, field: str):
    lst_info = info.split(';')
    max = 0.0
    for info in lst_info:
        if field in info:
            curr = float(info.split('=')[1].replace(' ', '').split('|')[0])
            if curr > max:
                max = curr
            return max


# **********************************************************************************************************************

parser = ap.ArgumentParser()
parser.add_argument('-c', '--collapsed_coords', help='Path to collapsed off-target coordinates', type=str,
                    required=True)
parser.add_argument('-o', '--output', help='Path to output file', type=str, required=True)
parser.add_argument('-s', '--single', help='List of single chromosomes (eg. chrY)', type=str, nargs='+', required=True)
args = parser.parse_args()

if os.stat(args.collapsed_coords).st_size == 0:
    pd.DataFrame().to_csv(args.output, sep='\t', header=None, index=None)
else:
    df_c = pd.read_csv(args.collapsed_coords, sep='\t', header=None)
    df_c_filter = df_c.copy()
    df_c_filter.columns = ['CHROM', 'SPOS', 'EPOS', 'SEQ', 'SCORE', 'STRAND', 'PAM', 'INFO']
    df_c_filter = df_c_filter[df_c_filter['INFO'].apply(
        lambda x: ('Type=Gene' in x) and (get_exp(x, 'OR=') > EXPRESSED) and (
                get_exp(x, 'ER=') > EXPRESSED))]  # The filter is applied only on expressed genes
    to_be_filtered = []
    for chr in args.single:
        df_c_tmp = df_c_filter[df_c_filter['CHROM'] == chr]
        to_be_filtered.extend(list(df_c_tmp.index))
    df_c = df_c.drop(to_be_filtered, axis=0)
    df_c.to_csv(args.output, sep='\t', index=None, header=None)
