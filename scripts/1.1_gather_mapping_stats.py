#!/usr/bin/env python3

'''
Summarizes mapping statistics from STAR output
'''

# **********************************************************************************************************************

import pandas as pd
import argparse as ap
import os

# **********************************************************************************************************************

def parse_star(star_out: str):
    ser = pd.Series(name=os.path.basename(star_out))
    with open(star_out, 'r') as s:
        s = s.readlines()
        s = [l.split('|') for l in s]
    for l in s[:-1]:
        if len(l) > 1:
            ser[l[0].rstrip().lstrip()] = l[1].rstrip().lstrip()
    return ser


# **********************************************************************************************************************

parser = ap.ArgumentParser()
parser.add_argument('-i', help='Path to stats output from STAR', type=str, nargs='+')
parser.add_argument('-o', help='Path to the output summary file', type=str)
args = parser.parse_args()

df_stats_samples = pd.DataFrame()
for path_to_stats_table in args.i:
    ser_stats = parse_star(path_to_stats_table)
    df_stats_samples = df_stats_samples.join(ser_stats, how='outer')
df_stats_samples.to_excel(args.o, sheet_name='Mapping stats')
