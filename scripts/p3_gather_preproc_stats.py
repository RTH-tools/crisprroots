#!/usr/bin/env python3

'''
Summarizes mapping statistics from STAR output
'''
import pandas as pd
import argparse as ap
import os
import glob

parser = ap.ArgumentParser()
parser.add_argument('-i',
                    help='Path to the stats folder . The script expect this folder to contain one tsv file with suffix stats.tsv for each analysis whose stats have been reported. This tsv folder should contain as raws the names of the samples (with column index "Sample") and as columns the statistics computed. In case of stats on forwards-reverse reads the name of the sample should contain R1-R2.',
                    type=str)
parser.add_argument('-o',
                    help='Path to the output summary file. Can be more than one file in case of paired reads results in searate raws (one output for stats on R1-R2 of each sample, another for stats on each samples).',
                    type=str)
args = parser.parse_args()

df_stats_samples = pd.DataFrame()
df_stats_R1_R2 = pd.DataFrame()

for path_to_stats_table in glob.glob(os.path.join(args.i, '*')):
    stats_filename = os.path.basename(path_to_stats_table)
    df_stats = pd.read_csv(path_to_stats_table, index_col='Sample', sep='\t').sort_index()
    print('Report %s contains stats on paired end reads' % stats_filename)
    df_stats.columns = ['_'.join([col, stats_filename[:-4]]) for col in df_stats.columns]
    df_stats.drop([col for col in df_stats.columns if 'total_sequences' not in col], axis=1, inplace=True)
    df_stats_R1_R2 = pd.concat([df_stats_R1_R2, df_stats], axis=1)
df_stats_R1_R2.to_excel(args.o, index_label='Sample', sheet_name='Preprocessing stats')
