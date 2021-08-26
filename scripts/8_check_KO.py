#!/usr/bin/env python
'''
Verifies the status of on-target knockouts by analyzing the differential expression data
'''

# **********************************************************************************************************************

import pandas as pd
import argparse as ap

FLOAT_THRESHOLD_NO_EXPRESSION = 3

# **********************************************************************************************************************

parser = ap.ArgumentParser()
parser.add_argument('-d', '--deseq2', help='Path to DESeq2 ouput file', type=str, required=True)
parser.add_argument('-o', '--output', help='Path to output file', type=str, required=True)
parser.add_argument('-KO', help='Transcript(s) IDs targeted', type=str, nargs='+', required=True)
args = parser.parse_args()

df_report = pd.DataFrame()
df = pd.read_csv(args.deseq2, index_col=0, sep='\t')
samples = df.columns[6:]
df.index = [x.split('|')[0] for x in df.index]
for x in args.KO:
    assert x in df.index, 'Gene ID %s was not found in the given annotations' % x
    for i in samples:
        if df.at[x, i] > FLOAT_THRESHOLD_NO_EXPRESSION:
            if df.at[x, 'log2FoldChange'] > 0 and df.at[x, 'padj'] < 0.05:
                df_report.at[i, '%s' % x] = 'Est. reads: %.2f\nUpregulated in Edited\n%.2f L2FC, %.2E P-adj' % (
                    df.at[x, i],df.at[x, 'log2FoldChange'], df.at[x, 'padj'])
            elif df.at[x, 'log2FoldChange'] < 0 and df.at[x, 'padj'] < 0.05:
                df_report.at[i, '%s' % x] = 'Est. reads: %.2f\nDownregulated in Edited\n%.2f L2FC, %.2E P-adj' % (
                    df.at[x, i],df.at[x, 'log2FoldChange'], df.at[x, 'padj'])
            else:
                df_report.at[i, '%s' % x] = 'Est. reads: %.2f\nNo expression difference\n%.2f L2FC, %.2E P-adj' % (
                    df.at[x, i],df.at[x, 'log2FoldChange'], df.at[x, 'padj'])
        else:
            df_report.at[i, '%s' % x] = 'Absent %.2f est. reads' % df.at[x, i]
with pd.ExcelWriter(args.output, engine='openpyxl') as writer:
    df_report.to_excel(writer, sheet_name='Summary KO', index_label='Sample')
