#!/usr/bin/env python3

'''
This script takes in input DESeq2 output and outputs 4 files containins:
- strongly downregulated genes
- strongly upregulated genes
- genes with insufficient expression
- expressed-genes with no extreme differential expression

'''

# **********************************************************************************************************************

import pandas as pd
import argparse as ap
import os

EXPRESSED = 10
L2FC = 0.5
PADJ = 0.01

# **********************************************************************************************************************

parser = ap.ArgumentParser()
parser.add_argument('-o', '--output', help="Path to output folder", type=str, required=True)
parser.add_argument('-r', '--report', help="Path to report folder", type=str, required=True)
parser.add_argument('-KO', help="Gene(s) IDs targeted", type=str, nargs='*', default=[])
parser.add_argument('-d', '--DESeq2_out', help="")
parser.add_argument('-p', '--phenotypes', help='Path to phenotypes table', type=str)
args = parser.parse_args()

no_expression = open(os.path.join(args.output, 'genes_NotDENotExpressed.tsv'), 'w')
flagged_de_up = open(os.path.join(args.output, 'genes_DEUp.tsv'), 'w')
flagged_de_down = open(os.path.join(args.output, 'genes_DEDown.tsv'), 'w')
expressed = open(os.path.join(args.output, 'genes_NotDEExpressed.tsv'), 'w')
for x in [no_expression, flagged_de_up, flagged_de_down, expressed]:
    x.write('#Minimum mean reads required for expression = %f\n'
            '#Min L2FC for differentially expressed genes = %f\n'
            'Gene\t'
            'L2FC\t'
            'Basemean exp. Original\t'
            'Basemean exp. Edited\n' % (EXPRESSED, L2FC))

df_p = pd.read_csv(args.phenotypes, sep='\t', index_col='Sample_ID')  # dataframe of phenotypes
original = list(df_p[df_p['Condition'] == 'Original'].index)
edited = list(df_p[df_p['Condition'] == 'Edited'].index)

dds = pd.read_csv(args.DESeq2_out, sep='\t', index_col=0)
dds['mean original'] = dds[original].mean(axis=1)
dds['mean edited'] = dds[edited].mean(axis=1)
dds.to_excel(os.path.join(args.report, 'DESeq2_differential_expression.xlsx'), sheet_name='DESeq2')

for id, line in dds.iterrows():
    if id.split('|')[0] not in args.KO:
        if line['baseMean'] < EXPRESSED:
            no_expression.write(id + '\t' + 'NA' + '\t' + str(line['mean original']) + '\t' + str(line[
                                                                                                      'mean edited']) + '\n')  # possible DE analysis false negatives (CRISPR edits impedes expression after editing, of a gene that was previously non expressed)
        else:
            text = id + '\t' + str(line['log2FoldChange']) + '\t' + str(line['mean original']) + '\t' + str(
                line['mean edited']) + '\n'
            if line['padj'] < PADJ:
                if line['log2FoldChange'] > L2FC:
                    flagged_de_up.write(text)
                elif line['log2FoldChange'] < -L2FC:
                    flagged_de_down.write(text)
                else:
                    expressed.write(text)
            else:
                expressed.write(text)
    else:
        print('Removing %s from off-target candidates output as it is one of the intended KO genes\n' % (
            id))
