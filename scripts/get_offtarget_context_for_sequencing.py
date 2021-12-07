#!/usr/bin/env python3

'''
This is a bonus script, not part of the pipeline. Can be used to obtain the sequences around a possible off-target in
the variant-aware genome. Useful to design the off-target sanger sequencing validation.
'''

# **********************************************************************************************************************

import pandas as pd
import argparse as ap
import os
import subprocess

# **********************************************************************************************************************

parser = ap.ArgumentParser()
parser.add_argument('-ot', '--off_target', help='Path to table with off-targets', required=True, type=str)
parser.add_argument('-r', '--risk', help='Risk class of the off-targets', required=True, type=str,
                    choices=['CRITICAL', 'MAJOR TYPE 1', 'MAJOR TYPE 2', 'MAJOR TYPE 3', 'MINOR'])
parser.add_argument('-t', '--type', help='Only vairant-based, expression-based or both types of off-targets',
                    required=True, type=str,
                    choices=['Expression-based', 'Variant-based', 'Both'])
parser.add_argument('-n', '--n_context', help='Number of nucleotides of left-right context', required=True, type=int)
parser.add_argument('-g', '--genome', help='Genome in fasta format', type=str, required=True)
parser.add_argument('-rk', '--filter_repeatmask', help='Remove repeatmasked off-targets', action='store_true')
parser.add_argument('-ksnp', '--filter_known_snp', help='Remove off-targets related to known snps', action='store_true')
parser.add_argument('-top', help='Output only the top off-targets sorted by lowest to highest binding energy', type=int)
parser.add_argument('-outf', help='Output folder', required=True, type=str)
parser.add_argument('-s', '--strandness',
                    help='Use to have strand-specific sequences, otherwhise everything will be in 5p-3p?',
                    action='store_true')
parser.add_argument('-os', '--only_stats', help='Only print off-target stats', action='store_true')
args = parser.parse_args()

if args.type == 'Variant-based':
    df = pd.read_excel(args.off_target, sheet_name=args.type)[
        ['COORDINATES (1-based inclusive)', 'RISK', 'PASS', 'Repeatmask', 'dbSNP', 'DeltaG_B']]
elif args.type == 'Expression-based':
    df = pd.read_excel(args.off_target, sheet_name=args.type)[
        ['COORDINATES (1-based inclusive)', 'RISK', 'PASS', 'Repeatmask', 'DeltaG_B']]
else:
    df1 = pd.read_excel(args.off_target, sheet_name='expression-based')[
        ['COORDINATES (1-based inclusive)', 'RISK', 'PASS', 'Repeatmask', 'DeltaG_B']]
    df2 = pd.read_excel(args.off_target, sheet_name='variant-based')[
        ['COORDINATES (1-based inclusive)', 'RISK', 'PASS', 'Repeatmask', 'dbSNP', 'DeltaG_B']]
    df = df1.append(df2)
print('there are %i off-targets candidates (OC)' % len(df))
print('Removing %i OC with low risk or not passing thresholds' % len(
    df[(df['RISK'] != args.risk) | (df['PASS'] != 'YES')]))
df = df[(df['RISK'] == args.risk) & (df['PASS'] == 'YES')]
print('There are %i off-targets candidates (OC)' % len(df))
if args.filter_repeatmask:
    print('Removing %i OC overlapping repeatmasked region' % len(df[~df['Repeatmask'].isna()]))
    df = df[df['Repeatmask'].isna()]
    print('There are %i off-targets candidates (OC)' % len(df))
if args.filter_known_snp:
    if 'dbSNP' in df.columns:
        print('Removing %i OC overlapping repeatmasked region' % len(df[~df['dbSNP'].isna()]))
        df = df[df['dbSNP'].isna()]
        print('There are %i off-targets candidates (OC)' % len(df))
    else:
        print('Skipping filter dbSNP, not applicable to expression-based off-targets only')
print('Remaining OC fulfilling criteria=%i' % len(df))
if args.top:
    df.sort_values('DeltaG_B', inplace=True)
    print('Keep only top %i off-targets' % (args.top))
    df = df.head(n=args.top)
    print('There are %i off-targets candidates (OC)' % len(df))
print(df)
if not os.path.exists(args.outf):
    os.mkdir(args.outf)
bedfile = os.path.join(args.outf,
                       '%s_%s.bed' % (os.path.basename(args.off_target).split('.')[0], args.type))
with open(bedfile, 'w') as out:
    for i, l in df.iterrows():
        chr = l['COORDINATES (1-based inclusive)'].split(':')[0]
        start, end = l['COORDINATES (1-based inclusive)'].split(':')[1].split('-')
        strand = l['COORDINATES (1-based inclusive)'].split(':')[2]
        out.write('%s\t%i\t%i\t%s\t.\t%s\n' % (
            chr, int(start) - args.n_context - 1, int(end) + args.n_context, l['COORDINATES (1-based inclusive)'],
            strand))
fastafile = os.path.join(args.outf, '%s_%s.fa' % (
    os.path.basename(args.off_target).split('.')[0], args.type))
if not args.only_stats:
    if args.strandness:
        subprocess.check_output(
            ['bedtools', 'getfasta', '-fi', args.genome, '-bed', bedfile, '-fo', fastafile, '-name', '-s'])
    else:
        subprocess.check_output(
            ['bedtools', 'getfasta', '-fi', args.genome, '-bed', bedfile, '-fo', fastafile, '-name'])
