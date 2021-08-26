#!/usr/bin/env python3

'''
Uses the outputs of the variant-based screening and the expression-based screening to generate a unique final report.
Flags for overlapping repeatmasked regions and dbSNP entries are added.
Off-targets are split into categories: minor, medium type 1-2-3, critical.
A PASS flag is added to signal off-targets that fulfill the requirements for minimum binding energy and max mismatches.
In a previous version, the maximum number of mismatches/unused bases was also considered. Code for this is commented.
'''

# **********************************************************************************************************************

import pandas as pd
import argparse as ap
from pandas import errors
import os


# **********************************************************************************************************************

def count_mismatches(x: str, l=None):
    x = x.split('_')
    n_m = 0  # number of mismatches
    n_b = 0  # number of bulges
    n_x = 0  # number unused
    for qt, pattern in enumerate(x):  # qt=0 if query, =1 if target
        if l:
            pattern = pattern[-l:]
            count_m = pattern.count('M')
            count_b = pattern.count('B')
            count_x = pattern.count('X')
        else:
            len_binding = len(pattern.lstrip('X'))  # how long the actual binding is
            count_m = pattern.count('M')
            count_b = pattern.count('B')
            count_x = len_gRNA - len_binding if len_binding < len_gRNA else 0  # to account for bindings that are shorter than the gRNA-DNA interaction
        n_m = count_m if n_m == 0 else n_m
        n_b = n_b + count_b
        n_x = count_x if qt == 1 else n_x
    return n_m + n_b + n_x


def map_risk(x):
    if x['N. SEED MISMATCH'] == 0:  # and x['N. MISMATCH/UNUSED'] <= args.total_mismatch_tolerance
        return 'CRITICAL'
    elif x[
        'N. SEED MISMATCH'] <= args.seed_mismatch_tolerance:  # and x['N. MISMATCH/UNUSED'] <= args.total_mismatch_tolerance
        return 'MAJOR TYPE 2'
    else:
        return 'MINOR'


def is_ontarget_gene(x):
    if (type(x['GENOMIC FEATURES']) == str) and any([k in x['GENOMIC FEATURES'] for k in args.ko]):
        return True  # the variant is located at the gene targeted for knockout but not at the expected position
    else:
        return False


def map_risk_de(x):
    if x['N. SEED MISMATCH'] == 0 and 'DEDown' in x[
        'DE_EVENTS']:  # and x['N. MISMATCH/UNUSED'] <= args.total_mismatch_tolerance
        return 'CRITICAL'
    elif x['N. SEED MISMATCH'] == 0 and 'DEUp' in x[
        'DE_EVENTS']:  # and x['N. MISMATCH/UNUSED'] <= args.total_mismatch_tolerance
        return 'MAJOR TYPE 1'
    elif x['N. SEED MISMATCH'] == 0 and x['N. MISMATCH/UNUSED'] == 0 and ' Exp' not in x['DE_EVENTS']:
        return 'MAJOR TYPE 3'
    elif x['N. SEED MISMATCH'] <= args.seed_mismatch_tolerance and ('DEUp' in x['DE_EVENTS'] or 'DEDown' in x[
        'DE_EVENTS']):  # and x['N. MISMATCH/UNUSED'] <= args.total_mismatch_tolerance
        return 'MAJOR TYPE 2'
    # elif x['N. SEED MISMATCH']<= args.seed_mismatch_tolerance and x['N. MISMATCH/UNUSED'] <= args.total_mismatch_tolerance:
    #    return 'MINOR'
    else:
        return 'MINOR'


def get_DE_event(events):
    str = ''
    for x in events.split('|'):
        if 'DE=DEUp' in x:
            str = str + 'DEUp|'
        if 'DE=DEDown' in x:
            str = str + 'DEDown|'
        if 'DE=NotDEExpressed' in x:
            str = str + 'NotDE, Exp|'
        if 'DE=NotDENotExpressed' in x:
            str = str + 'NotDE, NotExp|'
    return str[:-1]


def add_flags(df: pd.DataFrame, fr: str, fs: str):
    for x, col in [(fr, 'Repeatmask'), (fs, 'dbSNP')]:
        df[col] = ''
        if x == 'no':
            continue
        try:
            df_flags = pd.read_csv(x, sep='\t', header=None)
        except errors.EmptyDataError:
            continue
        for id in df_flags[3]:
            if id in df.index:
                df.at[id, col] = 'Yes'
    return df


def add_genes_overlap(df: pd.DataFrame, genes_overlap: str):
    def get_overlapping_genes(df):
        info = []
        for id, df_gene in df.groupby('gene_id'):
            info.append('Source=%s;GeneID=%s;GeneName=%s;Type=%s' % (
                df_gene[5].unique()[0], id, df_gene['gene_name'].unique()[0], ','.join(df_gene[6].unique())))
        return '|'.join(info)

    def get_feature(info, field):
        info = info.replace(' ', '').replace('"', '').split(';')
        for i in info:
            if field in i:
                return i.replace(field, '')

    df['Overlapping genes'] = ''
    try:
        df_overlaps = pd.read_csv(genes_overlap, sep='\t', header=None)
    except errors.EmptyDataError:
        return pd.DataFrame()
    else:
        df_overlaps = df_overlaps[df_overlaps[12] != '.']
        df_overlaps['gene_id'] = df_overlaps[12].apply(lambda x: get_feature(x, 'gene_id'))
        df_overlaps['gene_name'] = df_overlaps[12].apply(lambda x: get_feature(x, 'gene_name'))
        for id, df_features in df_overlaps.groupby(3):
            if id in df.index:
                df.at[id, 'GENOMIC FEATURES'] = get_overlapping_genes(df_features)
        return df


def filter_pass(x: pd.Series, int_max_mm_seed: int, float_max_deltagb: float):
    if x['N. SEED MISMATCH'] > int_max_mm_seed or x[
        'DeltaG_B'] > float_max_deltagb:  # x['N. MISMATCH/UNUSED'] > int_max_mm_total or
        return 'NO'
    else:
        return 'YES'


def from_0_to_1_based_coords(x):
    chr, startend, sign = x.split(':')
    start, end = startend.split('-')
    return ':'.join([chr, '-'.join([str(int(start) + 1), str(int(end) + 1)]), sign])


# **********************************************************************************************************************

parser = ap.ArgumentParser()
parser.add_argument('-oft', help='Table with potential off-targets from expression-based off-target screening',
                    type=str, metavar='<file>',
                    required=True)
parser.add_argument('-ofv', help='Table with potential off-targets from variant-based off-target screening', type=str,
                    metavar='<file>', required=True)
parser.add_argument('-s', '--seed', help='Length of the seed region', type=int, metavar='<int>', default=0)
parser.add_argument('-g', '--gRNA_sequence', help='gRNA sequence', metavar='<str>', type=str, required=True)
parser.add_argument('--extend_binding', help='Max size allowed for extend_binding', metavar='<int>', type=int,
                    default=0)
parser.add_argument('-cp', '--cut_position',
                    help='Number of nucleotides from the 3p end of the predicted off-target coordinates to the cleavage position. '
                         'The coordinates should refer to the forward DNA strand.', metavar='<int>', type=int,
                    required=True, default=-3)
parser.add_argument('-sm', '--seed_mismatch_tolerance', help='Max mismatches allowed in seed', type=int,
                    metavar='<int>', default=1)
# parser.add_argument('-tm', '--total_mismatch_tolerance', help='Max number of mismatches in gRNA-DNA binding', type=int, metavar='<int>', default=1)
parser.add_argument('-deltagb', '--float_max_deltagb',
                    help='Maximum weighted binding energy allowed in gRNA-DNA binding', type=float)
parser.add_argument('-frt', '--flag_repeatmask_transcripts',
                    help='Path to bed file with off-targets from transcriptome expression analysis which are flagged because of overlap with repeatmasked regions',
                    type=str, required=True)
parser.add_argument('-frv', '--flag_repeatmask_variants',
                    help='Path to bed file with off-targets from variant-based analysis which are flagged because of overlap with repeatmasked regions',
                    type=str, required=True)
parser.add_argument('-fsv', '--flag_dbsnp_variants',
                    help='Path to bed file with off-targets from variant-based analysis which are flagged because of overlap with SNPs reported in dbSNP',
                    type=str, required=True)
parser.add_argument('-vg', '--vars_genes', help='bed file of coordinates of genes intersecting variants', type=str,
                    required=True)
parser.add_argument('-o', help='Path to output file', required=True, type=str)
parser.add_argument('-ko', help='Transcript(s) IDs targeted', required=True, type=str, nargs='+')
parser.add_argument('-type', help='Type of edit: KI (knockin) or KO (knockout)', required=True, type=str,
                    choices=['KI', 'KO'])
args = parser.parse_args()

len_gRNA = len(args.gRNA_sequence)
if os.stat(args.ofv).st_size == 0:
    df_v = pd.DataFrame(
        columns=['COORDINATES (1-based inclusive)', 'OFF-TARGET+CONTEXT', 'PAM', 'BINDING_STRUCT(Q-T)', 'EVENT', 'TLOD',
                 'DeltaG_B', 'dbSNP', 'Repeatmask', 'GENOMIC FEATURES'])
else:
    df_v = pd.read_csv(args.ofv, sep='\t', index_col='VARIANT ID')
    df_v = add_flags(df=df_v, fr=args.flag_repeatmask_variants, fs=args.flag_dbsnp_variants)
    df_v = add_genes_overlap(df=df_v, genes_overlap=args.vars_genes)
if os.stat(args.oft).st_size == 0:
    df_t = pd.DataFrame(
        columns=['COORDINATES (1-based inclusive)', 'OFF-TARGET+CONTEXT', 'PAM', 'BINDING_STRUCT(Q-T)', 'DE_EVENTS',
                 'DeltaG_B', 'Repeatmask', 'GENOMIC FEATURES'])
else:
    df_t = pd.read_csv(args.oft, sep='\t', index_col='ID')
    df_t = add_flags(df=df_t, fr=args.flag_repeatmask_transcripts, fs='no')
if len(df_v) > 0:
    INVERT_STRAND = {'+': '-', '-': '+'}
    df_v['COORDINATES (0-based inclusive)'] = df_v['CHROM'] + ':' + df_v['START_BINDING'].astype('int').astype(
        'str') + '-' + df_v['END_BINDING'].apply(lambda x: str(int(x) - 1)) + ':' + df_v['STRAND_BINDING'].apply(
        lambda x: INVERT_STRAND[x])  # the strand is inverted to reflect the strand of the PAM, not of the binding site;  # coordinates are 0 based inclusive (-1 because on UCSC the last position is included [start:end]=[1,2] is 2 nts, not 1)
    df_v['COORDINATES (1-based inclusive)'] = df_v['COORDINATES (0-based inclusive)'].apply(
        lambda x: from_0_to_1_based_coords(x))
    df_v['CUT_SITE'] = df_v['CUT_SITE'].astype('int')
    df_v = df_v[['COORDINATES (1-based inclusive)', 'OFF-TARGET+CONTEXT', 'PAM', 'BINDING_STRUCT(Q-T)', 'EVENT', 'TLOD',
                 'DeltaG_B', 'dbSNP', 'Repeatmask', 'GENOMIC FEATURES']]
    df_v['N. SEED MISMATCH'] = df_v['BINDING_STRUCT(Q-T)'].apply(lambda x: count_mismatches(x=x, l=args.seed))
    df_v['N. MISMATCH/UNUSED'] = df_v['BINDING_STRUCT(Q-T)'].apply(lambda x: count_mismatches(x=x))
    df_v['RISK'] = df_v.apply(lambda x: map_risk(x), axis=1)
    if args.type == 'KO':
        df_v = df_v[df_v.apply(lambda x: not is_ontarget_gene(x), axis=1)].copy()
    df_v['PASS'] = df_v[['N. SEED MISMATCH', 'DeltaG_B']].apply(
        lambda x: filter_pass(x=x, int_max_mm_seed=args.seed_mismatch_tolerance,
                              float_max_deltagb=args.float_max_deltagb), axis=1)
    df_v['DeltaG_B'] = df_v['DeltaG_B'].apply(lambda x: '%.2f' % x)
    df_v.sort_values(["RISK", "DeltaG_B"], ascending=(True, True), inplace=True)
    df_v.drop('TLOD', axis = 1, inplace=True)
if len(df_t) > 0:
    def adjust_coords(ser_pos_strand, info):  # Coordinates are 0 based
        if info == 'Start':  # Start position on + stran
            if ser_pos_strand['Strand off-target PAM'] == '+':
                return str(ser_pos_strand['Start cut region'] - args.cut_position - len(
                    args.gRNA_sequence) - args.extend_binding)
            else:
                return str(ser_pos_strand['Start cut region'] + args.cut_position)
        elif info == 'End':
            if ser_pos_strand['Strand off-target PAM'] == '+':
                return str(ser_pos_strand[
                               'End cut region'] - args.cut_position - 1)  # -1 because on UCSC the last position is included [start:end]=[1,2] is 2 nts, not 1
            else:
                return str(ser_pos_strand['End cut region'] + args.cut_position + len(
                    args.gRNA_sequence) + args.extend_binding - 1)  # -1 because on UCSC the last position is included [start:end]=[1,2] is 2 nts, not 1
        else:
            exit(1)


    df_t['COORDINATES (0-based inclusive)'] = df_t['Off-target candidate chr.'] + ':' + df_t[
        ['Start cut region', 'Strand off-target PAM']].apply(
        lambda ser_pos_strand: adjust_coords(ser_pos_strand, 'Start'), axis=1) + '-' + df_t[
                                                  ['End cut region', 'Strand off-target PAM']].apply(
        lambda ser_pos_strand: adjust_coords(ser_pos_strand, 'End'), axis=1) + ':' + df_t['Strand off-target PAM']
    df_t['COORDINATES (1-based inclusive)'] = df_t['COORDINATES (0-based inclusive)'].apply(
        lambda x: from_0_to_1_based_coords(x))
    df_t['DE_EVENTS'] = df_t['GENOMIC FEATURES'].apply(lambda x: get_DE_event(x))
    df_t = df_t[
        ['COORDINATES (1-based inclusive)', 'OFF-TARGET+CONTEXT', 'PAM', 'BINDING_STRUCT(Q-T)', 'DE_EVENTS', 'DeltaG_B',
         'Repeatmask', 'GENOMIC FEATURES']]
    df_t['N. SEED MISMATCH'] = df_t['BINDING_STRUCT(Q-T)'].apply(lambda x: count_mismatches(x=x, l=args.seed))
    df_t['N. MISMATCH/UNUSED'] = df_t['BINDING_STRUCT(Q-T)'].apply(lambda x: count_mismatches(x=x))
    df_t['RISK'] = df_t.apply(lambda x: map_risk_de(x), axis=1)
    df_t['PASS'] = df_t[['N. SEED MISMATCH', 'DeltaG_B']].apply(
        lambda x: filter_pass(x=x, int_max_mm_seed=args.seed_mismatch_tolerance,
                              float_max_deltagb=args.float_max_deltagb), axis=1)
    df_t['DeltaG_B'] = df_t['DeltaG_B'].apply(lambda x: '%.2f' % x)
    df_t.sort_values(['RISK', 'DeltaG_B'], ascending=(True, True), inplace=True)
with pd.ExcelWriter(args.o, engine='openpyxl') as writer:
    df_t.to_excel(writer, sheet_name='Expression-based', index_label='ID')
    df_v.to_excel(writer, sheet_name='Variants-based', index_label='ID')
