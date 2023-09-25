#!/usr/bin/env python3
'''
Scans the list of off-targets from CRISPRoff and re-evaluates their binding potential allowing for bulges.
'''

# **********************************************************************************************************************

import argparse as ap
from Bio import Seq
import pandas as pd
import subprocess
import os
from RisearchEval import get_best_match
import pickle

DICT_STRAND = {0: '-', 1: '+'}


# **********************************************************************************************************************

# Notes:
# crisproff coordinates are already zero-based exclusive-end, like twoBitToFa
# Note that in the link to UCSC of the crisproff output instead there are 1-based end-inclusive coordinates

def get_seq(chr: str, start: int, end: int, dict_chroms_lengths: dict, two_bit_genome: str):
    if chr in dict_chroms_lengths:
        if start >= 0 and end < dict_chroms_lengths[chr]:
            seq = subprocess.check_output(
                ['twoBitToFa', '-seq=%s' % chr, '-start=%i' % start, '-end=%i' % end,
                 '-noMask', two_bit_genome, 'stdout']).decode('utf-8')
            return ''.join(seq.split('\n')[1:])
    return ''


def search_bs(ser: pd.Series, strand: str, gRNA_length: int, extend_binding: int,
              dict_chroms_lengths: dict, cut_pos: int, two_bit_genome: str):
    cut_site = int(ser['Start cut region']) # No -1, already 0-based
    if strand == '+':
        return Seq.Seq(get_seq(chr=ser['Off-target candidate chr.'],
                               start=cut_site - cut_pos - gRNA_length - extend_binding,
                               end=cut_site - cut_pos,
                               dict_chroms_lengths=dict_chroms_lengths,
                               two_bit_genome=two_bit_genome)).reverse_complement()  # end pos includes binding site
    else:
        return Seq.Seq(
            get_seq(chr=ser['Off-target candidate chr.'],
                    start=cut_site + cut_pos,
                    end=cut_site + cut_pos + gRNA_length + extend_binding,
                    dict_chroms_lengths=dict_chroms_lengths,
                    two_bit_genome=two_bit_genome))

def is_intended_edit(variant: pd.Series, edits_lst: list):
    for x in range(variant['Start cut region'], variant['End cut region']+1):#start cut region=end cut region, 1-based from CRISPRoff, on-target coords are 1-based.
        for edit in edits_lst:
            if '%s\t%i\t%i' % (variant['Off-target candidate chr.'], x-1, x) == edit: #lifted edits are saved as [edit pos - 1, edit pos]
                return True
    return False

def extract_off_target_candidates(path_input_table: str, gRNA_length: int, cut_pos: int, twobit='', extend_binding=0):
    df = pd.read_csv(path_input_table, sep='\t', header=None, comment='#')  # do not include lines starting with '#'
    df.columns = ['Off-target candidate chr.', 'Start cut region', 'End cut region', 'gRNA (DNA)',
                  'Off-targeting score', 'Strand off-target PAM', 'PAM', 'GENOMIC FEATURES']
    df['PAM'] = df['PAM'].apply(lambda x: x.upper())
    df['gRNA (DNA)'] = df['PAM'].apply(lambda x: x.upper())
    to_remove_ontarget = []
    with open(args.on_target_pos, 'r') as edit_pos:
        edits = [x.rstrip('\n') for x in edit_pos.readlines()]
    if args.crisproff and extend_binding == 0:
        for id, line in df.iterrows():
            if is_intended_edit(variant=line, edits_lst=edits):
                to_remove_ontarget.append(id)
                continue
            df.at[id, 'DNA_BS'] = str(Seq.Seq(line['gRNA (DNA)']).reverse_complement())
            df.at[id, 'OFF-TARGET+CONTEXT'] = str(Seq.Seq(line['gRNA (DNA)']))
            df.at[id, 'ENDONUCLEASE_BS'] = line['PAM'][
                                           args.binding_sites_distance:] if args.binding_sites_distance > 0 else line[
                'PAM']
    else:  # this will include the bulges, if required
        dict_chroms_lengths = pickle.load(open(args.dict_genome_info, 'rb'))
        for id, line in df.iterrows():
            if is_intended_edit(variant=line, edits_lst=edits):
                to_remove_ontarget.append(id)
                continue
            dnabs = str(search_bs(ser=line, strand=line['Strand off-target PAM'],
                                  cut_pos=cut_pos, dict_chroms_lengths=dict_chroms_lengths,
                                  two_bit_genome=twobit,
                                  gRNA_length=gRNA_length,
                                  extend_binding=extend_binding))
            df.at[id, 'DNA_BS'] = dnabs
            df.at[id, 'ENDONUCLEASE_BS'] = line['PAM'][
                                           args.binding_sites_distance:] if args.binding_sites_distance > 0 else line[
                'PAM']
            df.at[id, 'OFF-TARGET+CONTEXT'] = str(Seq.Seq(dnabs).reverse_complement())
    return df.drop(to_remove_ontarget)


def get_pam_ratios(pams: list, ratios: list):
    if len(pams) != len(ratios):
        parser.error('The number of input PAM ratios must be equal to the number of PAMs')
    return {pams[i]: ratios[i] for i in range(len(pams))}


# **********************************************************************************************************************

parser = ap.ArgumentParser()
parser.add_argument('-i', '--input',
                    help='Table containing predicted off-target coordinates intersected with gene differential expression data',
                    type=str, required=True)
parser.add_argument('-g', '--gRNA_sequence', help='gRNA sequence', type=str, required=True)
parser.add_argument('-bd', '--binding_sites_distance',
                    help='Number of nucleotides between the end of the gRNA and the beginning of the endonuclease binding site',
                    type=int, required=True)
parser.add_argument('-s', '--seed_size', help='Size of the seed region', type=int, required=True)
parser.add_argument('-t', '--two_bit_genome', help='Genome in twobit format', type=str, default='')
parser.add_argument('--extend_binding',
                    help='Try to find a better binding site by extending the search on the DNA at the PAM-distal region',
                    type=int, default=0)
parser.add_argument('-o', '--output', help='Path to the folder to be used for output', type=str, default='tmp')
parser.add_argument('-cp', '--cut_position',
                    help='Number of nucleotides from the 3p end of the predicted off-target coordinates to the cleavage position. '
                         'The coordinates should refer to the forward DNA strand.',
                    type=int, required=True, default=-3)
parser.add_argument('--crisproff', help='Flag. Use if the input is the crisproff tsv output table',
                    action='store_true', required=False)
parser.add_argument('-d', '--dict_genome_info', help='Path to dictionary of chromosome length', type=str, required=True)
parser.add_argument('-rfo', '--rna_fold_output', help='Path to output of RNAfold', type=str, required=True)
parser.add_argument('-b', '--binding_sites', help='Sequence of the endonuclease binding site(s)',
                    type=str, nargs='+', required=True)
parser.add_argument('-br', '--binding_sites_ratios', help='List of weights for binding site(s)',
                    type=float, nargs='+', required=True)
parser.add_argument('-p', '--on_target_pos', help='Position(s) where on-target edits are expected (bed file).',
                    type=str, default="")

args = parser.parse_args()
if args.extend_binding != 0 and args.two_bit_genome == '':
    parser.error('--bulges_size requires --two_bit_genome')
if os.stat(args.input).st_size == 0:
    pd.DataFrame().to_csv(os.path.join(args.output, 'EvaluatedExpressionOffTargets.bed'), sep='\t', header=None,
                          index=None)
    pd.DataFrame().to_csv(os.path.join(args.output, 'EvaluatedExpressionOffTargets.tsv'), sep='\t', index_label='ID')
else:
    path_scripts = os.path.dirname(os.path.realpath(__file__))
    df = extract_off_target_candidates(path_input_table=args.input,
                                       twobit=args.two_bit_genome, extend_binding=args.extend_binding,
                                       cut_pos=args.cut_position,
                                       gRNA_length=len(args.gRNA_sequence))
    for elem in df[df['DNA_BS'] == ''].index:
        print(
            'Candidate off-target with genomic coordinates: %s:%i-%i(%s) not found in the given reference (2bit) genome.\n'
            % (df.iloc[elem]['Off-target candidate chr.'],
               df.iloc[elem]['Start cut region'],
               df.iloc[elem]['End cut region'],
               df.iloc[elem]['Strand off-target PAM'])
        )
    df = df[df['DNA_BS'] != '']
    pamratios = get_pam_ratios(pams=args.binding_sites, ratios=args.binding_sites_ratios)
    df = get_best_match(df=df, gRNA=args.gRNA_sequence, pamratios=pamratios, rna_fold_out=args.rna_fold_output,
                        seed_len=args.seed_size)
    os.path.join(args.output, '.tmp')
    df.to_csv(os.path.join(args.output, 'EvaluatedExpressionOffTargets.tsv'), sep='\t', index_label='ID')

    df_bed = df[['Off-target candidate chr.']].copy()
    df_bed['POS'] = df['Start cut region']
    df_bed['POS2'] = df['End cut region']
    df_bed['Name'] = list(df.index)
    df_bed.to_csv(os.path.join(args.output, 'EvaluatedExpressionOffTargets.bed'), sep='\t', index=None, header=None)
