#!/usr/bin/env python3

'''
Given a set of variants, the script identifies potential Cas9 binding sites and evaluates the gRNA-DNA binding pattern
with RIsearchSubopt.
Out of all possible off-target bindings linked to a variant, it preserves the one with lowest hybridization energy.
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


def get_seq_with_delimiters(chr: str, start: int, end: int, dict_chroms_lengths: dict):
    if chr in dict_chroms_lengths:
        seq = subprocess.check_output(
            ['twoBitToFa', '-seq=%s' % chr, '-start=%i' % max(0, start), '-end=%i' % min(dict_chroms_lengths[chr], end),
             '-noMask', args.two_bit_genome, 'stdout']).decode('utf-8')
        return ''.join(seq.split('\n')[1:]), max(0, start), min(dict_chroms_lengths[chr], end)


def search_bs_given_seq(ser: pd.Series, back: int, front: int, dict_chroms_lengths: dict, seq_region: str, sstart: int,
                        send: int):
    df_bss = pd.DataFrame()
    for i in range(back, front + 1, 1):
        for len_bs in binding_sites.keys():
            cut_site = ser['START'] + 1 + i  # have to do +1 because in twobit end position is excluded
            start_p = cut_site - args.cut_position - len(args.gRNA_sequence) - args.extend_binding  # PAM on plus strand
            end_p = cut_site - args.cut_position + args.binding_sites_distance + len_bs  # PAM on plus strand
            start_m = cut_site + args.cut_position - args.binding_sites_distance - len_bs  # PAM on minus strand
            end_m = cut_site + args.cut_position + len(args.gRNA_sequence) + args.extend_binding  # PAM on minus strand
            seq_gRNA_and_bs_p = Seq.Seq(seq_region[start_p - sstart:end_p - send])
            seq_gRNA_and_bs_m = Seq.Seq(seq_region[start_m - sstart:end_m - send]).reverse_complement()
            for strand, seq in enumerate([seq_gRNA_and_bs_p, seq_gRNA_and_bs_m]):
                if str(seq[len(seq) - len_bs:len(seq)]) in binding_sites[len_bs]:  # if there is a binding site
                    df_bss = df_bss.append(pd.Series({'CUT_SITE': cut_site,
                                                      'STRAND_BINDING': DICT_STRAND[strand],
                                                      'START_BINDING': start_p if strand == 0 else start_m + args.binding_sites_distance + len_bs,
                                                      'END_BINDING': end_p - args.binding_sites_distance - len_bs if strand == 0 else end_m,
                                                      'ENDONUCLEASE_BS': str(seq[len(seq) - len_bs:len(seq)]),
                                                      'PAM': str(seq[
                                                                 len(seq) - len_bs - args.binding_sites_distance:len(
                                                                     seq)]),
                                                      'OFF-TARGET+CONTEXT': str(
                                                          seq[:len(seq) - len_bs - args.binding_sites_distance]),
                                                      'DNA_BS': str(seq[:len(
                                                          seq) - len_bs - args.binding_sites_distance].reverse_complement())}),
                                           ignore_index=True)  # sr = search region, bs= binding site, E_bs = endonuclease binding site
    return df_bss


def evaluate_bs(df_bss: pd.DataFrame, pamratios: dict):
    if len(df_bss) > 0:
        return get_best_match(df=df_bss, gRNA=args.gRNA_sequence, pamratios=pamratios,
                              rna_fold_out=args.rna_fold_output, seed_len=args.seed_size)
    else:
        return pd.DataFrame()


def search_evaluate_bs(ser: pd.Series, dict_chroms_lengths: dict, pamratios: dict):
    seq_region, sstart, send = get_seq_with_delimiters(chr=ser['CHROM'],
                                                       start=ser['START'] - args.expand_search + min((
                                                               -args.cut_position - len(
                                                           args.gRNA_sequence) - args.extend_binding),
                                                           (
                                                                   args.cut_position - args.binding_sites_distance - max(
                                                               binding_sites.keys()))),
                                                       end=ser['START'] + 1 + abs(
                                                           min(0, ser['EVENTLENGTH'])) + args.expand_search + max((
                                                               -args.cut_position + args.binding_sites_distance + max(
                                                           binding_sites.keys())),
                                                           (
                                                                   args.cut_position + len(
                                                               args.gRNA_sequence) + args.extend_binding)) + 1,
                                                       dict_chroms_lengths=dict_chroms_lengths)
    if ser['EVENTLENGTH'] == 0:
        return evaluate_bs(search_bs_given_seq(ser=ser, back=-1 - args.expand_search, front=0 + args.expand_search,
                                               dict_chroms_lengths=dict_chroms_lengths, seq_region=seq_region,
                                               sstart=sstart, send=send), pamratios=pamratios)
    if ser['EVENTLENGTH'] > 0:
        return evaluate_bs(search_bs_given_seq(ser=ser, back=0 - args.expand_search, front=0 + args.expand_search,
                                               dict_chroms_lengths=dict_chroms_lengths, seq_region=seq_region,
                                               sstart=sstart, send=send), pamratios=pamratios)
    if ser['EVENTLENGTH'] < 0:
        return evaluate_bs(search_bs_given_seq(ser=ser, back=0 - args.expand_search,
                                               front=abs(ser['EVENTLENGTH']) + args.expand_search,
                                               dict_chroms_lengths=dict_chroms_lengths, seq_region=seq_region,
                                               sstart=sstart, send=send), pamratios=pamratios)


def get_chroms_lengths():
    return pickle.load(open(args.dict_genome_info, 'rb'))


def get_pam_ratios(pams: list, ratios: list):
    if len(pams) != len(ratios):
        parser.error('The number of input PAM ratios must be equal to the number of PAMs')
    return {pams[i]: ratios[i] for i in range(len(pams))}


def is_intended_edit(variant: pd.Series, edits_file: str):
    with open(edits_file, 'r') as edit_pos:
        edits = edit_pos.readlines()
    for x in range(variant['START'], variant['START'] + 1 + abs(min(0, variant['EVENTLENGTH'])) + args.expand_search):
        for edit in edits:
            if '%s\t%i\t%i' % (variant['CHROM'], x, x + 1) == edit.rstrip('\n'):
                return True
    return False


# **********************************************************************************************************************


parser = ap.ArgumentParser()
parser.add_argument('-v', '--variants_table', help='table containing identified variants', metavar='<file>',
                    type=str, required=True)
parser.add_argument('-cp', '--cut_position',
                    help='Number of nucleotides between the end of the guide RNA molecule and the cleavage position',
                    metavar='<int>', type=int, required=True)
parser.add_argument('-g', '--gRNA_sequence', help='gRNA sequence', metavar='<str>', type=str, required=True)
parser.add_argument('-b', '--binding_sites', help='Sequence of the endonuclease binding site(s)', metavar='<list<str>>',
                    type=str, nargs='+', required=True)
parser.add_argument('-bd', '--binding_sites_distance',
                    help='Number of nucleotides between the end of the gRNA and the beginning of the endonuclease binding site.',
                    metavar='<int>', type=int, required=True)
parser.add_argument('-t', '--two_bit_genome', help='Genome in twobit format', metavar='<file>', type=str,
                    required=True)
parser.add_argument('-s', '--seed_size', help='Size of the seed region', metavar='<int>', type=int, required=True)
parser.add_argument('--extend_binding',
                    help='Try to find a better binding site by extending the search on the DNA at the PAM-distal region',
                    metavar='<int>', type=int, default=0)
parser.add_argument('-o', '--output', help='Path to the folder to be used for output', metavar='<folder>', type=str,
                    default='tmp', required=True)
parser.add_argument('--expand_search',
                    help='Expand the search region for binding sites +- N nucleotides, where N is the given parameter. This implies that the cut site position considered is to N nucleotides shifted compared to the position in which the actual variant is called. This is allowed because a repair outcome could result in the same nucleotide as it was before cleavage.',
                    type=int, default=2)
parser.add_argument('-p', '--on_target_pos', help='Position(s) where on-target edits are expected (bed file).',
                    type=str, default="")
parser.add_argument('-d', '--dict_genome_info', help='Path to dictionary of chromosome length', type=str)
parser.add_argument('-rfo', '--rna_fold_output', help='Path to output of RNAfold', type=str)
parser.add_argument('-br', '--binding_sites_ratios', help='List of weights for binding site(s)',
                    metavar='<list<float>>',
                    type=float, nargs='+', required=True)
args = parser.parse_args()

# create output folders
if not os.path.exists(args.output):
    os.mkdir(args.output)
pamratios = get_pam_ratios(pams=args.binding_sites, ratios=args.binding_sites_ratios)
# get path scripts
path_scripts = os.path.dirname(os.path.realpath(__file__))
# get reference genome info
dict_chroms_lengths = get_chroms_lengths()
# parse variants table (bcf)
if os.stat(args.variants_table).st_size == 0:
    df_out = pd.DataFrame(
        columns=['CHROM', 'START', 'END', 'DNA_BS', 'OFF-TARGET+CONTEXT', 'EVENT', 'PAM', 'CUT_SITE',
                 'STRAND_BINDING', 'START_BINDING', 'END_BINDING', 'BINDING_STRUCT(Q-T)', 'DeltaG_B', 'TLOD'])
    df_out.to_csv(os.path.join(args.output, 'EvaluatedVariantsOffTargets.tsv'), sep='\t', index_label='VARIANT ID')
    open(os.path.join(args.output, 'EvaluatedVariantsOffTargets.bed'), 'w').close()
else:
    df_v = pd.read_table(args.variants_table, sep='\t', header=None)
    df_v.columns = ['CHROM', 'START', 'END', 'EVENT', 'EVENTLENGTH']
    df_v['TLOD'] = df_v['EVENT'].apply(lambda x: x.split('tlod=')[1])
    df_v['EVENT'] = df_v['EVENT'].apply(lambda x: x[:x.find('|tlod=')])
    df_v['START'] = df_v['START'].astype(int)
    df_v['END'] = df_v['END'].astype(int)
    lst_engs = []
    df_out = pd.DataFrame()
    # create dictionary of binding sites lengths
    binding_sites = {}
    for x in args.binding_sites:  # split binding sites by their length
        if len(x) in binding_sites:
            binding_sites[len(x)].append(x)
        else:
            binding_sites[len(x)] = [x]
    for id, line in df_v.iterrows():
        # get binding sites and energies
        #print(id)
        if not is_intended_edit(variant=line, edits_file=args.on_target_pos):
            df_bss = search_evaluate_bs(ser=line, dict_chroms_lengths=dict_chroms_lengths, pamratios=pamratios)
            if len(df_bss) > 0:
                df_bss['DNA_BS'] = df_bss['DNA_BS'].apply(lambda x: x[::-1])
                # save best binding site (min. free energy)
                df_bss.sort_values('DeltaG_B', inplace=True)
                df_out = df_out.append(pd.Series(line.append(df_bss.iloc[0]), name=id))
                lst_engs.append(df_bss.iloc[0]['DeltaG_B'])
        else:
            print('Removing %s:%i-%i from off-target candidates output as it is one of the on-target mutations\n' % (
                line['CHROM'], line['START'], line['END']))
    # output best bindings at each position
    if len(df_out) > 0:
        df_out.sort_values('DeltaG_B', inplace=True)
        df_out = df_out[
            ['CHROM', 'START', 'END', 'DNA_BS', 'OFF-TARGET+CONTEXT', 'EVENT', 'PAM', 'CUT_SITE', 'STRAND_BINDING',
             'START_BINDING', 'END_BINDING', 'BINDING_STRUCT(Q-T)', 'DeltaG_B', 'TLOD']]
    else:
        df_out = pd.DataFrame(
            columns=['CHROM', 'START', 'END', 'DNA_BS', 'OFF-TARGET+CONTEXT', 'EVENT', 'PAM', 'CUT_SITE',
                     'STRAND_BINDING', 'START_BINDING', 'END_BINDING', 'BINDING_STRUCT(Q-T)', 'DeltaG_B', 'TLOD'])
    df_out['START'] = df_out['START'].astype(int)
    df_out['END'] = df_out['END'].astype(int)
    df_out.to_csv(os.path.join(args.output, 'EvaluatedVariantsOffTargets.tsv'), sep='\t', index_label='VARIANT ID')

    df_bed = df_out[['CHROM']].copy()
    df_bed['START'] = df_out['START']
    df_bed['END'] = df_out['END']
    df_bed['Name'] = list(df_out.index)
    df_bed.to_csv(os.path.join(args.output, 'EvaluatedVariantsOffTargets.bed'), sep='\t', index=None, header=None)
