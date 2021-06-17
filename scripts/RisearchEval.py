#!/usr/bin/env python3

'''
Given a gRNA and a target evaluates possible bindings between the two.
The bindings are forced to start at the 3' of the target and end at the 3' of the gRNA and at the 5' of the target.
The target is sliced by removing one nucleotide at the 3' end at a time such that all suboptimals are evaluated.
Among all suboptimals it returns the one with minimum binding energy.
'''

# **********************************************************************************************************************

import pandas as pd
import subprocess
import os
from Bio import Seq
import numpy as np

SUGIMOTO_MATRIX = 'su95'


def get_best_match(df: pd.DataFrame, gRNA: str, rna_fold_out: str, pamratios: dict, seed_len: int):
    rnafold = get_rnafold(path_RNA_fold=rna_fold_out)
    for bs, line in df.iterrows():
        # get best RIsearch interaction terminating at the PAM
        min_dgb = np.nan
        min_struct = ''
        if line['DNA_BS'] != '':
            for i in range(0, len(line['DNA_BS']) - 1, ):
                target = str(line['DNA_BS'])[:len(str(line['DNA_BS'])) - i]
                if len(target) < seed_len:
                    continue
                risearch_out = subprocess.check_output(
                    ['RIsearch1', '-Q', gRNA, '-T', target, '-m', SUGIMOTO_MATRIX, '-f', '5000', '-w',
                     'CRISPR_20nt_5p_3p']).decode()
                risearch_out = risearch_out.split('\n')
                rnadna = float(risearch_out[-2].split(':')[1].split(' ')[1])
                tmp_bs = '_'.join([risearch_out[-3].split('\t')[0], 'X' * i + risearch_out[-3].split('\t')[
                    1]])  # Structure of the binding in the format QueryRNAStruct_TargetDNAStruct
                dna_dna = get_DNA_DNA_energy(Seq.Seq(line['DNA_BS']).reverse_complement(),
                                             len(risearch_out[-3].split('\t')[1].replace('X', '')))
                tmp_dgb = get_crisproff_deltagb(rnadna=rnadna, dnadna=dna_dna, rnafold=rnafold,
                                                PAM=line['ENDONUCLEASE_BS'], pamratios=pamratios)
                if np.isnan(min_dgb) or tmp_dgb < min_dgb:
                    min_dgb = tmp_dgb
                    min_struct = tmp_bs
            df.at[bs, 'BINDING_STRUCT(Q-T)'] = min_struct
            df.at[bs, 'DeltaG_B'] = min_dgb
    return df


def get_rnafold(path_RNA_fold: str):
    with open(path_RNA_fold, 'r') as infile:
        return float(infile.readlines()[-1].split(' ')[-1].replace('(', '').replace(')', ''))


def get_DNA_DNA_energy(seq, len_DNA_interaction):
    # the content of this function is taken from the CRISPRoff pipeline (Alkan et al. Genome Biology 2018)
    RI_REV_NT_MAP = {'-': '', 'a': 'T', 'A': 'T', 'c': 'G', 'C': 'G', 'g': 'C', 'G': 'C',
                     't': 'A', 'T': 'A', 'u': 'A', 'U': 'A', 'n': 'N', 'N': 'N'}

    RI_DNA_DNA_NN = {'AA': {'TT': -1.00}, 'TT': {'AA': -1.00}, 'AT': {'TA': -0.88}, 'TA': {'AT': -0.58},
                     'CA': {'GT': -1.45}, 'TG': {'AC': -1.45}, 'GT': {'CA': -1.44}, 'AC': {'TG': -1.44},
                     'CT': {'GA': -1.28}, 'AG': {'TC': -1.28}, 'GA': {'CT': -1.30}, 'TC': {'AG': -1.30},
                     'CG': {'GC': -2.17}, 'GC': {'CG': -2.24}, 'GG': {'CC': -1.84}, 'CC': {'GG': -1.84},
                     'NA': {'NT': 0.05}, 'NT': {'NA': 0.05}, 'NC': {'NG': 0.05}, 'NG': {'NC': 0.05},
                     'AN': {'TN': 4.35}, 'TN': {'AN': 4.35}, 'GN': {'CN': 3.95}, 'CN': {'GN': 3.95},
                     'NN': {'NN': 0.66}}
    if len_DNA_interaction == 0:
        return 0.0
    else:
        seq = seq.upper()[-len_DNA_interaction:]
        energy = [0.0] * len(seq)
        for i in range(1, len(seq)):
            energy[i] = float(RI_DNA_DNA_NN[seq[i - 1] + seq[i]][RI_REV_NT_MAP[seq[i - 1]] + RI_REV_NT_MAP[seq[i]]])
        return sum(energy)


def get_crisproff_deltagb(rnadna: float, dnadna: float, rnafold: float, PAM: str, pamratios: dict):
    if PAM in pamratios:
        return pamratios[PAM] * (rnadna - dnadna - rnafold)
    else:
        PAM_tmp = PAM[:]
        while len(PAM_tmp) > 0:
            PAM_tmp = PAM_tmp[:-1]
            if PAM_tmp in pamratios:
                return pamratios[PAM_tmp] * (rnadna - dnadna - rnafold)
        exit('PAM binding site %s not found in the dictionary of PAMs and PAM ratios' % PAM)
