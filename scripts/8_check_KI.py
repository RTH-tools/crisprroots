#!/usr/bin/env python3

'''
Verifies the status of on-target edits by analyzing the mapping pileup
'''

# **********************************************************************************************************************

import pandas as pd
import argparse as ap
import subprocess as sp
import os


# **********************************************************************************************************************

def str_get_pielup_stats(int_nreads: int, str_nt: str, mutation: str):
    str_nt_no_indel = []
    int_indels = 0
    i = 0
    while i < len(str_nt):
        nt = str_nt[i]
        if nt != '+' and nt != '-' and nt != '$' and nt != '^':
            str_nt_no_indel.append(nt)
        else:
            if nt == '+' or nt == '-':
                int_indels += 1
                int_skipped_poss = int(
                    str_nt[i + 1])  # this number will say how many positions we need to remove from the pileup
                i = i + int_skipped_poss + 1  # +1 for skipping the number of indels (eg. +4ATTA)
                str_nt_no_indel = str_nt_no_indel[
                                  :-1]  # remove lats character of str_nt_no_indel, because it is connected with the indel event
        i += 1
    int_mutations = len([nt for nt in str_nt_no_indel if nt.upper() == mutation])
    lst_other_mutations = ['A', 'T', 'G', 'C']  # Do not consider N (unknown nt in read) as a possible mutation
    if mutation != '<>' and mutation != 'N':
        lst_other_mutations.remove(mutation)
    int_other_events = len([nt for nt in str_nt_no_indel if nt.upper() in lst_other_mutations]) + int_indels
    int_skip = len([nt for nt in str_nt_no_indel if nt == '>' or nt == '<'])
    if mutation != '<>':
        stats = 'Ref=%i; %s=%i; Skip=%i; Other=%i; \n%s' % \
                (int_nreads - int_other_events - int_mutations - int_skip, mutation, int_mutations, int_skip,
                 int_other_events, str_nt)
    else:
        stats = 'Ref=%i; Skip=%i; \n%s' % \
                (int_nreads - int_other_events - int_mutations - int_skip, int_skip, str_nt)
    return stats


def set_comment(ser_mutations: pd.Series, intron: list, lst_splice_acc: list, lst_splice_don: list,
                lst_positions: list):  # note that position is the list of originally given positions, before adding splice acceptors/donors and introns
    comment = ''
    for i, pos in enumerate(
            lst_positions):  # lst positions does not contain additional pos for upstream/downstream splice sites
        if ser_mutations[pos] != 'No read mapping at position %s' % pos:
            lst_str_ps = ser_mutations[pos].split('; ')  # pileup statistics
            ref = int(lst_str_ps[0][lst_str_ps[0].find('=') + 1:])
            mut = int(lst_str_ps[1][lst_str_ps[1].find('=') + 1:])
            skip = int(lst_str_ps[2][lst_str_ps[2].find('=') + 1:])
            other = int(lst_str_ps[3][lst_str_ps[3].find('=') + 1:])
            tot_reads = ref + mut + skip + other
            if intron[i]:
                str_chr = pos[:pos.find(':') + 1]
                int_pos = int(pos[pos.find(':') + 1:])
                if lst_splice_acc[i]:
                    lst_str_f = ser_mutations[str_chr + str(int_pos + 1)].split('; ')  # pileup statistics forward 1 nt
                    skip_f = int(lst_str_f[1][lst_str_f[1].find('=') + 1:])
                    if (mut > 0 or skip_f > 0) and (
                            skip == skip_f):  # if some has terminated then we would have a split like in the ref
                        comment = comment + 'Homozygous mutation at %s.\n' % pos
                    elif (mut > 0 or skip_f > 0) and (skip > skip_f):  # some skip has terminated
                        comment = comment + 'Heterozygous mutation at %s.\n' % pos
                    elif (mut == 0 and skip_f == 0 and skip == tot_reads):
                        comment = comment + 'No mutation detected at %s.\n' % pos
                    else:
                        comment = comment + 'Could not assess mutation status at %s.\n' % pos
                elif lst_splice_don[i]:
                    lst_str_b = ser_mutations[str_chr + str(int_pos - 1)].split('; ')  # pileup statistics back 1 nt
                    skip_b = int(lst_str_b[1][lst_str_b[1].find('=') + 1:])
                    if (mut > 0 or skip_b > 0) and (
                            skip == skip_b):  # if the skip just started then the splice donor is still there
                        comment = comment + 'Homozygous mutation at %s.\n' % pos
                    elif (mut > 0 and skip_b > 0) and (skip > skip_b):  # some skip has just started
                        comment = comment + 'Heterozygous mutation at %s.\n' % pos
                    elif (mut == 0 and skip_b == 0 and skip == tot_reads):
                        comment = comment + 'No mutation detected at %s.\n' % pos
                    else:
                        comment = comment + 'Could not assess mutation status at %s.\n' % pos
                else:
                    lst_str_f = ser_mutations[str_chr + str(int_pos + 1)].split('; ')  # pileup statistics forward 1 nt
                    ref_f = int(lst_str_f[0][lst_str_f[0].find('=') + 1:])
                    skip_f = int(lst_str_f[1][lst_str_f[1].find('=') + 1:])
                    lst_str_b = ser_mutations[str_chr + str(int_pos - 1)].split('; ')  # pileup statistics back 1 nt
                    ref_b = int(lst_str_b[0][lst_str_b[0].find('=') + 1:])
                    skip_b = int(lst_str_b[1][lst_str_b[1].find('=') + 1:])
                    if (skip > 0 or skip_b > 0 or skip_f > 0) and (mut == 0 and ref_f == 0 and ref_b == 0):
                        comment = comment + 'No mutation is detected at %s.\n' % pos
                    elif ((skip_b > 0 and skip_f > 0 and skip > 0)
                          and (mut > 0 or ref_f > 0 or ref_b > 0)):
                        comment = comment + 'Heterozygous mutation at %s.\n' % pos
                    elif (skip == 0 or (skip > 0 and ((skip_b != skip) or (skip_f != skip)))):
                        comment = comment + 'Homozygous mutation at %s.\n' % pos
                    else:
                        comment = comment + 'Could not assess mutation status at %s.\n' % pos
                if other > 0:
                    comment = comment + 'Other events (skipping, mutations) also present at %s\n' % pos
            else:
                if mut == 0 and ref > 0:
                    comment = comment + 'No mutation is detected at %s.\n' % pos
                elif mut > 0 and ref > 0:
                    comment = comment + 'Heterozygous mutation at %s.\n' % pos
                elif ref == 0 and mut > 0:
                    comment = comment + 'Homozygous mutation at %s.\n' % pos
                else:
                    comment = comment + 'Could not assess mutation status at %s.\n' % pos
                if skip > 0 or other > 0:
                    comment = comment + 'Other events (skipping, mutations) also present at %s\n' % pos
        else:
            return 'No read mapping at position %s' % str_pos
    return comment[:-1]


# **********************************************************************************************************************

parser = ap.ArgumentParser()
parser.add_argument('-p', help='Position(s) to analyze (eg. chr1:2356981 chr1:2356985).', required=True, type=str,
                    nargs='+')
parser.add_argument('-m',
                    help='Mutation to search for, given in the same order as the corresponding positions specified in parameter -p',
                    required=True, type=str, nargs='+', choices=['A', 'T', 'G', 'C', 'N', 'a', 't', 'g', 'c', 'n'])
parser.add_argument('-r', help='Fasta file of the reference genome.', required=True, type=str)
parser.add_argument('-b', help='Alignment file(s) in bam format. The correspornding index should also be present.',
                    required=True, type=str, nargs='+')
parser.add_argument('-o', help='Output file', required=True, type=str)
parser.add_argument('-sa', '--splice_acceptor',
                    help='The mutation affects a spice acceptor. For each mutation, use value 1 as "yes", 0 as "no" (eg. 1 0 0)',
                    type=int, nargs='+', required=True)
parser.add_argument('-sd', '--splice_donor',
                    help='The mutation affects a spice donor. For each mutation, use value 1 as "yes", 0 as "no" (eg. 0 0 1)',
                    type=int, nargs='+', required=True)
parser.add_argument('-i', '--intron',
                    help='The mutation affects an intron. For each mutation, use value 1 as "yes", 0 as "no" (eg. 0 1 0 )',
                    type=int, nargs='+', required=True)
args = parser.parse_args()

if not len(args.splice_acceptor) == len(args.p) == len(args.splice_donor) == len(args.intron) == len(args.m):
    raise parser.error(
        'The number of given positions has to match the number of mutations, splice_acceptor, splice_donor and intron elements')
for x in range(len(args.p)):
    if (args.splice_acceptor[x] or args.splice_donor[x]) and not args.intron:
        raise parser.error('Mutations in --splice_acceptors or --splice_donors --intron to be turned on.')
    if args.splice_acceptor[x] and args.splice_donor[x]:
        raise parser.error('--splice_acceptor and --splice_donor are mutually exclusive.')
for x in args.splice_acceptor:
    if x not in [0, 1]:
        raise parser.error('--splice_acceptor argument must be a list of 0 and 1')
for x in args.splice_donor:
    if x not in [0, 1]:
        raise parser.error('--splice_donor argument must be a list of 0 and 1')
    for x in args.intron:
        if x not in [0, 1]:
            raise parser.error('--intron argument must be a list of 0 and 1')

list_str_positions = args.p[:]
list_str_mutations = args.m[:]
for i in range(len(args.p)):
    if args.intron[i]:
        str_chr = args.p[i][:args.p[i].find(':') + 1]
        int_pos = int(args.p[i][args.p[i].find(':') + 1:])
        if args.splice_acceptor[i]:
            list_str_positions.append(str_chr + str(int_pos + 1))
            list_str_mutations.append('<>')
        elif args.splice_donor[i]:
            list_str_positions.append(str_chr + str(int_pos + -1))
            list_str_mutations.append('<>')
        else:  # in case a splice site results from the mutation we need to check up-down
            list_str_positions.append(str_chr + str(int_pos - 1))
            list_str_mutations.append('<>')
            list_str_positions.append(str_chr + str(int_pos + 1))
            list_str_mutations.append('<>')

index_isolates = [os.path.basename(os.path.dirname(x)) for x in args.b]
df = pd.DataFrame(index=index_isolates, columns=list_str_positions + ['Status'])

for path_alignment in args.b:
    id = os.path.basename(os.path.dirname(path_alignment))
    for i, str_pos in enumerate(list_str_positions):
        try:
            cmd = ['samtools', 'mpileup', '--reference', args.r, '-r',
                   str_pos + '-' + str(int(str_pos[str_pos.find(':') + 1:])), path_alignment]
            lst_out = sp.check_output(cmd).decode('utf-8').split('\t')
            if len(lst_out) > 1:
                df.loc[id][str_pos] = str_get_pielup_stats(int_nreads=int(lst_out[3]), str_nt=lst_out[4],
                                                           mutation=list_str_mutations[i].upper())
            else:
                df.loc[id][str_pos] = 'No read mapping at position %s' % str_pos
        except:
            df.loc[id][str_pos] = 'No read mapping at position %s' % str_pos

    df.loc[id]['Status'] = set_comment(ser_mutations=df.loc[id], intron=args.intron,
                                       lst_splice_acc=args.splice_acceptor, lst_splice_don=args.splice_donor,
                                       lst_positions=args.p)
with pd.ExcelWriter(args.o, engine='openpyxl') as writer:
    df.to_excel(writer, sheet_name='Summary KI', index_label='Sample')
