#!/usr/bin/env python3

'''
Given a table of possible off targets in bed format this script uses the coordinates and strand information to
return the cut positions, also in bed format.
If the input is a CRISPRoff file, then the coordinates of the off-targets include their PAMs, therefore the PAM is removed.
If the input is a bed file with predictions of off-targets, this is not necessary.
Input bed format should contain the coordinates of the off-target (NOT including the PAM) in the following format:
Chromosome, Start position, End position, - , - , Strand
The strand should reflect the strand of the PAM site, and is not optional.
'''

# **********************************************************************************************************************

import pandas as pd
import argparse as ap
# import numpy as np
import os

CRISPRoff_PAM_size = 3


# **********************************************************************************************************************

def make_cut_position_bed(path_input: str, path_output: str, cut_position: int, gRNA: str, report_folder: str,
                          is_crisproff=False, webserver=False):
    if is_crisproff:
        if webserver:
            df = pd.read_table(path_input, sep='\t', comment='#', usecols=[0, 1, 2, 3, 4, 5])
        else:
            df = pd.read_table(path_input, sep='\t', comment='#', header=None, usecols=[0, 1, 2, 3, 4, 5])
        df.columns = ['Chr', 'StartPos', 'EndPos', 'OffTarget-seq', 'CRISPRoff', 'Strand']
        #########CRISPRoff has 0-based coordinates#############
        df.loc[df['Strand'] == '+', 'StartPos'] = df.loc[df[
                                                             'Strand'] == '+', 'EndPos'] - CRISPRoff_PAM_size + cut_position  # Need to modify StartPos before EndPos otherwise StartPos will be updated twice
        df.loc[df['Strand'] == '+', 'EndPos'] = df.loc[
                                                    df['Strand'] == '+', 'EndPos'] - CRISPRoff_PAM_size + cut_position
        df.loc[df['Strand'] == '-', 'EndPos'] = df.loc[
                                                    df['Strand'] == '-', 'StartPos'] + CRISPRoff_PAM_size - cut_position
        df.loc[df['Strand'] == '-', 'StartPos'] = df.loc[df[
                                                             'Strand'] == '-', 'StartPos'] + CRISPRoff_PAM_size - cut_position  # Need to modify StartPos after EndPos otherwise EndPos will be updated twice
        df['PAM'] = df['OffTarget-seq'].apply(lambda x: x[-3:])
        df['OffTarget-seq'] = df['OffTarget-seq'].apply(lambda x: x[:-3])
        df = df[df['CRISPRoff'] > 0]
    else:  # here the pam is never included.
        df = pd.read_table(path_input, sep='\t', comment='#', header=None)
        df.columns = ['Chr', 'StartPos', 'EndPos', 'OffTarget-seq', 'CRISPRoff', 'Strand']
        df.loc[df['Strand'] == '+', 'StartPos'] = df.loc[df[
                                                             'Strand'] == '+', 'EndPos'] + cut_position  # Need to modify StartPos before EndPos otherwise StartPos will be updated twice
        df.loc[df['Strand'] == '+', 'EndPos'] = df.loc[df['Strand'] == '+', 'EndPos'] + cut_position
        df.loc[df['Strand'] == '-', 'EndPos'] = df.loc[df['Strand'] == '-', 'StartPos'] - cut_position
        df.loc[df['Strand'] == '-', 'StartPos'] = df.loc[df[
                                                             'Strand'] == '-', 'StartPos'] - cut_position  # Need to modify StartPos after EndPos otherwise EndPos will be updated twice
    if len(df[df['OffTarget-seq'] == gRNA]) > 1:
        with open(os.path.join(report_folder, 'Alert.txt'), 'w') as out:
            out.write('The gRNA %s has %i targets in the variated genome.\n '
                      'This alert could be an artifact in case a mutation in the patient genome is affecting the PAM site.\n'
                      'In that case, please provide the gRNA with the PAM including the genomic variation to avoid this warning.\n'
                      'You can easily check if this is the case by looking into the output file of CRISPRoff.\n'
                      'If the gRNA+PAM you provided has position "NA" in the genome then look for the match in the ouput that contains your gRNA but a different PAM.' % (
                          gRNA, len(df[df['OffTarget-seq'] == gRNA])))
    else:
        df = df[df['OffTarget-seq'] != gRNA]
    df.dropna(inplace=True)
    df['StartPos'] = df['StartPos'].astype(int)
    df['EndPos'] = df['EndPos'].astype(int)
    df.to_csv(path_output, sep='\t', header=False, index=False)


# **********************************************************************************************************************


parser = ap.ArgumentParser()
parser.add_argument('-i', '--input',
                    help="Path to input file. The input file contains the coordinates of predicted off-targets in bed format. Coordinates are in 5'-3' direction and must not include the PAM site, unless flag CRISPRoff is specified",
                    type=str, required=True)
parser.add_argument('-o', '--output', help='Path to output file.', type=str, required=True)
parser.add_argument('-cp', '--cut_position',
                    help='Number of nucleotides from the 3p end of the predicted off-target coordinates to the cleavage position. '
                         'The coordinates should refer to the forward DNA strand.',
                    type=int, required=True)
parser.add_argument('--crisproff', help='Flag. Use if the input is the crisproff tsv output table', action='store_true',
                    required=False)
parser.add_argument('--webserver', help='Flag. Use if the input is from the CRISPRoff webserver', action='store_true',
                    required=False)
parser.add_argument('-g', help='gRNA sequence', type=str, required=True)
parser.add_argument('-f', help='Path to report folder', type=str, required=True)
args = parser.parse_args()

make_cut_position_bed(path_input=args.input, path_output=args.output, cut_position=args.cut_position,
                      is_crisproff=args.crisproff, gRNA=args.g, report_folder=args.f, webserver=args.webserver)
