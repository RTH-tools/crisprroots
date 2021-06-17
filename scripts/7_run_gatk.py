#!/usr/bin/env python3

'''
Executes GATK Mutect2
'''

# **********************************************************************************************************************

import os
import subprocess
import argparse as ap

# **********************************************************************************************************************

parser = ap.ArgumentParser()
parser.add_argument('-r', '--reference', help='Path to reference', required=True, type=str)
parser.add_argument('-bqst', help='Base quality score threshold', required=True, type=int)
parser.add_argument('-cd', help='Callable depth', required=True, type=int)
parser.add_argument('-mqb', help='Min base quality score', required=True, type=int)
parser.add_argument('-o', help='Path to vcf output', required=True, type=str)
parser.add_argument('-ed', '--edited', help='Path to mapping files of edited samples', required=True, type=str,
                    nargs='+')
parser.add_argument('-or', '--original', help='Path to mapping files of original samples', required=True, type=str,
                    nargs='+')
parser.add_argument('-O', help='Path to filtered vcf output', required=True, type=str)

args = parser.parse_args()
cmd1 = ['gatk', 'Mutect2', '--reference', args.reference, '--base-quality-score-threshold', str(args.bqst),
        '--callable-depth', str(args.cd), '--min-base-quality-score', str(args.mqb), '-O', args.o]
for i in args.edited:
    cmd1.extend(['-I', i])
for o in args.original:
    cmd1.extend(['-I', o, '-normal', os.path.basename(o).split('.')[0]])  # remove .bam

cmd2 = ['gatk', 'FilterMutectCalls', '-V', args.o, '-O', args.O, '-R', args.reference]

subprocess.run(cmd1)
subprocess.run(cmd2)
