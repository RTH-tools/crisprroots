#!/usr/bin/env python3

'''
Generate RNAfold output for the gRNA
'''

# **********************************************************************************************************************

import argparse as ap
from subprocess import Popen

# **********************************************************************************************************************

parser = ap.ArgumentParser()
parser.add_argument('-of', '--out_RNAfold', help='Path to RNAfold output file', type=str, required=True)
parser.add_argument('-og', '--out_gRNA', help='Path to gRNA fasta file', type=str, required=True)
parser.add_argument('-g', '--gRNA', help='gRNA sequence', type=str, required=True)
args = parser.parse_args()

with open(args.out_gRNA, 'w') as out:
    out.write('>gRNA\n%s\n' % (args.gRNA))
with open(args.out_RNAfold, 'w') as out:
    p = Popen(['RNAfold', '-i', args.out_gRNA], stdout=out)
    p.wait()
print('RNAfold: completed')
