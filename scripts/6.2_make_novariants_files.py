#!/usr/bin/env python3

'''
Generates empty variant files for step 6.2 in case no variants to the reference need to be called
'''

# **********************************************************************************************************************

import os
import argparse as ap

# **********************************************************************************************************************

parser = ap.ArgumentParser()
parser.add_argument('-v', help='Path to intersected variants', type=str)
args = parser.parse_args()

if not os.path.exists(os.path.dirname(args.v)):
    os.makedirs(os.path.dirname(args.v))
with open(args.v, 'wb') as out:
    out.write(
        '##fileformat=VCFv4.\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE'.encode())
