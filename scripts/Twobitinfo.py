#!/usr/bin/env python3

'''
Get chomosome lengths from two bit file with twoBitInfo function as a dictionary of chromosome-length pairs
'''

# **********************************************************************************************************************

import pickle
import subprocess
import sys

# **********************************************************************************************************************

dict_chroms_lengths = {}
for line in subprocess.check_output(['twoBitInfo', sys.argv[1], 'stdout']).decode('utf-8').split('\n')[:-1]:
    id, len = line.split('\t')
    dict_chroms_lengths[id] = int(len)
pickle.dump(dict_chroms_lengths, open(sys.argv[2], 'wb'))
