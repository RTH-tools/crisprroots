#!/usr/bin/env python3
import argparse as ap
import os
import warnings

parser = ap.ArgumentParser()
parser.add_argument('-c', '--CRISPRroots', help='Global path to CRISPRroots folder', type=str, required=True)
parser.add_argument('-d', '--data_directory', help='Global path to the directory containing the data, if different from'
                                                   '<current_directory>/QPRT_DEL268T_chr16_10M-40M', type=str,
                    required=False)
parser.add_argument('-r', '--resources_directory', help='Global path to the directory containing the genome resources, '
                                                        'if different from <current_directory>/resources', type=str,
                    required=False)
args = parser.parse_args()

CRISPRroots_test_directory = os.path.dirname(os.path.abspath(__file__))
config_content = []
path_to_config= 'QPRT_DEL268T_chr16_10M-40M/config.yaml' if not args.data_directory else os.path.join(args.data_directory, 'config.yaml')
with open(path_to_config, 'r+') as conf:
    for line in conf:
        line = line.replace('<set_global_path_to_data_folder>', (
            os.path.join(CRISPRroots_test_directory,'QPRT_DEL268T_chr16_10M-40M') if not args.data_directory else args.data_directory))
        line = line.replace('<set_global_path_to_resources_folder>', (
            os.path.join(CRISPRroots_test_directory,'resources') if not args.resources_directory else args.resources_directory))
        line = line.replace('<set_global_path_to_CRISPRroots>', args.CRISPRroots.rstrip('/'))
        if '//' in line:
            warnings.warn('One or more of the paths generated contain double forward slashes (//). These are converted to a single forward slash (/).')
            line=line.replace('//','/')
        config_content.append(line)
    conf.seek(0)
    conf.writelines(config_content)
print("The configuration file has been updated.")
