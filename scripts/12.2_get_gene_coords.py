#!/usr/bin/env python3

'''
This script retrieves the promoters and genes coordinates of all genes specified in the file given as --diffexp_genes.
'''

# **********************************************************************************************************************

import pandas as pd
import argparse as ap
import os
import pickle


# **********************************************************************************************************************


def get_chroms_lengths():
    return pickle.load(open(args.dict_genome_info, 'rb'))


def get_feature(info, field):
    info = info.replace(' ', '').replace('"', '').split(';')
    for i in info:
        if field in i:
            return i.replace(field, '')


def retrieve_coordinates():
    dict_chroms_lengths = get_chroms_lengths()
    df = pd.read_table(args.annotations, delimiter='\t', header=None, comment='#')
    df_genes = df[df[2] == 'gene'].drop_duplicates()[[0, 1, 3, 4, 6, 8]]
    df_genes.columns = ['Chr', 'Source', 'Start', 'End', 'Strand', 'Info']
    df_genes['gene_id'] = df_genes['Info'].apply(lambda x: get_feature(info=x, field='gene_id'))
    df_genes.set_index('gene_id', inplace=True)
    df_genes.drop('Info', axis=1, inplace=True)
    lst_df_diffexp_genes = []
    for path_diffexp_genes in args.diffexp_genes:
        if os.stat(path_diffexp_genes).st_size == 0:
            df_diffexp_gene = pd.DataFrame()
        else:
            df_diffexp_gene = pd.read_table(path_diffexp_genes, sep='\t', skiprows=2, index_col='Gene')
        df_diffexp_gene['gene_name'] = [x.split('|')[1].replace('NA;', '').split(';')[0] for x in
                                        list(df_diffexp_gene.index)]
        df_diffexp_gene.index = [x.split('|')[0] for x in list(df_diffexp_gene.index)]
        df_diffexp_gene = df_diffexp_gene.join(df_genes)
        df_diffexp_gene['Expression'] = os.path.basename(path_diffexp_genes).split('_')[1][:-4]
        lst_df_diffexp_genes.append(df_diffexp_gene)
    with open(args.output, 'w') as output_handle:
        for df_genes in lst_df_diffexp_genes:
            for id, line in df_genes.iterrows():
                if line['Chr'] in dict_chroms_lengths.keys():  # Only search primary reference sequence
                    if line['Strand'] == '+':
                        # PROMOTER CASE
                        if line['Start'] - 1 - 1 - args.length <= 0:  # Ensamble has INCLUSIVE coordinates STARTING AT 1
                            start = 0
                        else:
                            start = line['Start'] - 1 - 1 - args.length
                        output_handle.write(
                            '%s\t%i\t%i\tSource=%s;GeneID=%s;GeneName=%s;Type=Promoter;DE=%s;L2FC=%f;OR=%f;ER=%f\t.\t%s\n' % (
                                line['Chr'],  # Chromosome
                                start,  # Start position (0-based, inclusive)
                                line['Start'] - 1 - 1,
                                # End position (0-based, last position of promoter is = last position of gene -1)
                                line['Source'],
                                id,  # Gene ID
                                line['gene_name'],  # Gene ID
                                line['Expression'],  # Expression level
                                line['L2FC'],  # Expression fold change between Edited and Original
                                line['Basemean exp. Original'],  # Expression level
                                line['Basemean exp. Edited'],  # Expression level
                                line['Strand']  # strand
                            ))
                    else:
                        # PROMOTER CASE
                        if line['End'] + args.length >= dict_chroms_lengths[line['Chr']]:
                            end = dict_chroms_lengths[line['Chr']]  # Ensamble has INCLUSIVE coordinates STARTING AT 1
                        else:
                            end = line['End'] + args.length
                        output_handle.write(
                            '%s\t%i\t%i\tSource=%s;GeneID=%s;GeneName=%s;Type=Promoter;DE=%s;L2FC=%f;OR=%f;ER=%f\t.\t%s\n' % (
                                line['Chr'],  # Chromosome
                                line['End'],  # Start position of promoter (0-based, inclusive)
                                end,  # End position (0-based, last position of gene)
                                line['Source'],
                                id,  # Gene ID
                                line['gene_name'],  # Gene ID
                                line['Expression'],  # Expression level
                                line['L2FC'],  # Expression fold change between Edited and Original
                                line['Basemean exp. Original'],  # Gene exp
                                line['Basemean exp. Edited'],  # Gene exp
                                line['Strand']  # strand
                            )
                        )
                    # GENE CASE
                    output_handle.write(
                        '%s\t%i\t%i\tSource=%s;GeneID=%s;GeneName=%s;Type=Gene;DE=%s;L2FC=%f;OR=%f;ER=%f\t.\t%s\n' % (
                            line['Chr'],  # Chromosome
                            line['Start'] - 1,  # Start position (0-based, inclusive)
                            line['End'] - 1,  # End position (0-based, last position of gene, inclusive)
                            line['Source'],
                            id,  # Gene ID
                            line['gene_name'],  # Gene ID,
                            line['Expression'],  # Expression level
                            line['L2FC'],  # Expression fold change between Edited and Original
                            line['Basemean exp. Original'],  # Gene exp
                            line['Basemean exp. Edited'],  # Gene exp
                            line['Strand']  # strand
                        )
                    )


# **********************************************************************************************************************

parser = ap.ArgumentParser()
parser.add_argument('-a', '--annotations', help='Annotations file in gtf format', type=str, required=True)
parser.add_argument('-de', '--diffexp_genes',
                    help='Table containing strongly differentially expressed genes to be analyzed',
                    type=str, required=True, nargs='+')
parser.add_argument('-o', '--output', help='Output file name', type=str, required=True)
parser.add_argument('-l', '--length', help='Length of the promoter region', type=int,
                    default=1000)
parser.add_argument('-d', '--dict_genome_info', help='Path to dictionary of chromosome length', type=str)

args = parser.parse_args()
retrieve_coordinates()
