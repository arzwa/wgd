#!/usr/bin/python3.5
"""
Arthur Zwaenepoel - 2017

This is part of the `wgd` package for whole genome duplication analysis.
"""
from .utils import read_fasta, process_gene_families_ortho_mcl, get_sequences
import os
import logging


def generate_fasta_files(gene_families, fasta_files, output_dir='gfs_fasta_files'):
    """
    Generates multi-fasta files as input for an MSA program based on gene families and species fasta files

    :param gene_families: typical OrthoMCL output file, one family per line, with tab-separated gene IDs
    :param fasta_files: comma-separated list of files
    :param output_dir: output directory
    """
    fasta_files_list = fasta_files.strip().split(",")
    gene_family_dict = process_gene_families_ortho_mcl(gene_families)
    sequence_dict = {}

    for fasta_file in fasta_files_list:
        print(fasta_file)
        sequence_dict.update(read_fasta(fasta_file))

    sequences= get_sequences(gene_family_dict, fasta_files)

    for family in sequences.keys():
        logging.debug('Writing {}.fasta'.format(family))
        target_path_fasta = os.path.join(output_dir, '{}.fasta'.format(family))
        with open(target_path_fasta, 'w') as o:
            for gene in sequences[family].keys():
                logging.debug('Writing {0} to {1}.fasta'.format(gene, family))
                o.write('>{}\n'.format(gene))
                o.write(sequences[family][gene])
                o.write('\n')

    return



def ete_build_tree():
    pass
