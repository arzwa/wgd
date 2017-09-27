#!/usr/bin/python3.5
"""
Arthur Zwaenepoel
"""
import os
import logging
import re
import warnings


def get_gfs_for_species(gene_family_dict, gene_pattern):
    """
    Get non-singleton gene families for a species of interest
    :param gene_family_dict: gene family dictionary
    :param species: species of interest
    :return: dictionairy with gene families
    """
    gene_pattern = re.compile(gene_pattern)
    gf_dict = {}

    for key in gene_family_dict.keys():
        if gene_pattern.search(' '.join(gene_family_dict[key])) is not None:
            gf_dict[key] = gene_family_dict[key]

    return gf_dict


def get_sequences(paralog_dict, sequences):
    """
    Fetch sequences from a fasta file or sequence dict and put them in a two level dictionairy
    {gene_family: {gene: seq, gene: seq, ...}, ...}
    :param paralog_dict:
    :return: two-level dictionairy
    """

    if not type(sequences) == dict:
        sequences = read_fasta(sequences, split_on_pipe=True)

    paralog_sequence_dict = {}

    for family in paralog_dict:
        paralog_sequence_dict[family] = {}
        for gene in paralog_dict[family]:
            if gene not in sequences.keys():
                warnings.warn("Gene {} in gene families but not in protein fasta!".format(gene))
            else:
                paralog_sequence_dict[family][gene] = sequences[gene]

    return paralog_sequence_dict


def process_gene_families(gene_family_file, ignore_prefix=False):
    """
    Processes a raw gene family file as e.g. from OrthoMCL into a generic dictionary structure
    OrthoMCL raw file consists of one gene family per line, including tab separated gene IDs,
    without gene family ID.
    """
    gene_family_dict = {}
    ID = 1

    with open(gene_family_file, 'r') as f:
        for line in f:
            genes = line.strip().split("\t")
            if ignore_prefix:
                if '|' in genes[0]:
                    genes = [gene.split('|')[1] for gene in genes]
            gene_family_dict["GF_{}".format(ID)] = genes
            ID += 1

    return gene_family_dict


def check_dirs(tmp_dir, output_dir, prompt, preserve):
    """
    Check directories needed
    :param tmp_dir: tmp directory
    :param output_dir: output directory
    :param prompt: prompt for overwrites (boolean)?
    :param preserve: preserve MSA files (boolean)?
    :return: nothing
    """
    # Check the tmp directory
    if tmp_dir:
        if os.path.exists(tmp_dir):
            if prompt:
                overwrite = input("tmp directory {} already exists. Overwrite? [y/n]: ".format(tmp_dir))
                if overwrite == 'y':
                    os.system('rm -r {}'.format(os.path.join(tmp_dir)))
                else:
                    print('EXIT')
                    return
            else:
                os.system('rm -r {}'.format(os.path.join(tmp_dir)))
        os.mkdir(tmp_dir)

    # Check the output directory
    if output_dir:
        if os.path.exists(output_dir):
            if prompt:
                overwrite = input("Output directory {} already exists. Overwrite? [y/n]: ".format(output_dir))
                if overwrite == 'y':
                    os.system('rm -r {}'.format(os.path.join(output_dir)))
                else:
                    print('EXIT')
                    return
            else:
                os.system('rm -r {}'.format(os.path.join(tmp_dir)))
        os.mkdir(output_dir)

    # If data should be preserved, make/clean directories
    if preserve:
        if 'msa' not in os.listdir(output_dir):
            os.mkdir(os.path.join(output_dir, 'msa'))

        if 'codeml' not in os.listdir(output_dir):
            os.mkdir(os.path.join(output_dir, 'codeml'))


def read_fasta(fasta_file, prefix=None, split_on_pipe=True, split_on_whitespace=True, raw=False):
    """
    Generic fastafile reader
    Returns a dictionairy {ID: sequence, ID: sequence, ...}.
    """
    sequence_dict = {}

    with open(fasta_file, 'r') as f:
        if raw:
            return f.read()

        genes = f.read().split(">")
        for gene in genes:
            ID = gene.split("\n")[0]
            if ID != '':
                if split_on_pipe:
                    ID = ID.split("|")[0].strip()
                if split_on_whitespace:
                    ID = ID.split()[0]
                if prefix and prefix != '':
                    ID = prefix + '|' + ID
                sequence = "".join(gene.split("\n")[1:])
                sequence = sequence.replace('*', '')
                sequence_dict[ID] = sequence

    if '' in sequence_dict.keys():
        del sequence_dict['']
    return sequence_dict


def translate_cds(sequence_dict):
    """
    Just another CDS to protein translater

    :param sequence_dict: dictionary with gene IDs and CDS sequences
    :return: dictionary with gene IDs and proteins sequences
    """
    aa_dict = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '', 'TAG': '',
        'TGC': 'C', 'TGT': 'C', 'TGA': '', 'TGG': 'W',
    }
    protein_dict = {}

    for key, val in sequence_dict.items():
        aa_seq = ''

        for i in range(0,len(val),3):
            if val[i:i+3] not in aa_dict.keys():
                logging.warning('Invalid codon {0:>3} in sequence {1}'.format(val[i:i+3], key))
            else:
                aa_seq += aa_dict[val[i:i+3]]
        protein_dict[key] = aa_seq

    return protein_dict


def write_fasta(seq_dict, output_file):
    """
    Write a sequence dictionary to a fasta file

    :param seq_dict:
    :param output_file:
    """
    with open(output_file, 'w') as o:
        for key, val in seq_dict.items():
            o.write('>' + key + '\n')
            o.write(val + '\n')


def prefix_fasta(prefix, input_file, output_file):
    d = read_fasta(input_file)
    with open(output_file, 'w') as o:
        for key, val in d.items():
            o.write('>' + prefix + '|' + key + '\n')
            o.write(val + '\n')


def prefix_mcl(prefixes_dict, input_file, output_file):
    with open(output_file, 'w') as o:
        with open(input_file, 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                new_line = []
                for gene in line:
                    prefix = regex_matcher(gene, prefixes_dict)
                    if prefix:
                        new_line.append(prefix + '|' + gene )
                    else:
                        logging.warning("No prefix match for gene {}".format(gene))
                o.write('\t'.join(new_line))
                o.write('\n')


def prefix_multi_fasta(prefixes_dict, input_file, output_file):
    d = read_fasta(input_file)
    with open(output_file, 'w') as o:
        for key, val in d.items():
            prefix = regex_matcher(key, prefixes_dict)
            if prefix:
                o.write('>' + prefix + '|' + key + '\n')
                o.write(val + '\n')
            else:
                logging.warning("No prefix match for gene {}".format(key))


def regex_matcher(gene, r_dict):
    for key in r_dict.keys():
        if r_dict[key].match(gene):
            return key
    return None


def filter_one_vs_one_families(gene_families, s1, s2):
    """
    Filter one-vs-one ortholog containing families for two given species.

    :param gene_families:
    :param s1:
    :param s2:
    :return:
    """
    to_delete = []
    for key, val in gene_families.items():
        count = 0
        for gene in val:
            prefix=gene.split('|')[0]
            if prefix == s1 or prefix == s2:
                count += 1
        if count != 2:
            to_delete.append(key)
    for k in to_delete:
        del gene_families[k]
    return gene_families