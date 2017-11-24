#!/usr/bin/python3.5
"""
Arthur Zwaenepoel
"""
import os
import logging
import re
import warnings
import random
import json
from progressbar import ProgressBar


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
            gene_family_dict["GF_{:06d}".format(ID)] = genes
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

    with ProgressBar(max_value=len(sequence_dict.keys())+1) as pb:
        j = 0
        for key, val in sequence_dict.items():
            j += 1
            aa_seq = ''
            for i in range(0,len(val),3):
                if val[i:i+3] not in aa_dict.keys():
                    logging.debug('Invalid codon {0:>3} in sequence {1}'.format(val[i:i+3], key))
                else:
                    aa_seq += aa_dict[val[i:i+3]]
            protein_dict[key] = aa_seq
            pb.update(j)

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


def get_number_of_sp(genes):
    """
    Get the number of unique species in a list of gene IDs.
    Will be approximate since it is based on the assumption that
    the leading non-digit part identifies the species
    (which is not always the case e.g. mitochondrial genes etc.)
    """
    p = re.compile('\D*')
    gene_set = set()
    for gene in genes:
        m = p.match(gene)
        if m:
            gene_set.add(m.group())
    return len(list(gene_set))


def check_genes(genes, ids):
    for gene in genes:
        for id in ids:
            if gene.startswith(id):
                return True


def _random_color():
    """
    Generate a random hex color
    """
    def r(): return random.randint(0, 255)
    return '#%02X%02X%02X' % (r(), r(), r())


class Genome:
    """
    Class that represents a structural annotation.
    Collects several nice data structures for a genome and parsers for various
    genomic data file formats (e.g. gff, fasta, ...)
    """

    def __init__(self):
        """
        Genome.genome: dictionary with a full representation
        Genome.gene_lists: dictionary with ordered lists of genes per
        chromosome (as for I-ADHoRe)
        """
        self.parent_file = None
        self.genome = {}
        self.gene_lists = {}
        self.colors = {}

    def parse_plaza_gff(self, gff_file, keyword='mRNA', id_string='Parent'):
        """
        Parse a PLAZA annotation file into a genome dictionary

        :param gff_file: input gff (PLAZA style)
        :param keyword: keyword for elements to parse out
        """
        self.parent_file = gff_file

        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                line = line.strip().split('\t')

                if line[2] == keyword:
                    chromosome = line[0]
                    start = line[3]
                    stop = line[4]
                    orientation = line[6]
                    gene_l = line[8].split(';')
                    gene_dict = {x.split('=')[0]: x.split('=')[1] for x in gene_l if len(x.split('=')) == 2}

                    if chromosome not in self.genome:
                        self.genome[chromosome] = {}
                        self.gene_lists[chromosome] = []
                        self.colors[chromosome] = _random_color()

                    self.genome[chromosome][gene_dict[id_string]] = {
                        'orientation': orientation, 'start': start, 'stop': stop}
                    self.gene_lists[chromosome].append(
                        (gene_dict[id_string], orientation, start, stop))
        return

    def karyotype_json(self, out_file='genome.json'):
        """
        Generate karyotype data file in json format (as per Circos.js/d3.js)
        """
        karyotype = []
        for chrom in self.gene_lists.keys():
            # approximate chromosome length
            coordinates = [int(x[2]) for x in self.gene_lists[chrom]]
            coordinates += [int(x[3]) for x in self.gene_lists[chrom]]
            length = max(coordinates) - min(coordinates)
            karyotype.append({'id': chrom, 'label': chrom,
                              'color': self.colors[chrom], 'len': length})

        if out_file:
            with open(out_file, 'w') as f:
                json.dump(karyotype, f)

        else:
            return json.dumps(karyotype)
