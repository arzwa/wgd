#!/usr/bin/python3.5
"""
Arthur Zwaenepoel
"""
from .ks_distribution import analyse_family
from .utils import read_fasta, process_gene_families_ortho_mcl, get_sequences
from .utils import get_gfs_for_species, check_dirs, translate_cds
import pandas as pd
import numpy as np
import logging
import asyncio
import re
import os
from scipy.stats import t


def _mean_for_genes(data_frame, gene_pattern):
    """
    Get the mean omega value for all genes of a species
    :param data_frame:
    :param gene_pattern:
    :return:
    """
    if data_frame.empty:
        return None

    genes = set([gene for gene in list(data_frame['Paralog1']) + list(data_frame['Paralog2'])
                 if gene_pattern.match(gene)])

    results = pd.DataFrame(columns=['Omega', 'n_pairs', 'std', 'CI95_lower', 'CI95_upper'],
                           index=genes)

    # calculate means etc. for genes
    for gene in genes:
        data = list(data_frame[data_frame['Paralog1'] == gene]['Omega']) + \
               list(data_frame[data_frame['Paralog2'] == gene]['Omega'])
        mean = sum(data)/len(data)
        sigma = np.sqrt(sum([(mean-x)**2 for x in data])/len(data))
        ci = t.interval(0.95, len(data) - 1, loc=mean, scale=sigma)
        results.set_value(gene, 'Omega', mean)
        results.set_value(gene, 'n_pairs', len(data))
        results.set_value(gene, 'std', sigma)
        results.set_value(gene, 'CI95_lower', ci[0])
        results.set_value(gene, 'CI95_upper', ci[1])

    return results


class PositiveSelection():
    """
    Class for performing positive selection screen.
    Highly similar to :py:class:`KsDistribution`
    Required parameters:

    :param species: species id
    :param gene_pattern: gene pattern (regex for parsing out species specfic genes from gene families)
    :param gene_families: gene families file
    :param nucleotide: fasta file containing all nucleotide sequences (CDS) for all species.
    :param protein: fasta file containing all protein sequences for all species.
    """

    def __init__(self, gene_pattern, gene_families, nucleotide, protein, species='pos',
                 muscle_path='muscle', codeml_path='codeml', tmp_dir='./', output_dir='./'):
        """
        init for positive selection class.

        :param gene_pattern: gene pattern (regex for parsing out species specfic genes from gene families)
        :param gene_families: gene families file
        :param nucleotide: fasta file containing all nucleotide sequences (CDS) for all species
        :param protein: fasta file containing all protein sequences for all species
        :param species: species id
        :param muscle_path: path to muscle executable
        :param codeml_path: path to codeml executable
        :param tmp_dir: path to tmp directory
        :param output_dir: path to output directory
        """
        self.species = species
        self.gene_pattern = re.compile(gene_pattern)
        self.gene_families = gene_families
        self.nucleotide_file = nucleotide
        self.protein_file = protein
        self.output_dir = output_dir
        self.muscle_path = muscle_path
        self.codeml_path = codeml_path
        self.tmp_dir = tmp_dir
        self.gf_dict = None
        self.paralogs = None
        self.nucleotide = None
        self.protein = None
        self.mixture_models = None
        self.kde = None
        self.counts = 0
        self.results = None

    def positive_selection(self, check=True, preserve=False, prompt=True, times=1):
        """
        Calculates for each non-singleton gene of a species of interest the mean Kn/Ks (omega) value.
        This is performed by performing codeml analysis.

        :param check: boolean, check directories if present/empty?
        :param preserve: boolean, preserve intermediate files?
        :param prompt: boolean, prompt for possible overwriting?
        :return: data frame
        """
        logging.debug('Pre-processing sequence and gene family data')
        self.gf_dict = process_gene_families_ortho_mcl(self.gene_families)
        species_gfs = get_gfs_for_species(self.gf_dict, self.gene_pattern)
        self.nucleotide = get_sequences(species_gfs, self.nucleotide_file)

        if not self.protein_file:
            self.protein = translate_cds(self.nucleotide)
        else:
            self.protein = get_sequences(species_gfs, self.protein_file)


        if check:
            logging.debug("Checking directories")
            check_dirs(self.tmp_dir, self.output_dir, prompt, preserve)

        logging.info('Started analysis in parallel')
        loop = asyncio.get_event_loop()
        tasks = [loop.run_in_executor(None, analyse_family, family, self.protein[family],
                                      self.nucleotide, self.tmp_dir, self.muscle_path,
                                      self.codeml_path, preserve, times) for family in self.protein.keys()]

        loop.run_until_complete(asyncio.gather(*tasks))
        logging.info('Analysis done')

        if preserve:
            logging.info('Moving files (preserve=True)')
            os.system(" ".join(['mv', os.path.join(self.tmp_dir, '*.msa*'), os.path.join(self.output_dir, 'msa')]))
            os.system(" ".join(['mv', os.path.join(self.tmp_dir, '*.codeml'), os.path.join(self.output_dir, 'codeml')]))

        else:
            logging.info('Removing files (preserve=False)')
            os.system(" ".join(['rm', os.path.join(self.tmp_dir, '*.msa*')]))

        logging.info('Making results data frame')
        results_frame = pd.DataFrame(columns=['Omega', 'n_pairs', 'std', 'CI95_lower', 'CI95_upper'])

        counts = 0
        for f in os.listdir(self.tmp_dir):
            if f[-3:] == '.Ks':
                counts += 1
                df = pd.read_csv(os.path.join(self.tmp_dir, f), index_col=0)
                results_frame = pd.concat([results_frame, _mean_for_genes(df, self.gene_pattern)])

        return results_frame
