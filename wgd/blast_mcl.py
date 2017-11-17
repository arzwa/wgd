#!/usr/bin/python3.5
"""
Arthur Zwaenepoel

Markov clustering (MCL) related functions
"""

from .utils import process_gene_families, translate_cds, read_fasta
import os
import subprocess
import re
import logging
import pandas as pd
import numpy as np


def _ava_blast_to_abc(ava_file, regex='.+', prefix=None, col_1=0, col_2=1, col_3=2):
    """
    Convert tab separated all-vs-all (ava) Blast results to an input graph for ``mcl``.

    :param prefix: prefix for genes of interest (optional), e.g. 'arat' for 'arat|AT1G001290'.
    :param regex: regular expression for parsing genes of interest (optional), e.g. 'AT.+' for 'AT1G001290'.
    :param ava_file: all-vs-all Blast results file (tab separated)
    :param pipe: boolean, piped prefix or not? (see species parameter)
    :param col_1: column in file for gene 1 (starts from 0)
    :param col_2: column in file for gene 2
    :param col_3: column in file for e-value (weight in graph)
    :return: graph in ``abc`` format for ``mcl``
    """

    if not prefix:
        if not regex:
            regex = '.+'
        p = re.compile(regex)

    graph = []
    with open(ava_file, 'r') as f:
        for line in f:
            line = line.split("\t")
            gene_1, gene_2, e = line[col_1], line[col_2], line[col_3]

            # when prefixed with species id (e.g. arat|AT1G00120)
            if prefix:
                gene_1 = gene_1.split('|')
                gene_2 = gene_2.split('|')

                if gene_1[0] == gene_2[0] == prefix:
                    graph.append([gene_1[1], gene_2[1], str(e)])

            # when not prefixed but regex patter defined
            else:
                if p.match(gene_1) and p.match(gene_2):
                    graph.append([gene_1, gene_2, str(e)])

    return graph


def run_mcl_ava(input_file, regex='.+', prefix=None, tmp_dir='./', output_file='mcl.out', inflation=2,
                return_dict=False, preserve=True, **kwargs):
    """
    Run ``mcl`` on all-vs-all Blast results for a species of interest.
    Note if the parameter ``output_file`` is not given and the parameter ``return_dict`` is set to True,
    only a python dictionary is returned and no file is written.

    :param prefix: prefix for genes of interest (optional), e.g. 'arat' for 'arat|AT1G001290'.
    :param regex: regular expression for parsing genes of interest (optional), e.g. 'AT.+' for 'AT1G001290'.
    :param input_file: all-vs-all Blast input file
    :param tmp_dir: directory to store temporary/intermediate files
    :param output_file: output_file (optional)
    :param inflation: inflation factor for ``mcl``
    :param return_dict: boolean, return results as dictionary?
    :param preserve: boolean, preserve tmp/intermediate files?
    :param kwargs: other arguments for :py:func:`_ava_blast_to_abc`
    :return: results as output file or as gene family dictionary
    """
    # get all-vs-all results in abc format graph for mcl
    logging.info('Making input graph for mcl')
    graph = _ava_blast_to_abc(input_file, regex=regex, prefix=prefix, **kwargs)

    tmp_file = os.path.basename(input_file) + '.abc'
    tmp_file = os.path.join(tmp_dir, tmp_file)

    with open(tmp_file, 'w') as o:
        o.write("\n".join(["\t".join(x) for x in graph]))

    # configuration

    # run mcl pipeline
    logging.info('Started MCL clustering (mcl)')
    command = ['mcxload', '-abc', tmp_file, '--stream-mirror', '--stream-neg-log10',
               '-o', tmp_file + '.mci', '-write-tab', tmp_file + '.tab']
    logging.debug(" ".join(command))
    completed = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logging.debug(completed.stderr.decode('utf-8'))

    command = ['mcl', tmp_file + '.mci', '-I', str(inflation), '-o', tmp_file + '.mci.I' + str(inflation*10)]
    logging.debug(" ".join(command))
    completed = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logging.debug(completed.stderr.decode('utf-8'))

    command = ['mcxdump', '-icl', tmp_file + '.mci.I' + str(inflation*10), '-tabr', tmp_file + '.tab',
               '-o', output_file]
    logging.debug(" ".join(command))
    completed = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logging.debug(completed.stderr.decode('utf-8'))

    # remove temporary files
    if not preserve:
        os.system("rm {}*".format(tmp_file))

    # return output as desired
    if return_dict:
        results = process_gene_families(output_file)

        # if dict asked and no output file given by user, remove output file
        if output_file == 'mcl.out':
            os.remove(output_file)

        return results

    return output_file


def all_v_all_blast(query, db, output_directory, eval_cutoff=1e-10):
    """
    Perform all-versus-all Blastp.

    :param seq_file: fasta file with protein sequences
    :param output_directory: output directory
    :param eval_cutoff: e-value cut off for Blastp results
    :return: all-`versus`-all Blastp results
    """

    logging.info("Making Blastdb")
    subprocess.run(['makeblastdb', '-in', db, '-dbtype', 'prot'])

    logging.info("Running Blastp")
    outfile = os.path.join(output_directory, os.path.basename(db) + '.out.blast')
    subprocess.run(['blastp', '-db', db, '-query', query, '-evalue', str(eval_cutoff), '-outfmt', '6',
                    '-out', outfile])
    logging.info("All versus all Blastp done")

    logging.info("Reformatting output")
    os.system(" ".join(['cut', '-f1,2,11', outfile, '> ava.tsv'])) # this did'nt seem to work with suprocess.run ?
    subprocess.run(['mv', 'ava.tsv', outfile])
    subprocess.run(['rm', db + '.phr', db + '.pin', db + '.psq'])

    return outfile


# REWRITE: MAINTAIN BOTH IN ORDER TO PREVENT BROKEN STUFF --------------------------------------------------------------

def ava_blast_to_abc_2(ava_file, col_1=0, col_2=1, col_3=2):
    """
    Convert tab separated all-vs-all (ava) Blast results to an input graph for ``mcl``.

    :param ava_file: all-vs-all Blast results file (tab separated)
    :param col_1: column in file for gene 1 (starts from 0)
    :param col_2: column in file for gene 2
    :param col_3: column in file for e-value (weight in graph)
    :return: graph in ``abc`` format for ``mcl``
    """
    graph = []
    with open(ava_file, 'r') as f:
        for line in f:
            line = line.split("\t")
            graph.append([line[col_1], line[col_2], line[col_3]])
    return graph


def run_mcl_ava_2(graph, output_dir='./', inflation=2, output_file='out.mcl', preserve=False, return_dict=False):
    """
    Run ``mcl`` on all-vs-all Blast results for a species of interest.
    Note if the parameter ``output_file`` is not given and the parameter ``return_dict`` is set to True,
    only a python dictionary is returned and no file is written.

    :param graph: input graph (abc format)
    :param output_dir: output_dir
    :param output_file: output_file (optional)
    :param inflation: inflation factor for ``mcl``
    :param return_dict: boolean, return results as dictionary?
    :param preserve: boolean, preserve tmp/intermediate files?
    :return: results as output file or as gene family dictionary
    """
    # get all-vs-all results in abc format graph for mcl

    tmp_file = 'mcl.abc'
    tmp_file = os.path.join(output_dir, tmp_file)
    output_file = os.path.join(output_dir, output_file)

    with open(tmp_file, 'w') as o:
        o.write("\n".join(["\t".join(x) for x in graph]))

    # configuration

    # run mcl pipeline
    logging.info('Started MCL clustering (mcl)')
    command = ['mcxload', '-abc', tmp_file, '--stream-mirror', '--stream-neg-log10',
               '-o', tmp_file + '.mci', '-write-tab', tmp_file + '.tab']
    logging.debug(" ".join(command))
    completed = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logging.debug(completed.stderr.decode('utf-8'))

    command = ['mcl', tmp_file + '.mci', '-I', str(inflation), '-o', tmp_file + '.mci.I' + str(inflation*10)]
    logging.debug(" ".join(command))
    completed = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logging.debug(completed.stderr.decode('utf-8'))

    command = ['mcxdump', '-icl', tmp_file + '.mci.I' + str(inflation*10), '-tabr', tmp_file + '.tab',
               '-o', output_file]
    logging.debug(" ".join(command))
    completed = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logging.debug(completed.stderr.decode('utf-8'))

    # remove temporary files
    if not preserve:
        os.system("rm {}*".format(tmp_file))

    # return output as desired
    if return_dict:
        results = process_gene_families(output_file)
        return results

    return output_file
