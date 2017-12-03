#!/usr/bin/python3.5
"""
Arthur Zwaenepoel - 2017

Markov clustering (MCL) related functions
"""

from .utils import process_gene_families
import os
import subprocess
import uuid
import logging


def all_v_all_blast(query, db, output_directory='./', output_file='blast.tsv', eval_cutoff=1e-10, n_threads=4):
    """
    Perform all-versus-all Blastp.
    Runs a blast of ``query`` vs. ``db``.

    :param query: query sequences fasta file
    :param db: database sequences fasta file
    :param output_directory: output directory
    :param output_file: output file name
    :param eval_cutoff: e-value cut-off
    :param: n_threads: number of threads to use for blastp
    :return: blast file path
    """
    logging.info("Making Blastdb")
    subprocess.run(['makeblastdb', '-in', db, '-dbtype', 'prot'])

    logging.info("Running Blastp")
    subprocess.run(['blastp', '-db', db, '-query', query, '-evalue', str(eval_cutoff), '-outfmt', '6',
                    '-num_threads', str(n_threads), '-out', os.path.join(output_directory, output_file)])
    logging.info("All versus all Blastp done")

    logging.info("Reformatting output")
    # this didn't seem to work with suprocess.run ?
    tmp = str(uuid.uuid4())
    os.system(" ".join(['cut', '-f1,2,11', os.path.join(output_directory, output_file), '>', tmp]))
    subprocess.run(['mv', tmp, os.path.join(output_directory, output_file)])
    subprocess.run(['rm', db + '.phr', db + '.pin', db + '.psq'])

    return os.path.join(output_directory, output_file)


def get_one_v_one_orthologs_rbh(blast_file, output_dir):
    """
    Get one-vs-one orthologs (using RBHs).
    Implemented for an arbitrary number of species. note that every gene ID in the blast file
    should be prefixed with a species ID e.g. ``ath|AT1G01000``.

    :param blast_file: all vs. all blastp results, gene IDs should be prefixed
    :param output_dir: output directory
    :return: the last output file that was written
    """
    one_v_one_orthologs = {}
    with open(blast_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            sp_1, gene_1 = line[0].split('|')
            sp_2, gene_2 = line[1].split('|')
            e = float(line[2])

            if sp_1 == sp_2:  # putative paralog
                continue

            comb = '_'.join(sorted([sp_1, sp_2]))
            if comb not in one_v_one_orthologs.keys():
                one_v_one_orthologs[comb] = {}

            if gene_1 not in one_v_one_orthologs[comb].keys():
                one_v_one_orthologs[comb][gene_1] = (gene_2, e)
            elif e < one_v_one_orthologs[comb][gene_1][1]:
                one_v_one_orthologs[comb][gene_1] = (gene_2, e)

            if gene_2 not in one_v_one_orthologs[comb].keys():
                one_v_one_orthologs[comb][gene_2] = (gene_1, e)
            elif e < one_v_one_orthologs[comb][gene_2][1]:
                one_v_one_orthologs[comb][gene_2] = (gene_1, e)

    last = None
    for comb, d in one_v_one_orthologs.items():
        logging.info('Writing one vs one orthologs to {}.tsv'.format(comb))
        with open(os.path.join(output_dir, '{}.ovo.tsv'.format(comb)), 'w') as o:
            for key, val in d.items():
                if val[0] not in d.keys():
                    logging.warning('Gene {} not found in dictionary strangely enough?'.format(val[0]))
                if key == d[val[0]][0]:  # RBH
                    o.write('{0}\t{1}\n'.format(key, val[0]))
                else:
                    logging.debug('Best hit for {0} is {1} but for {1} is {2}'.format(key, val, d[val[0]][0]))
        last = os.path.join(output_dir, '{}.ovo.tsv'.format(comb))
    return last


def ava_blast_to_abc(ava_file, col_1=0, col_2=1, col_3=2):
    """
    Convert tab separated all-vs-all (ava) Blast results to an input graph in ``abc`` format for ``mcl``.

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
    tmp_file = str(uuid.uuid4())
    tmp_file = os.path.join(output_dir, tmp_file)
    output_file = os.path.join(output_dir, output_file)

    # get all-vs-all results in abc format graph for mcl
    with open(tmp_file, 'w') as o:
        o.write("\n".join(["\t".join(x) for x in graph]))

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
