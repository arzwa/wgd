#!/usr/bin/python3.5
"""
--------------------------------------------------------------------------------

Copyright (C) 2018 Arthur Zwaenepoel

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Contact: arzwa@psb.vib-ugent.be

--------------------------------------------------------------------------------
"""

from .utils import process_gene_families, log_subprocess, uniq_id
import os
import subprocess
import logging


def all_v_all_blast(query, db, output_directory='./', output_file='blast.tsv',
                    eval_cutoff=1e-10, n_threads=4):
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
    command = ['blastp', '-db', db, '-query', query, '-evalue',
               str(eval_cutoff), '-outfmt', '6', '-num_threads', str(n_threads),
               '-out', os.path.join(output_directory, output_file)]
    logging.info(' '.join(command))
    sp = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    log_subprocess('blastp', sp)
    logging.info("All versus all Blastp done")

    # remove database and other stuff
    subprocess.run(['rm', db + '.phr', db + '.pin', db + '.psq'])

    return os.path.join(output_directory, output_file)


def get_one_v_one_orthologs_rbh(blast_file, output_dir):
    """
    Get one-vs-one orthologs (using RBHs). Implemented for an arbitrary number
    of species. note that every gene ID in the blast file should be prefixed
    with a species ID e.g. ``ath|AT1G01000``.

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
            e = float(line[10])

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
    seen = set()
    for comb, d in one_v_one_orthologs.items():
        logging.info('Writing one vs one orthologs to {}.tsv'.format(comb))
        with open(os.path.join(output_dir, '{}.ovo.tsv'.format(comb)),
                  'w') as o:
            for key, val in d.items():
                if val[0] not in d.keys():
                    logging.warning('Gene {} not found in dictionary strangely '
                                    'enough?'.format(val[0]))
                if key == d[val[0]][0]:  # RBH
                    kv = sorted([key, val[0]])
                    id_ = "\t".join(kv)
                    if id_ not in seen:  # prevent duplicate entries
                        o.write('{0}\n'.format(id_))
                        seen.add(id_)
                else:
                    logging.debug('Best hit for {0} is {1} but for {1} is {2}'
                                  ''.format(key, val, d[val[0]][0]))
        last = os.path.join(output_dir, '{}.ovo.tsv'.format(comb))
    return last


def ava_blast_to_abc(ava_file, col_1=0, col_2=1, col_3=10):
    """
    Convert tab separated all-vs-all (ava) Blast results to an input graph in
    ``abc`` format for ``mcl``.

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


def run_mcl_ava(graph, output_dir='./', inflation=2, output_file='out.mcl',
                preserve=False, return_dict=False):
    """
    Run ``mcl`` on all-vs-all Blast results for a species of interest.
    Note if the parameter ``output_file`` is not given and the parameter
    ``return_dict`` is set to True, only a python dictionary is returned and no
    file is written.

    :param graph: input graph (abc format)
    :param output_dir: output_dir
    :param output_file: output_file (optional)
    :param inflation: inflation factor for ``mcl``
    :param return_dict: boolean, return results as dictionary?
    :param preserve: boolean, preserve tmp/intermediate files?
    :return: results as output file or as gene family dictionary
    """
    tmp_file = uniq_id()
    tmp_file = os.path.join(output_dir, tmp_file)
    output_file = os.path.join(output_dir, output_file)

    # get all-vs-all results in abc format graph for mcl
    with open(tmp_file, 'w') as o:
        o.write("\n".join(["\t".join(x) for x in graph]))

    # run mcl pipeline
    logging.info('Started MCL clustering (mcl)')
    command = ['mcxload', '-abc', tmp_file, '--stream-mirror',
               '--stream-neg-log10',
               '-o', tmp_file + '.mci', '-write-tab', tmp_file + '.tab']
    logging.debug(" ".join(command))
    completed = subprocess.run(command, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    log_subprocess('mcxload', completed)

    command = ['mcl', tmp_file + '.mci', '-I', str(inflation), '-o',
               tmp_file + '.mci.I' + str(inflation * 10)]
    logging.debug(" ".join(command))
    completed = subprocess.run(command, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    log_subprocess('mcl', completed)

    command = ['mcxdump', '-icl', tmp_file + '.mci.I' + str(inflation * 10),
               '-tabr', tmp_file + '.tab',
               '-o', output_file]
    logging.debug(" ".join(command))
    completed = subprocess.run(command, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    log_subprocess('mcxdump', completed)

    # remove temporary files
    if not preserve:
        os.system("rm {}*".format(tmp_file))

    # return output as desired
    if return_dict:
        results = process_gene_families(output_file)
        return results

    return output_file



