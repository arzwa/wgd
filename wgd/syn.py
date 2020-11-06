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

The co-linearity analysis functions currently only support intragenomic
analyses. It can be a bit finicky, but that's mainly due to I-ADHoRe and people
messing up the GFF format, not me! (definitely also me)
"""
# TODO: multiple genomes

import os
import subprocess
import logging
import pandas as pd
from operator import itemgetter


def gff_parser(gff_file, feature='mRNA', gene_attribute='Parent'):
    # This is a GFF parser with as sole purpose to write out gene list files
    # for I-ADHoRe 3.0. It should work with GFF3, GFF2 and GTF files and follows
    # standard nomenclature (wikipedia). It assumes that the ninth column (the
    # attributes section) is a semicolon separated list of attributes in
    # 'attribute_key=value' format.
    all_features = set()
    genome = {}

    with open(gff_file, 'r') as f:
        for i, line in enumerate(f):

            # ignore comments
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')

            if len(line) != 9:
                raise IndexError

            if line[2] == feature:
                sequence = line[0]
                start = line[3]
                end = line[4]
                strand = line[6]
                if strand != "+" and strand != "-":
                    logging.warning("No strand information for line {}".format(
                        i+1))
                    continue

                # get the feature attributes
                attributes = line[8].split(';')
                attributes = {
                    x.split('=')[0]: x.split('=')[1] for x in attributes
                    if len(x.split('=')) == 2
                }

                if gene_attribute not in attributes:
                    logging.warning('Attribute {0} not found in GFF line {1}'
                                    ''.format(gene_attribute, i))
                    continue

                # store information (why?, if the sole purpose is to write out
                # gene lists we might as well write them here? However if the
                # gff file would not be sorted (which is seldomly the case), it
                # is better to store and sort first
                if sequence not in genome:
                    genome[sequence] = []
                genome[sequence].append((
                    attributes[gene_attribute],
                    int(start), int(end),
                    strand
                ))

                all_features.add(attributes[gene_attribute])

    # sort the features on every sequence
    for v in genome.values():
        v.sort(key=itemgetter(1))

    return genome, all_features


# WRITE FILES AND CONFIG -------------------------------------------------------
def write_gene_lists(genome, output_dir='gene_lists'):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    for k, v in genome.items():
        with open(os.path.join(output_dir, k + '.lst'), 'w') as o:
            for feature in v:
                o.write(feature[0] + feature[-1] + '\n')


def _write_gene_lists(genome, output_dir='gene_lists'):
    """
    Write out the gene lists

    :param genome: Genome object 
    :param output_dir: output directory
    """
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # meanwhile, store all genes in a set to later add as singletons
    all_genes = set()

    for chromosome in genome.genome.keys():
        with open(os.path.join(output_dir, chromosome + '.lst'), 'w') as o:
            for gene in genome.gene_lists[chromosome]:
                o.write(gene[0] + gene[1] + '\n')
                all_genes.add(gene[0])

    return all_genes


def write_families_file(families, all_genes, output_file='families.tsv'):
    """
    Write out families file

    :param families: gene families
    :param all_genes: set object with all genes for the species of interest
    :param output_file: output file name
    :return: nada
    """
    counter = 1
    genes_seen = set()

    with open(families, 'r') as f:
        with open(output_file, 'w') as o:
            for line in f:
                line = line.strip().split('\t')
                for gene in line:
                    o.write(gene + '\t' + str(counter) + '\n')
                    genes_seen.add(gene)
                counter += 1

            # add genes not seen to the families file as singletons
            # I-ADHoRe throws an error when there are genes in the gene
            # lists that are not in the blast table/families file
            rest = all_genes - genes_seen
            for gene in list(rest):
                o.write(gene + '\t' + str(counter) + '\n')
                counter += 1
    return


def write_config_adhore(
        gene_lists, families, config_file_name='i-adhore.conf',
        genome='genome', output_path='i-adhore_out', gap_size=30,
        cluster_gap=35, q_value=0.75, prob_cutoff=0.01, anchor_points=3,
        alignment_method='gg2', level_2_only='false', table_type='family',
        multiple_hypothesis_correction='FDR', visualize_ghm='false',
        visualize_alignment='true'):
    """
    Write out the config file for I-ADHoRe. See I-ADHoRe manual for information
    on parameter settings.

    :param gene_lists: directory with gene lists per chromosome
    :param families: file with gene to family mapping
    :param config_file_name: name for the config file
    :param genome: genome name
    :param output_path: output path name
    :param gap_size: see I-ADHoRe 3.0 documentation
    :param cluster_gap: see I-ADHoRe 3.0 documentation
    :param q_value: see I-ADHoRe 3.0 documentation
    :param prob_cutoff: see I-ADHoRe 3.0 documentation
    :param anchor_points: see I-ADHoRe 3.0 documentation
    :param alignment_method: see I-ADHoRe 3.0 documentation
    :param level_2_only: see I-ADHoRe 3.0 documentation
    :param table_type: see I-ADHoRe 3.0 documentation
    :param multiple_hypothesis_correction: see I-ADHoRe 3.0 documentation
    :param visualize_ghm: see I-ADHoRe 3.0 documentation
    :param visualize_alignment: see I-ADHoRe 3.0 documentation
    :return: configuration file see I-ADHoRe 3.0 documentation
    """
    with open(config_file_name, 'w') as o:
        o.write('genome= {}\n'.format(genome))

        counter = 1
        for l in sorted(os.listdir(gene_lists)):
            o.write("{0} {1}/{2}\n".format(l[:-4], gene_lists, l))
            counter += 1

        o.write('blast_table= {}\n'.format(families))
        o.write('output_path= {}\n'.format(output_path))
        o.write('gap_size= {}\n'.format(gap_size))
        o.write('q_value= {}\n'.format(q_value))
        o.write('cluster_gap= {}\n'.format(cluster_gap))
        o.write('prob_cutoff= {}\n'.format(prob_cutoff))
        o.write('anchor_points= {}\n'.format(anchor_points))
        o.write('alignment_method= {}\n'.format(alignment_method))
        o.write('level_2_only= {}\n'.format(level_2_only))
        o.write('table_type= {}\n'.format(table_type))
        o.write('multiple_hypothesis_correction= {}\n'.format(
                multiple_hypothesis_correction))
        o.write('visualizeGHM= {}\n'.format(visualize_ghm))
        # o.write('visGPairs= {}\n'.format(output_path))
        o.write('visualizeAlignment= {}\n'.format(visualize_alignment))

    return


def get_anchor_pairs(anchors, ks_distribution=None,
                     out_file='anchors_ks.csv'):
    """
    Get anchor pairs and their corresponding Ks values (if provided)

    :param anchors: anchorpoints.txt output from I-ADHoRe 3.0
    :param ks_distribution: Ks distribution dataf rame
    :param out_file: output file name
    :return: pandas dataframe(s): anchors and data frame
    """
    anchors = anchors[['gene_x', 'gene_y']]
    ids = anchors.apply(lambda x: '__'.join(sorted(x)), axis=1)
    
    if type(ks_distribution) == pd.DataFrame:
        ks_anchors = ks_distribution.loc[
            ks_distribution.index.intersection(ids)]

        if out_file:
            ks_anchors.to_csv(out_file, sep='\t')

        return ks_distribution, ks_anchors

    else:
        if out_file:
            anchors.to_csv(out_file, sep='\t')

        return anchors


def segments_to_chords_table(segments_file, genome, output_file='chords.tsv'):
    """
    Create chords table for visualization in a chord diagram. Uses the
    segments.txt output of I-ADHoRe. Chords are defined by a source chromosome
    and a target chromosome and begin and end coordinates for each chromosome
    respectively.

    TODO: the length of each syntenic block should be included in the table as
    well with length defined as number of genes, not physical length.

    :param segments_file: pat to the I-ADHoRe segments.txt output file
    :param genome: a :func:`gff_parser.Genome object`
    :param output_file: output file name
    """
    segments = pd.read_csv(segments_file, sep='\t', index_col=0)
    segments = segments.groupby('multiplicon')
    chromosomes = segments['list'].apply(list)
    first = segments['first'].apply(list)
    last = segments['last'].apply(list)

    chords = []
    for i in range(1, len(chromosomes) + 1):
        for j in range(len(chromosomes[i])):
            for k in range(j + 1, len(chromosomes[i])):
                source = str(chromosomes[i][j])
                target = str(chromosomes[i][k])

                if source not in genome.genome.keys():
                    logging.warning(
                        'Chromosome ID `{}` not found'.format(source))
                elif target not in genome.genome.keys():
                    logging.warning(
                        'Chromosome ID `{}` not found'.format(target))
                else:
                    d = {'source_id': source,
                         'source_1': genome.genome[source][first[i][j]][
                             'start'],
                         'source_2': genome.genome[source][last[i][j]]['stop'],
                         'target_id': target,
                         'target_1': genome.genome[target][first[i][k]][
                             'start'],
                         'target_2': genome.genome[target][last[i][k]]['stop'],
                         'label': '{0}-{1};{2}-{3}'.format(first[i][j],
                                                           last[i][j],
                                                           first[i][k],
                                                           last[i][k]),
                         'color': genome.colors[source],
                         'source_length': abs(
                             int(genome.genome[source][last[i][j]]['stop']) -
                             int(genome.genome[source][first[i][j]]['start']))}
                    chords.append(d)

    df = pd.DataFrame.from_records(chords)
    df.to_csv(output_file, sep='\t')


# RUNNING EXTERNAL SOFTWARE ----------------------------------------------------

def run_adhore(config_file):
    """
    Run I-ADHoRe for a given config file

    :param config_file: path to I-ADHoRe configuration file
    """
    completed = subprocess.run(
            ['i-adhore', config_file], stderr=subprocess.PIPE,
            stdout=subprocess.PIPE)
    logging.warning(completed.stderr.decode('utf-8'))
    logging.info(completed.stdout.decode('utf-8'))
    return
