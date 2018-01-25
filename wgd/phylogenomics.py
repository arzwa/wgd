#!/usr/bin/python3
"""
Arthur Zwaenepoel - 2017

This is part of the `wgd` package for whole genome duplication analysis.

Not operative yet, not sure if it will be completed.
A good algorithm for phylogenomic WGD inference is GRAMPA (Gregg et al. (2017)).
So it seems there is not much need for other algorithms/implementations here.
"""
from .utils import read_fasta, process_gene_families, get_sequences
from ete3 import Tree
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
    gene_family_dict = process_gene_families(gene_families)
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


def pug(paralogs, trees, species_tree, outgroup_list):
    """
    The PUG, or a PUG-like, algorithm for relative dating of WGDs (McKain et al. 2015).
    """
    count_dict = {}

    for tree_file in trees:
        tree = Tree(tree_file)
        tree = root_to_outgroup(tree, outgroup_list)
        for gene in tree.get_leaves():
            if gene.name in paralogs.keys():
                if paralogs[gene.name] not in [x.name for x in tree.get_leaves()]:
                    logging.info('Paralog not in tree for tree {}!'.format(tree_file))
                    continue

                # we have a tree of interest
                lca = tree.get_common_ancestor([gene.name, paralogs[gene.name]])
                sister = tree.get_common_ancestor([gene.name, paralogs[gene.name]]).get_sisters()
                sister_sp = species_tree.get_common_ancestor([gene.name.split('_')[0],
                                                              paralogs[gene.name].split('_')[0]]).get_sisters()

                if len(sister) > 1:
                    logging.info("Seems like there's a polytomy in tree {}".format(tree_file))

                dup_clade = set(get_species(lca))
                sister_clade = set(get_species(sister[0]).keys())
                sister_sp = set(sister_sp)

                # (1) Check bootstrap support

                # (2) A minimum of two taxa has to be present in the child clade of the duplication node
                if len(list(dup_clade)) < 2:
                    logging.info('To few taxa in child clade for tree {}'.format(tree_file))
                    continue

                # (3) None of the taxa in the child clade of the duplication node can be found
                # in the sister lineage of the duplication node
                if dup_clade & sister_clade != set():
                    logging.info('Child clade taxa found in sister lineage for tree {}'.format(tree_file))
                    continue

                # (4) The sister lineage of the duplication node should have at least one taxon from the
                # sister lineage to the hypothesized WGD node in the species tree
                if sister_clade & sister_sp != set():
                    logging.info("No overlap between sister lineage in species \
                                 tree and sister lineage in gene tree {}".format(tree_file))
                    continue

                lca_name = '-'.join(sorted(list(get_species(lca).keys())))
                if lca_name not in count_dict.keys():
                    count_dict[lca_name] = 0
                count_dict[lca_name] += 1

    return count_dict


def get_species(tree):
    """
    Get species from an OrthoFinder gene tree (prefixed with 'species ID' + '_')
    """
    return {x.name.split('_')[0]: x.name for x in tree.get_leaves()}


def root_to_outgroup(gene_tree, outgroup_list):
    """
    Root a gene tree to the farthest possible outgroup based on species list
    (farthest first in list). Currently will not work with rooting to a non-
    leaf node.
    """
    sp_in_gene_tree = get_species(gene_tree)
    for sp in outgroup_list:
        if sp in sp_in_gene_tree.keys():
            gene_tree.set_outgroup(sp_in_gene_tree[sp])
            return gene_tree


def paralog_dict(paralog_df, prefix=''):
    """
    Get a dictionary with all paralogs referring to each other.
    """
    paralogs = {}
    for row in paralog_df.index:
        p1, p2 = row.split('-')
        paralogs[prefix + p1] = prefix + p2
        paralogs[prefix + p2] = prefix + p1
    return paralogs
