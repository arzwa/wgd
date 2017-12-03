"""
Arthur Zwaenepoel - 2017
"""
from .utils import read_fasta, write_fasta
from ete3 import Tree
import subprocess
import logging
import numpy as np


def write_sequential_phyml(sequence_dict, output_file):
    """
    Write a multiple sequence alignment in sequential format (e.g. for PhyML)

    :param sequence_dict: sequence dictionary
    :param output_file: filename
    """
    first = True
    with open(output_file, 'w') as o:
        for id, seq in sequence_dict.items():
            if first:
                o.write('{0} {1}\n'.format(len(list(sequence_dict.keys())), len(seq)))
                first = False
            o.write('{0}\t\t{1}\n'.format(id, seq))
    return 0


def run_phyml(msa, phyml_path='phyml'):
    """
    Run PhyML on a protein multiple sequence alignment

    :param msa: file path to protein multiple sequence alignment in multifasta format
    :param phyml_path: path to phyml executable
    :return: path to the tree file
    """
    msa_phyml = msa + '.phyml'
    tree_path = msa + '.nw'
    write_sequential_phyml(read_fasta(msa), msa_phyml)
    command = [phyml_path, '-i', msa_phyml, '-q', '-d', 'aa']

    logging.debug('Running PhyML: {}'.format(' '.join(command)))
    subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run(['rm', msa_phyml, '{}_phyml_stats.txt'.format(msa_phyml)])
    subprocess.run(['mv', '{}_phyml_tree.txt'.format(msa_phyml), tree_path])

    return tree_path


def run_fasttree(msa, fasttree_path='FastTree'):
    """
    Run FastTree on a protein multiple sequence alignment

    :param msa: file path to protein multiple sequence alignment in multifasta format
    :param phyml_path: path to phyml executable
    :return: path to the tree file
    """
    tree_path = msa + '.nw'
    command = [fasttree_path, '-out', tree_path, msa]

    logging.debug('Running FastTree: {}'.format(' '.join(command)))
    subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    return tree_path


def phylogenetic_tree_to_cluster_format(tree, pairwise_estimates):
    """
    Convert a phylogenetic tree to a 'cluster' data structure as in ``fastcluster``.
    The first two columns indicate the nodes that are joined by the relevant node,
    the third indicates the distance (calculated from branch lengths in the case of a phylogenetic tree)
    and the fourth the number of leaves underneath the node.

    Example of the data structure (output from ``fastcluster``):
    [[   3.            7.            4.26269776    2.        ]
     [   0.            5.           26.75703595    2.        ]
     [   2.            8.           56.16007598    2.        ]
     [   9.           12.           78.91813609    3.        ]
     [   1.           11.           87.91756528    3.        ]
     [   4.            6.           93.04790855    2.        ]
     [  14.           15.          114.71302639    5.        ]
     [  13.           16.          137.94616373    8.        ]
     [  10.           17.          157.29055403   10.        ]]

    :param tree: newick tree file
    :return: clustering data structure, pairwise distances dictionary
    """
    id_map = {pairwise_estimates.index[i]: i for i in range(len(pairwise_estimates))}
    t = Tree(tree)

    # midpoint rooting
    midpoint = t.get_midpoint_outgroup()
    t.set_outgroup(midpoint)

    # algorithm for getting cluster data structure
    n = len(id_map)
    out = []
    pairwise_distances = {}
    for node in t.traverse('postorder'):
        if node.is_leaf():
            node.name = id_map[node.name]
            id_map[node.name] = node.name  # add identity map for renamed nodes to id_map for line below
            pairwise_distances[node.name] = {id_map[x.name]: node.get_distance(x) for x in t.get_leaves()}
        else:
            node.name = n
            n += 1
            children = node.get_children()
            out.append(
                [children[0].name, children[1].name, children[0].get_distance(children[1]), len(node.get_leaves())])
    return np.array(out), pairwise_distances



