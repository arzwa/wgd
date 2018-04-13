"""
`wgd` WGD analysis in python

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
"""
# TODO: other outlier detection approaches
# TODO: alignment stripping as in Vanneste 2013
# TODO: add functionality to use custom alignments

# IMPORTS
from .codeml import Codeml
from .alignment import MSA, multiple_sequence_aligment_nucleotide
from .alignment import get_pairwise_nucleotide_alignments, \
    pairwise_alignment_stats
from .utils import process_gene_families, get_sequences, write_fasta
from .phy import run_phyml, phylogenetic_tree_to_cluster_format, run_fasttree, \
    average_linkage_clustering
from operator import itemgetter
from joblib import Parallel, delayed
import numpy as np
import pandas as pd
import os
import logging
import asyncio
import subprocess


# HELPER FUNCTIONS -------------------------------------------------------------
def _weighting(pairwise_estimates, msa=None, method='alc'):
    """
    Wrapper for different weighting methods.
    The fastest method is average linkage clustering based on Ks values,
    Other methods use phylogenetic trees (phyml or fasttree)

    :param pairwise_estimates: results dictionary from
        :py:meth:`codeml.Codeml.run()` output
    :param msa: protein multiple sequence alignment file path (optional,
        not necessary for method=``alc``)
    :param method: method, one of ``alc|phyml|fasttree``
    :return: clustering data structure and pairwise leaf distances
        (not if method=``alc``)
    """
    if pairwise_estimates is None:
        return None, None, None

    if pairwise_estimates['Ks'].shape[0] < 2:
        return None, None, None

    pairwise_distances = None
    tree_path = None

    if method == 'phyml':
        # PhyML tree construction
        logging.debug('Constructing phylogenetic tree with PhyML')
        tree_path = run_phyml(msa)
        clustering, pairwise_distances = phylogenetic_tree_to_cluster_format(
                tree_path, pairwise_estimates['Ks'])

    elif method == 'fasttree':
        # FastTree tree construction
        logging.debug('Constructing phylogenetic tree with FastTree')
        tree_path = run_fasttree(msa)
        clustering, pairwise_distances = phylogenetic_tree_to_cluster_format(
                tree_path, pairwise_estimates['Ks'])

    else:
        # Average linkage clustering based on Ks
        logging.debug('Performing average linkage clustering on Ks values.')
        clustering = average_linkage_clustering(pairwise_estimates['Ks'])

    logging.debug('Clustering used for weighting: \n{}'.format(str(clustering)))
    return clustering, pairwise_distances, tree_path


def _calculate_weights(clustering, pairwise_estimates, pairwise_distances=None):
    """
    This is a patch for weight calculation in the pairwise approach

    :param clustering: clustering array
    :param pairwise_estimates: pairwise Ks estimates array
    :param pairwise_distances: pairwise distances array
    :return: data frame
    """
    # None -> None
    if pairwise_estimates is None or clustering is None:
        return None

    # process the clustering structure to get weights
    leaves = pairwise_estimates.shape[0]
    nodes = {i: [i] for i in range(leaves)}
    weights = {}
    out = set()

    for x in range(clustering.shape[0]):
        node_1, node_2, distance = \
            clustering[x, 0], clustering[x, 1], clustering[x, 2] * 2
        grouping_node = leaves + x
        nodes[grouping_node] = nodes[node_1] + nodes[node_2]
        weight = 1 / (len(nodes[node_1]) * len(nodes[node_2]))

        # get for every pair the weight with outliers included and detect
        # outliers
        for i in nodes[node_1]:
            for j in nodes[node_2]:
                if pairwise_distances:
                    distance = pairwise_distances[i][j]
                pair = '_'.join(sorted([
                    pairwise_estimates.index[j], pairwise_estimates.index[i]
                ]))
                weights[pair] = {
                    'WeightOutliersIncluded': 0,
                    'WeightOutliersExcluded': 0,
                    'Outlier': 'TRUE',
                    'Distance': distance
                }
                weights[pair]['WeightOutliersIncluded'] = weight
                if pairwise_estimates.iloc[i, j] > 5:
                    out.add(grouping_node)

    # reweigh after outlier removal
    for node in out:
        nodes.pop(node)

    for x in range(clustering.shape[0]):
        node_1, node_2 = clustering[x, 0], clustering[x, 1]
        grouping_node = leaves + x
        if grouping_node in nodes and node_1 in nodes and node_2 in nodes:
            weight = 1 / (len(nodes[node_1]) * len(nodes[node_2]))
            for i in nodes[node_1]:
                for j in nodes[node_2]:
                    pair = '_'.join(sorted([
                        pairwise_estimates.index[j], pairwise_estimates.index[i]
                    ]))
                    weights[pair]['WeightOutliersExcluded'] = weight
                    weights[pair]['Outlier'] = 'FALSE'

    return pd.DataFrame.from_dict(weights, orient='index')


def _calculate_weighted_ks(clustering, pairwise_estimates,
                           pairwise_distances=None, family_id=None):
    """
    Calculate weighted Ks, Kn and w values following the procedure as outlined
    in Vanneste et al. (2013)

    :param clustering: clustering results as produced with
        :py:func:`_weighting`, can be by average linkage clustering or
        phylogenetic means (phyml, fasttree)
    :param pairwise_estimates: pairwise Ks, Kn and w estimates as produced by
        :py:meth:`codeml.Codeml.run_codeml`
    :return: a pandas data frame of weighted Ks, Kn and w values.
    """
    # None -> None
    if pairwise_estimates is None or clustering is None:
        return None

    # process the clustering structure to get weights
    leaves = pairwise_estimates['Ks'].shape[0]
    nodes = {i: [i] for i in range(leaves)}
    array = []
    out = set()

    for x in range(clustering.shape[0]):
        node_1, node_2, distance = clustering[x, 0], clustering[x, 1], \
                                   clustering[x, 2] * 2
        grouping_node = leaves + x
        nodes[grouping_node] = nodes[node_1] + nodes[node_2]
        weight = 1 / (len(nodes[node_1]) * len(nodes[node_2]))

        for i in nodes[node_1]:
            for j in nodes[node_2]:
                if pairwise_distances:
                    distance = pairwise_distances[i][j]
                array.append([
                    pairwise_estimates['Ks'].index[i],
                    pairwise_estimates['Ks'].index[j],
                    family_id.split("_")[-1],
                    weight,
                    pairwise_estimates['Ks'].iloc[i, j],
                    pairwise_estimates['Ka'].iloc[i, j],
                    pairwise_estimates['Omega'].iloc[i, j],
                    distance
                ])

                if pairwise_estimates['Ks'].iloc[i, j] > 5:
                    out.add(grouping_node)

    df1 = pd.DataFrame(array, columns=['Paralog1', 'Paralog2', 'Family',
                                       'WeightOutliersIncluded',
                                       'Ks', 'Ka', 'Omega', 'Distance'])

    # reweigh after outlier removal
    for node in out:
        nodes.pop(node)

    reweight = []
    for x in range(clustering.shape[0]):
        node_1, node_2 = clustering[x, 0], clustering[x, 1]
        grouping_node = leaves + x
        if grouping_node in nodes and node_1 in nodes and node_2 in nodes:
            weight = 1 / (len(nodes[node_1]) * len(nodes[node_2]))
            for i in nodes[node_1]:
                for j in nodes[node_2]:
                    reweight.append([pairwise_estimates['Ks'].index[i],
                                     pairwise_estimates['Ks'].index[j],
                                     weight, 'FALSE'])

    df2 = pd.DataFrame(reweight, columns=['Paralog1', 'Paralog2',
                                          'WeightOutliersExcluded', 'Outlier'])
    data_frame = pd.merge(df1, df2, how='outer', on=['Paralog1', 'Paralog2'])
    data_frame['Outlier'] = data_frame['Outlier'].fillna('TRUE')

    return data_frame


def add_alignment_stats(df, stats, l, l_stripped):
    """
    Add alignment statistics to the data frame

    :param df: pandas data frame
    :param stats: stats dict see :py:func:`alignment.pairwise_alignment_stats`
    :param l: alignment length
    :param l_stripped: stripped alignment length
    :return: data frame
    """
    identity = []
    coverage = []
    for row in df.index:
        paralogs = sorted([df.loc[row]['Paralog1'], df.loc[row]['Paralog2']])
        identity.append(stats[paralogs[0]][paralogs[1]][0])
        coverage.append(stats[paralogs[0]][paralogs[1]][1])
    df['AlignmentIdentity'] = pd.Series(identity, index=df.index)
    df['AlignmentCoverage'] = pd.Series(coverage, index=df.index)
    df['AlignmentLength'] = pd.Series([l] * len(df.index), index=df.index)
    df['AlignmentLengthStripped'] = pd.Series([l_stripped] * len(df.index),
                                              index=df.index)
    return df


# ANALYSE WHOLE FAMILY ---------------------------------------------------------
def analyse_family(
        family_id, family, nucleotide, tmp='./', muscle='muscle',
        codeml='codeml', prank='prank', preserve=False, times=1,
        min_length=100, method='alc', aligner='muscle',
        output_dir='./out'
):
    """
    Wrapper function for the analysis of one paralog family. Performs alignment
    with :py:meth:`alignment.MSA.run_aligner` and codeml analysis with
    :py:meth:`codeml.Codeml.run_codeml`. Subsequently also clustering with
    :py:func:`_average_linkage_clustering` is performed and weighted Ks, Kn and
    w values are calculated using :py:func:`_calculate_weighted_ks`.

    :param family_id: gene family id
    :param family: dictionary with sequences of paralogs
    :param nucleotide: nucleotide (CDS) sequences dictionary
    :param tmp: tmp directory
    :param muscle: muscle path
    :param codeml: codeml path
    :param prank: prank path
    :param preserve: preserve intermediate files
    :param times: number of times to perform ML estimation of Ks, Ka and omega
        values
    :param min_length: minimum length of the stripped multiple sequence
        alignment
    :param method: weighting method, from fast to slow: ``alc``, ``fasttree``,
        ``phyml``
    :param aligner: alignment program
    :return: ``csv`` file with results for the paralog family of interest
    """
    # pre-processing -----------------------------------------------------------
    if os.path.isfile(os.path.join(tmp, family_id + '.Ks')):
        logging.info('Found {}.Ks in tmp directory, will use this'
                     ''.format(family_id))
        return

    if len(list(family.keys())) < 2:
        logging.debug("Skipping singleton gene family {}.".format(family_id))
        return

    logging.info('Performing analysis on gene family {}'.format(family_id))

    # multiple sequence alignment ----------------------------------------------
    align = MSA(muscle=muscle, tmp=tmp, prank=prank)
    if aligner == 'prank':
        logging.debug('Aligner is prank, will perform codon alignment')
        family = {k: nucleotide[k] for k in family.keys() if len(nucleotide[k])
                  % 3 == 0}

    logging.debug('Performing multiple sequence alignment ({0}) on gene family '
                  '{1}.'.format(aligner, family_id))
    msa_path_protein = align.run_aligner(
            family, file=family_id, aligner=aligner)
    msa_out = multiple_sequence_aligment_nucleotide(
            msa_path_protein, nucleotide, min_length=min_length, aligner=aligner
    )
    if not msa_out:
        logging.warning('Did not analyze gene family {0}, stripped MSA too '
                        'short (< {1}).'.format(family_id, min_length))
        return

    else:
        msa_path, l, l_stripped, alignment_stats = msa_out

    # Calculate Ks values (codeml) ---------------------------------------------
    codeml = Codeml(codeml=codeml, tmp=tmp, id=family_id)
    logging.debug('Performing codeml analysis on gene family {}'
                  ''.format(family_id))
    results_dict, codeml_out = codeml.run_codeml(os.path.basename(msa_path),
                                                 preserve=preserve, times=times)

    # Calculate weights according to method ------------------------------------
    clustering, pairwise_distances, tree_path = _weighting(
            results_dict, msa=msa_path_protein, method=method)

    if clustering is not None:
        out = _calculate_weighted_ks(
                clustering, results_dict, pairwise_distances=pairwise_distances,
                family_id=family_id)
        out = add_alignment_stats(out, alignment_stats, l, l_stripped)
        out.to_csv(os.path.join(tmp, family_id + '.Ks'))

    # preserve or remove data --------------------------------------------------
    if preserve:
        subprocess.run(['mv', msa_path, os.path.join(output_dir, 'msa')])
        subprocess.run(['mv', codeml_out, os.path.join(output_dir, 'codeml')])
        if tree_path:
            subprocess.run(['mv', tree_path, os.path.join(output_dir, 'trees')])
    else:
        subprocess.run(['rm', msa_path])
        subprocess.run(['rm', codeml_out])
        if tree_path:
            subprocess.run(['rm', tree_path])


# PAIRWISE ANALYSIS ------------------------------------------------------------
def analyse_family_pairwise(
        family_id, family, nucleotide, tmp='./', muscle='muscle',
        codeml='codeml', prank='prank', preserve=False, times=1,
        min_length=100, method='alc', aligner='muscle',
        output_dir='./out'
):
    """
    Perform Ks analysis for one gene family using teh pairwise approach. The
    pairwise analysis approach is the one most commonly used in the literature,
    e.g. Vanneste et al. (2013), Vanneste et al. (2015), Barker et al. (2008),
    Li et al. (2015), Maere et al. (2005) and many more. Implementation details
    differ for many publications, and their are several implementation choices
    that are non-trivial. A main point of importance is the weighting of
    multiple pairwise Ks estimates based on a phylogeny (or proxy thereof), with
    some authors that do not weight (naive pairwise Ks distributions). Some
    authors use phylogentic trees, oter use hierarchical clustering for the
    weighting. Lastly some authors perform reweighting after outlier removal,
    others don't. Here the approach of Vanneste et al. (2013) and later papers
    is largely followed:

    (1) A multiple sequence alignment at the protein level is constructed for a
    full paralogous family.

    (2) For every pair of sequences in the multiple sequence alignment (i.e.
    n*(n-1)/2 pairs for a family of n members), strip all gaps in the pairwise
    alignment.

    (3) Then back-translate the sequences to obtain a codon alignment. Discard
    the sequence pair if the alignment is too short (min_length, default=100).

    (4) Use codeml to perform ML estimation of Ks, Ka and Ka/Ks. Here the
    empirical codon frequencies of the pairwise alignment with the F3x4 method.

    (5) A matrix is constructed with all pairwise estimates. Entries that do not
    have a valid Ks estimate (because the alignment was too short) are removed
    from the matrix. This is done by iteratively removing the columns (and rows)
    with most NaN values. These sequences are also removed from the aligment.

    (6) A weighting method is applied. Either a phylogenetic tree is constructed
    for the faily or the Pairwise Ks matrix is clustered using average linkage
    clustering. Every Ks estimate for a duplication node is then weighted by the
    number of estimates that are present for that node such that the different
    estimates for the same duplication event result in  weight of 1 in total.

    (7) Outliers are removed, where outliers are naively defined as pairs with
    Ks estimates > 5. Outlier removal essentially consists in putting the weight
    of these pairs to 0 and recalculating the weights for the remaining pairs in
    the family.

    :param family_id: gene family ID
    :param family: gene family sequence dictionary
    :param nucleotide: nucleotide sequences
    :param tmp: tmp directory
    :param muscle: muscle executable
    :param codeml: codeml executable
    :param prank: prank executable
    :param preserve: preserve codeml, alignment and tree results
    :param times: number of times to perform codeml estimation
    :param min_length: minimum gap-stripped alignment length to consider
    :param method: weighting method
    :param aligner: alignment method
    :param output_dir: output directory
    :return: nada
    """
    # pre-processing -----------------------------------------------------------
    if os.path.isfile(os.path.join(tmp, family_id + '.Ks')):
        logging.info(
                'Found {}.Ks in tmp directory, will use this'.format(family_id)
        )
        return

    if len(list(family.keys())) < 2:
        logging.debug("Skipping singleton gene family {}.".format(family_id))
        return

    logging.info('Performing analysis on gene family {}'.format(family_id))

    # multiple sequence alignment ----------------------------------------------
    align = MSA(muscle=muscle, tmp=tmp, prank=prank)
    if aligner == 'prank':
        logging.debug('Aligner is prank, will perform codon alignment')
        family = {
            k: nucleotide[k] for k in family.keys()
            if len(nucleotide[k]) % 3 == 0
        }

    logging.debug(
            'Performing multiple sequence alignment ({0}) on gene family '
            '{1}.'.format(aligner, family_id)
    )
    msa_path_protein, msa_dict = align.run_aligner(
            family, file=family_id, aligner=aligner, return_dict=True)
    pairwise_alignments = get_pairwise_nucleotide_alignments(
            msa_path_protein, nucleotide, min_length=min_length, aligner=aligner
    )

    # for every pair perform codeml analysis -----------------------------------
    codeml_out_string = ''
    family_dict = {}
    i = 1
    ks_mat = {g: {h: np.nan for h in family.keys()} for g in family.keys()}
    for pair in pairwise_alignments:
        pairwise_msa, g1, g2, l_stripped, l_unstripped, stats = pair
        codeml_ = Codeml(codeml=codeml, tmp=tmp, id=family_id + '.' + str(i))
        logging.debug(
                'Performing codeml analysis on gene family {}'.format(family_id)
        )
        results_dict, codeml_out = codeml_.run_codeml(
                os.path.basename(pairwise_msa), preserve=preserve, times=times
        )

        # remove pairwise alignment
        os.remove(pairwise_msa)

        # preserve codeml
        if preserve:
            with open(codeml_out, 'r') as f:
                codeml_out_string += f.read()
        os.remove(codeml_out)

        # store results
        family_dict['_'.join(sorted([g1, g2]))] = {
            'Paralog1': g1, 'Paralog2': g2, 'Family': family_id,
            'Ks': results_dict['Ks'][g1][g2],
            'Ka': results_dict['Ka'][g1][g2],
            'Omega': results_dict['Omega'][g1][g2],
            'AlignmentLengthStripped': l_stripped,
            'AlignmentLength': l_unstripped,
            'AlignmentCoverage': stats[1],
            'AlignmentIdentity': stats[0],
        }

        # also keep a Ks distance matrix
        ks_mat[g1][g2] = results_dict['Ks'][g1][g2]
        ks_mat[g2][g1] = results_dict['Ks'][g1][g2]
        ks_mat[g1][g1] = 0
        ks_mat[g2][g2] = 0

    ks_mat = pd.DataFrame.from_dict(ks_mat)
    logging.debug('Ks matrix: \n{}'.format(ks_mat))

    family_df = pd.DataFrame.from_dict(family_dict, orient='index')

    # filter the Ks matrix -----------------------------------------------------
    # As in the pairwise method, not every pair necessarily has a Ks estimate
    # (since filtering is on pairs within families, not entire families). Some
    # genes should be filtered out, as they might occur in pairs which have a
    # ks estimate and pairs which have not got one. In principle, this filtering
    # is not necessary, unless the weighting method is average linkage
    # clustering, in which case an incomplete distance matrix would be attained.
    # The strategy here is to remove iteratively the column/row with the most NA
    # entries until no more NA entries are in the Ks matrix. This will give the
    # most parsimonious filtered matrix
    while sum(ks_mat.isnull().sum()) != 0:
        logging.debug('Filtering matrix')
        ks_mat.drop(ks_mat.isnull().sum().idxmax(), axis=0, inplace=True)
        ks_mat.drop(ks_mat.isnull().sum().idxmax(), axis=1, inplace=True)

    logging.debug('Filtered Ks matrix: \n{}'.format(ks_mat))

    # filter the alignment -----------------------------------------------------
    # remove the sequences that were filtered out from the Ks matrix from the
    # alignment as well.
    to_del = [k for k in msa_dict.keys() if k not in ks_mat.index]
    for k in to_del:
        del msa_dict[k]

    write_fasta(msa_dict, msa_path_protein)

    # weighting ----------------------------------------------------------------
    # TODO: when only two family members, this is a bit of an overkill, but
    # otherwise the distance between the gene pair is not computed correctly.
    # Also, apparently PhyML breaks on families of two members.
    clustering, pairwise_distances, tree_path = _weighting(
            {'Ks': ks_mat}, msa=msa_path_protein, method=method
    )
    weights = _calculate_weights(clustering, ks_mat, pairwise_distances)
    pd.merge(
            family_df, weights, left_index=True, right_index=True, how='outer'
    ).to_csv(os.path.join(tmp, family_id + '.Ks'))

    # preserve -----------------------------------------------------------------
    if preserve:
        base = os.path.basename(msa_path_protein)
        os.rename(msa_path_protein, os.path.join(output_dir, 'msa', base))
        if tree_path:
            base = os.path.basename(tree_path)
            os.rename(tree_path, os.path.join(output_dir, 'trees', base))
        codeml_path = os.path.join(output_dir, 'codeml', family_id + '.codeml')
        with open(codeml_path, 'w') as f:
            f.write(codeml_out_string)
    else:
        os.remove(msa_path_protein)
        if tree_path:
            os.remove(tree_path)

    return


def ks_analysis_one_vs_one(
        nucleotide_sequences, protein_sequences,
        gene_families, tmp_dir='./tmp',
        output_dir='./ks.out', muscle_path='muscle',
        codeml_path='codeml', prank_path='prank',
        aligner='muscle', preserve=True, times=1,
        n_threads=4, async=False, min_length=100
):
    """
    Calculate a Ks distribution for one vs. one orthologs.

    :param nucleotide_sequences: sequence dictionary
    :param protein_sequences: protein sequence dictionary
    :param paralogs: file with paralog families
    :param tmp_dir: tmp directory
    :param output_dir: output directory
    :param muscle_path: path to muscle executable
    :param codeml_path: path to codeml executable
    :param prank_path: path to prank executable
    :param preserve: preserve intermediate results (muscle, codeml)
    :param times: number of times to perform codeml analysis
    :param ignore_prefixes: ignore prefixes in paralog/gene family file
        (e.g. in ath|AT1G45000, ath| will be ignored)
    :param async: use asyncio library for parallelization
    :param min_length: minimum MSA length
    :param aligner: aligner to use (muslce|prank)
    :param n_threads: number of CPU cores to use
    :return: data frame
    """
    # Filter families with one vs one orthologs for the species of interest. ---
    gene_families = process_gene_families(gene_families, ignore_prefix=False)
    protein = get_sequences(gene_families, protein_sequences)

    # preserve intermediate data if asked --------------------------------------
    if preserve:
        # msa and codeml results
        if not os.path.isdir(os.path.join(output_dir, 'msa')):
            os.mkdir(os.path.join(output_dir, 'msa'))
        if not os.path.isdir(os.path.join(output_dir, 'codeml')):
            os.mkdir(os.path.join(output_dir, 'codeml'))

    # start analysis -----------------------------------------------------------
    logging.info('Started analysis in parallel')
    if async:
        loop = asyncio.get_event_loop()
        tasks = [loop.run_in_executor(
                None, analyse_family, family, protein[family],
                nucleotide_sequences, tmp_dir, muscle_path,
                codeml_path, prank_path, preserve, times,
                min_length, 'alc', aligner, output_dir)
            for family in protein.keys()]
        loop.run_until_complete(asyncio.gather(*tasks))

    else:
        Parallel(n_jobs=n_threads)(
                delayed(analyse_family)(
                        family, protein[family],
                        nucleotide_sequences,
                        tmp_dir,
                        muscle_path,
                        codeml_path,
                        prank_path,
                        preserve,
                        times,
                        min_length,
                        'alc',
                        aligner,
                        output_dir) for family
                in protein.keys())
    logging.info('Analysis done')

    logging.info('Making results data frame')
    results_frame = pd.DataFrame(columns=['Paralog1', 'Paralog2', 'Family',
                                          'WeightOutliersIncluded', 'Ks', 'Ka',
                                          'Omega'])

    # count the number of analyzed pairs ---------------------------------------
    counts = 0
    for f in os.listdir(tmp_dir):
        if f[-3:] == '.Ks':
            counts += 1
            df = pd.read_csv(os.path.join(tmp_dir, f), index_col=0)
            results_frame = pd.concat([results_frame, df])
    results_frame.index = list(range(len(results_frame.index)))

    logging.info('Removing tmp directory')
    os.system('rm -r {}'.format(tmp_dir))

    # rename the index of the data_frame to gene1_gene2 (alphabetically) -------
    new_index = results_frame[['Paralog1', 'Paralog2']].apply(
            lambda x: '_'.join(sorted(x)), axis=1)
    results_frame.index = new_index

    return results_frame


def ks_analysis_paranome(
        nucleotide_sequences, protein_sequences, paralogs,
        tmp_dir='./tmp', output_dir='./ks.out',
        muscle_path='muscle', codeml_path='codeml',
        prank_path='prank', preserve=True, times=1,
        ignore_prefixes=False, n_threads=4, async=False,
        min_length=100, method='alc', aligner='muscle',
        pairwise=False, max_pairwise=10000
):
    """
    Calculate a Ks distribution for a whole paranome.

    :param nucleotide_sequences: sequence dictionary
    :param protein_sequences: protein sequence dictionary
    :param paralogs: file with paralog families
    :param tmp_dir: tmp directory
    :param output_dir: output directory
    :param muscle_path: path to muscle executable
    :param codeml_path: path to codeml executable
    :param prank_path: path to prank executable
    :param preserve: preserve intermediate results (muscle, codeml)
    :param times: number of times to perform codeml analysis
    :param ignore_prefixes: ignore prefixes in paralog/gene family file
        (e.g. in ath|AT1G45000, ath| will be ignored)
    :param async: use asyncio library for parallelization
    :param min_length: minimum MSA length
    :param method: method to use, from fast to slow: ``alc``, ``fasttree``,
        ``phyml``
    :param aligner: alignment program to use (muscle|prank)
    :param n_threads: number of CPU cores to use:
    :param pairwise: perform pairwise (instead of gene family wise) analysis
    :param max_pairwise: maximum number of pairwise combinations a family can
        have.
    :return: data frame
    """
    # ignore prefixes in gene families, since only one species -----------------
    paralogs = process_gene_families(paralogs, ignore_prefix=ignore_prefixes)
    protein = get_sequences(paralogs, protein_sequences)

    # preserve intermediate data if asked --------------------------------------
    if preserve:
        # msa and codeml results
        if not os.path.isdir(os.path.join(output_dir, 'msa')):
            os.mkdir(os.path.join(output_dir, 'msa'))
        if not os.path.isdir(os.path.join(output_dir, 'codeml')):
            os.mkdir(os.path.join(output_dir, 'codeml'))
        if method in ['fasttree', 'phyml']:
            if not os.path.isdir(os.path.join(output_dir, 'trees')):
                os.mkdir(os.path.join(output_dir, 'trees'))

    # sort family ids by family size -------------------------------------------
    sorted_families = sort_families_by_size(protein, pairwise, max_pairwise)

    # start analysis -----------------------------------------------------------
    logging.info('Started analysis in parallel (n_threads = {})'
                 ''.format(n_threads))
    if pairwise:
        analysis_function = analyse_family_pairwise
    else:
        analysis_function = analyse_family

    if async:
        loop = asyncio.get_event_loop()
        tasks = [loop.run_in_executor(
                None, analysis_function, family[0], protein[family[0]],
                nucleotide_sequences, tmp_dir, muscle_path,
                codeml_path, prank_path, preserve, times,
                min_length, method, aligner, output_dir)
            for family in sorted_families]
        loop.run_until_complete(asyncio.gather(*tasks))
    else:
        Parallel(n_jobs=n_threads)(delayed(analysis_function)(
                family[0], protein[family[0]],
                nucleotide_sequences,
                tmp_dir,
                muscle_path,
                codeml_path,
                prank_path,
                preserve,
                times,
                min_length,
                method,
                aligner,
                output_dir
        ) for family in sorted_families)
    logging.info('Analysis done')

    logging.info('Making results data frame')
    results_frame = pd.DataFrame(
            # is this necessary?
            columns=[
                'Paralog1', 'Paralog2', 'Family', 'WeightOutliersIncluded',
                'Ks', 'Ka', 'Omega'
            ]
    )

    # count the number of analyzed pairs and get data frame --------------------
    counts = 0
    for f in os.listdir(tmp_dir):
        if f[-3:] == '.Ks':
            counts += 1
            df = pd.read_csv(os.path.join(tmp_dir, f), index_col=0)
            results_frame = pd.concat([results_frame, df])
    results_frame.index = list(range(len(results_frame.index)))

    logging.info('Removing tmp directory')
    os.system('rm -r {}'.format(tmp_dir))

    # rename the index of the data_frame to gene1_gene2 (alphabetically) -------
    new_index = results_frame[['Paralog1', 'Paralog2']].apply(
            lambda x: '_'.join(sorted(x)), axis=1)
    results_frame.index = new_index

    return results_frame


def sort_families_by_size(
        families, pairwise=False, max_pairwise=10000
):
    """
    Returns a list of non-singleton gene families ordered by size and apply some
    filters

    :param families: nested gene family dictionary {family: {gene: sequence}}
    :param pairwise: pairwise analysis
    :param max_pairwise: maximum number of pairwise combinations a family may
        have
    :return: list of tuples [(family id, size)] sorted by size
    """
    sorted_families = []
    for k, v in families.items():
        if len(v.keys()) > 1:
            sorted_families.append((k, len(v.keys())))

    sorted_families = sorted(sorted_families, key=itemgetter(1), reverse=True)
    n_families = len(sorted_families)

    # filter for pairwise analysis
    if pairwise:
        sorted_families = [
            x for x in sorted_families
            if x[1] * (x[1] - 1) * 0.5 <= max_pairwise
        ]
        logging.warning(
                'Filtered out the {} largest gene families because n*(n-1)/2 > '
                '`max_pairwise`'.format(n_families - len(sorted_families))
        )

    return sorted_families
