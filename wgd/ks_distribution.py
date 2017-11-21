# Rewrite of the Ks distribution construction part of wgd
# Goals:
# 1. Simpler, assume standard data inputs generated with subtools of the wgd CLI.
# 2. One-to-one ortholog Ks distributions
# 3. No classes

# IMPORTS
from .codeml import Codeml
from .alignment import Muscle, multiple_sequence_aligment_nucleotide
from .utils import read_fasta, check_dirs, process_gene_families, get_sequences
from joblib import Parallel, delayed
import pandas as pd
import os
import fastcluster
import logging
import asyncio


# HELPER FUNCTIONS -----------------------------------------------------------------------------------------------------
def _average_linkage_clustering(pairwise_estimates):
    """
    Perform average linkage clustering using fastcluster.
    The first two columns of the output contain the node indices which are joined in each step.
    The input nodes are labeled 0, . . . , N - 1, and the newly generated nodes have the labels N, . . . , 2N - 2.
    The third column contains the distance between the two nodes at each step, ie. the
    current minimal distance at the time of the merge. The fourth column counts the
    number of points which comprise each new node.

    :param pairwise_estimates: dictionary with data frames with pairwise estimates of Ks, Ka and Ka/Ks
        (or at least Ks), as returned by :py:func:`analyse_family`.
    :return: average linkage clustering as performed with ``fastcluster.average``.
    """
    if pairwise_estimates is None:
        return None

    if pairwise_estimates['Ks'].shape[0] < 2:
        return None

    clustering = fastcluster.average(pairwise_estimates['Ks'])

    return clustering


def _calculate_weighted_ks(clustering, pairwise_estimates, family_id=None):
    """
    Calculate weighted Ks, Kn and w values following the procedure as outlined in Vanneste et al. (2013)

    :param clustering: average linkage clustering results as produced with :py:func:`_average_linkage_clustering`
    :param pairwise_estimates: pairwise Ks, Kn and w estimates as produced by :py:meth:`codeml.Codeml.run_codeml`
    :return: a pandas data frame of weighted Ks, Kn and w values.
    """
    # TODO: Check if reweighing is correct

    if pairwise_estimates is None:
        return None

    if clustering is None:
        return None

    leaves = pairwise_estimates['Ks'].shape[0]
    nodes = {i: [i] for i in range(leaves)}
    nodes_re = {i: [i] for i in range(leaves)}
    array = []
    out = set()

    for x in range(clustering.shape[0]):
        node_1, node_2 = clustering[x, 0], clustering[x, 1]
        grouping_node = leaves + x
        nodes[grouping_node] = nodes[node_1] + nodes[node_2]
        weight = 1 / (len(nodes[node_1]) * len(nodes[node_2]))

        for i in nodes[node_1]:
            for j in nodes[node_2]:
                array.append([
                    pairwise_estimates['Ks'].index[i],
                    pairwise_estimates['Ks'].index[j],
                    family_id.split("_")[-1],
                    weight,
                    pairwise_estimates['Ks'].iloc[i, j],
                    pairwise_estimates['Ka'].iloc[i, j],
                    pairwise_estimates['Omega'].iloc[i, j]
                ])

                if pairwise_estimates['Ks'].iloc[i, j] > 5:
                    out.add(grouping_node)

    df1 = pd.DataFrame(array, columns=['Paralog1', 'Paralog2', 'Family', 'WeightOutliersIncluded',
                                       'Ks', 'Ka', 'Omega'])

    # reweigh
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
                    reweight.append([pairwise_estimates['Ks'].index[i], pairwise_estimates['Ks'].index[j],
                                     weight, 'FALSE'])

    df2 = pd.DataFrame(reweight, columns=['Paralog1', 'Paralog2', 'WeightOutliersExcluded', 'Outlier'])
    data_frame = pd.merge(df1, df2, how='outer', on=['Paralog1', 'Paralog2'])
    data_frame['Outlier'] = data_frame['Outlier'].fillna('TRUE')

    return data_frame


def filter_species(results, s1, s2):
    """
    Filter out one vs one orthologs from codeml results for two species

    :param results:
    :param s1:
    :param s2:
    :return:
    """
    filtered = {}
    for key, val in results.items():
        for i in val.index:
            if not (i.split('|')[0] == s1 or i.split('|')[0] == s2):
                val = val.drop(i)
                val = val.drop(i, axis=1)
        filtered[key] = val
    return filtered


def analyse_family(family_id, family, nucleotide, tmp='./', muscle='muscle', codeml='codeml', preserve=False, times=1,
                   min_length=100):
    """
    Wrapper function for the analysis of one paralog family. Performs alignment with
    :py:meth:`alignment.Muscle.run_muscle` and codeml analysis with :py:meth:`codeml.Codeml.run_codeml`.
    Subsequently also clustering with :py:func:`_average_linkage_clustering` is performed and weighted Ks, Kn and w
    values are calculated using :py:func:`_calculate_weighted_ks`.

    :param family_id: gene family id
    :param family: dictionary with sequences of paralogs
    :param nucleotide: nucleotide (CDS) sequences dictionary
    :param tmp: tmp directory
    :param muscle: muscle path
    :param codeml: codeml path
    :param preserve: preserve intermediate files
    :return: ``csv`` file with results for the paralog family of interest
    """
    if os.path.isfile(os.path.join(tmp, family_id + '.Ks')):
        logging.info('Found {}.Ks in tmp directory, will use this'.format(family_id))
        return

    if len(list(family.keys())) < 2:
        logging.info("Skipping singleton gene family {}.".format(family_id))
        return

    logging.info('Performing analysis on gene family {}'.format(family_id))

    codeml = Codeml(codeml=codeml, tmp=tmp, id=family_id)
    muscle = Muscle(muscle=muscle, tmp=tmp)

    logging.debug('Performing multiple sequence alignment (MUSCLE) on gene family {}.'.format(family_id))
    msa_path = muscle.run_muscle(family, file=family_id)
    msa_path = multiple_sequence_aligment_nucleotide(msa_path, nucleotide, min_length=min_length)

    if not msa_path:
        logging.warning('Did not analyze gene family {0}, stripped MSA too short (< {1}).'.format(
            family_id, min_length))
        return

    # Calculate Ks values (codeml)
    logging.debug('Performing codeml analysis on gene family {}'.format(family_id))
    results_dict = codeml.run_codeml(msa_path, preserve=preserve, times=times)

    logging.debug('Performing average linkage clustering on Ks values for gene family {}.'.format(family_id))
    # Average linkage clustering based on Ks
    clustering = _average_linkage_clustering(results_dict)

    if clustering is not None:
        out = _calculate_weighted_ks(clustering, results_dict, family_id)
        out.to_csv(os.path.join(tmp, family_id + '.Ks'))


def analyse_family_ortholog(family_id, family, nucleotide, tmp='./', muscle='muscle', codeml='codeml',
                            preserve=False, times=1, min_length=100):
    """
    Wrapper function for the analysis of one paralog family. Performs alignment with
    :py:meth:`alignment.Muscle.run_muscle` and codeml analysis with :py:meth:`codeml.Codeml.run_codeml`.
    Subsequently also clustering with :py:func:`_average_linkage_clustering` is performed and weighted Ks, Kn and w
    values are calculated using :py:func:`_calculate_weighted_ks`.

    :param family_id: gene family id
    :param family: dictionary with sequences of paralogs
    :param nucleotide: nucleotide (CDS) sequences dictionary
    :param tmp: tmp directory
    :param muscle: muscle path
    :param codeml: codeml path
    :param preserve: preserve intermediate files
    :return: ``csv`` file with results for the paralog family of interest
    """
    if os.path.isfile(os.path.join(tmp, family_id + '.Ks')):
        logging.info('Found {}.Ks in tmp directory, will use this'.format(family_id))
        return

    if len(list(family.keys())) < 2:
        logging.info("Skipping singleton gene family {}.".format(family_id))
        return

    logging.info('Performing analysis on gene family {}'.format(family_id))

    codeml = Codeml(codeml=codeml, tmp=tmp, id=family_id)
    muscle = Muscle(muscle=muscle, tmp=tmp)

    logging.debug('Performing multiple sequence analysis (MUSCLE) on gene family {}.'.format(family_id))
    msa_path = muscle.run_muscle(family, file=family_id)
    msa_path = multiple_sequence_aligment_nucleotide(msa_path, nucleotide, min_length=min_length)

    if not msa_path:
        logging.warning('Did not analyze gene family {0}, stripped MSA too short (< {1}).'.format(
            family_id, min_length))
        return

    # Calculate Ks values (codeml)
    logging.debug('Performing codeml analysis on gene family {}'.format(family_id))
    results_dict = codeml.run_codeml(msa_path, preserve=preserve, times=times)

    # TODO: this is not necessary (overhead) but I just do it because it requires less programming
    logging.debug('Performing average linkage clustering on Ks values for gene family {}.'.format(family_id))
    # Average linkage clustering based on Ks
    clustering = _average_linkage_clustering(results_dict)

    if clustering is not None:
        out = _calculate_weighted_ks(clustering, results_dict, family_id)
        out.to_csv(os.path.join(tmp, family_id + '.Ks'))


def ks_analysis_one_vs_one(nucleotide_sequences, protein_sequences, gene_families, tmp_dir='./tmp',
                           output_dir='./ks.out', muscle_path='muscle', codeml_path='codeml',
                           preserve=True, times=1, n_cores=4, async=False, min_length=100):

    # Filter families with one vs one orthologs for the species of interest.
    gene_families = process_gene_families(gene_families, ignore_prefix=False)
    protein = get_sequences(gene_families, protein_sequences)

    # start analysis
    logging.info('Started analysis in parallel')
    if async:
        loop = asyncio.get_event_loop()
        tasks = [loop.run_in_executor(None, analyse_family_ortholog, family, protein[family],
                                      nucleotide_sequences, tmp_dir, muscle_path,
                                      codeml_path, preserve, times, min_length) for family in protein.keys()]
        loop.run_until_complete(asyncio.gather(*tasks))

    else:
        Parallel(n_jobs=n_cores)(delayed(analyse_family)(family, protein[family],
                                                         nucleotide_sequences, tmp_dir, muscle_path,
                                                         codeml_path, preserve, times,
                                                         min_length) for family in protein.keys())
    logging.info('Analysis done')

    # preserve intermediate data if asked
    if preserve:
        logging.info('Moving files (preserve=True)')
        if not os.path.isdir(os.path.join(output_dir, 'msa')):
            os.mkdir(os.path.join(output_dir, 'msa'))
        if not os.path.isdir(os.path.join(output_dir, 'codeml')):
            os.mkdir(os.path.join(output_dir, 'codeml'))
        os.system(" ".join(['mv', os.path.join(tmp_dir, '*.msa*'), os.path.join(output_dir, 'msa')]))
        os.system(" ".join(['mv', os.path.join(tmp_dir, '*.codeml'), os.path.join(output_dir, 'codeml')]))

    else:
        logging.info('Removing files (preserve=False)')
        os.system(" ".join(['rm', os.path.join(tmp_dir, '*.msa*')]))

    logging.info('Making results data frame')
    results_frame = pd.DataFrame(columns=['Paralog1', 'Paralog2', 'Family',
                                          'WeightOutliersIncluded', 'Ks', 'Ka', 'Omega'])

    # count the number of analyzed pairs
    counts = 0
    for f in os.listdir(tmp_dir):
        if f[-3:] == '.Ks':
            counts += 1
            df = pd.read_csv(os.path.join(tmp_dir, f), index_col=0)
            results_frame = pd.concat([results_frame, df])
    results_frame.index = list(range(len(results_frame.index)))

    logging.info('Removing tmp directory')
    os.system('rm -r {}'.format(tmp_dir))

    return results_frame


def ks_analysis_paranome(nucleotide_sequences, protein_sequences, paralogs, tmp_dir='./tmp', output_dir='./ks.out',
                         muscle_path='muscle', codeml_path='codeml', preserve=True, times=1,
                         ignore_prefixes=False, n_cores=4, async=False, min_length=100):
    """
    Calculate a Ks distribution for a whole paranome. Asyncio version

    :param nucleotide_sequences: CDS fasta file
    :param protein_sequences: protein fasta file
    :param paralogs: file with paralog families
    :param tmp_dir: tmp directory
    :param output_dir: output directory
    :param muscle_path: path to muscle executable
    :param codeml_path: path to codeml executable
    :param preserve: preserve intermediate results (muscle, codeml)
    :param times: number of times to perform codeml analysis
    :param ignore_prefixes: ignore prefixes in paralog/gene family file (e.g. in ath|AT1G45000, ath| will be ignored)
    :return: some files and plots in the chosen output directory
    """

    # ignore prefixes in gene families, since only one species
    paralogs = process_gene_families(paralogs, ignore_prefix=ignore_prefixes)
    protein = get_sequences(paralogs, protein_sequences)

    # start analysis
    logging.info('Started analysis in parallel (n_cores = {})'.format(n_cores))
    if async:
        loop = asyncio.get_event_loop()
        tasks = [loop.run_in_executor(None, analyse_family, family, protein[family],
                                      nucleotide_sequences, tmp_dir, muscle_path,
                                      codeml_path, preserve, times, min_length) for family in protein.keys()]

        loop.run_until_complete(asyncio.gather(*tasks))
    else:
        Parallel(n_jobs=n_cores)(delayed(analyse_family)(family, protein[family],
                                                         nucleotide_sequences, tmp_dir, muscle_path,
                                                         codeml_path, preserve, times,
                                                         min_length) for family in protein.keys())
    logging.info('Analysis done')

    # preserve intermediate data if asked
    if preserve:
        logging.info('Moving files (preserve=True)')
        if not os.path.isdir(os.path.join(output_dir, 'msa')):
            os.mkdir(os.path.join(output_dir, 'msa'))
        if not os.path.isdir(os.path.join(output_dir, 'codeml')):
            os.mkdir(os.path.join(output_dir, 'codeml'))
        os.system(" ".join(['mv', os.path.join(tmp_dir, '*.msa*'), os.path.join(output_dir, 'msa')]))
        os.system(" ".join(['mv', os.path.join(tmp_dir, '*.codeml'), os.path.join(output_dir, 'codeml')]))

    else:
        logging.info('Removing files (preserve=False)')
        os.system(" ".join(['rm', os.path.join(tmp_dir, '*.msa*')]))

    logging.info('Making results data frame')
    results_frame = pd.DataFrame(columns=['Paralog1', 'Paralog2', 'Family',
                                          'WeightOutliersIncluded', 'Ks', 'Ka', 'Omega'])

    # count the number of analyzed pairs
    counts = 0
    for f in os.listdir(tmp_dir):
        if f[-3:] == '.Ks':
            counts += 1
            df = pd.read_csv(os.path.join(tmp_dir, f), index_col=0)
            results_frame = pd.concat([results_frame, df])
    results_frame.index = list(range(len(results_frame.index)))

    logging.info('Removing tmp directory')
    os.system('rm -r {}'.format(tmp_dir))

    return results_frame
