#!/usr/bin/python3.5
"""
Arthur Zwaenepoel
"""
# TODO: maybe make one visualization for all Ks analyses using d3.js

# IMPORTS
from .codeml import Codeml
from .alignment import Muscle, multiple_sequence_aligment_nucleotide
from .utils import read_fasta, check_dirs, process_gene_families_ortho_mcl, get_sequences, translate_cds
from .modeling import mixture_model_bgmm, mixture_model_gmm, kernel_density_estimation, weighted_to_unweighted
from .html_report import write_report_ks
import pandas as pd
import numpy as np
import os
import re
import fastcluster
import logging
import plumbum as pb
import asyncio
import matplotlib
if not 'DISPLAY' in pb.local.env:
    matplotlib.use('Agg')  # use this backend when no X server
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)


def _get_paralogs(gene_family_dict, gene_pattern):
    """
    Get paralogs for a species from gene families

    :param gene_family_dict: gene family dictionary
    :param gene_pattern: gene pattern
    :return: dictionary with paralogs
    """
    if not gene_pattern:
        return gene_family_dict

    gene_pattern = re.compile(gene_pattern)
    paralog_dict = {}

    for key in gene_family_dict.keys():
        paralogs = [gene for gene in gene_family_dict[key] if gene_pattern.match(gene) is not None]
        if len(paralogs) > 1:
            paralog_dict[key] = paralogs

    return paralog_dict


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
    Calculate weighted Ks, Kn and w values following the procedure as outline in Vanneste et al. (2013)

    :param clustering: average linkage clustering results as produced with :py:func:`_average_linkage_clustering`
    :param pairwise_estimates: pairwise Ks, Kn and w estimates as produced by :py:meth:`codeml.Codeml.run_codeml`
    :return: a pandas data frame of weighted Ks, Kn and w values.
    """
    # TODO: Check if reweighting is correct

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

    # reweight
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


def analyse_family(family_id, family, nucleotide, tmp='./', muscle='muscle', codeml='codeml', preserve=False, times=1):
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
    msa_path = multiple_sequence_aligment_nucleotide(msa_path, nucleotide)

    if not msa_path:
        logging.warning('Did not analyze gene family {}, stripped MSA too short.'.format(family_id))
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


def _histogram_for_d3(data_frame, out_file='unweighted.csv'):
    """
    Generate a histogram for usage in d3 visualizations

    :param data_frame: data frame with Ks results as obtained using :py:class:`KsDistribution`
    :param out_file: output file name
    """
    df = pd.DataFrame(weighted_to_unweighted(data_frame))
    df.columns = ['Ks']
    df.to_csv(out_file)


# Ks DISTRIBUTION CLASS
class KsDistribution():
    """
    Class representing a Ks distribution and the data used for it.
    Note that Ks distribution analysis starts from CDS sequences and gene families (singe or multi-species).
    For performing all-`versus`-all Blastp, use :py:func:`wgd.mcl.all_v_all_blast`.
    For creating gene families from all-`versus`-all Blast results using MCL use :py:func:`wgd.mcl.run_mcl_ava`.
    """

    def __init__(self, species='species', gene_families=None, nucleotide=None, protein=None, gene_pattern=None,
                 muscle_path='muscle', codeml_path='codeml', tmp_dir='./', output_dir='./', results_file=None):
        """
        init for KsDistribution class.

        :param species: species id
        :param gene_families: gene families file (tab separated genes per line)
        :param nucleotide: nucleotide fasta file
        :param protein: protein fasta file
        :param gene_pattern: gene pattern (in case multi-species gene families provided)
        :param muscle_path: path to musce executable
        :param codeml_path: path to codeml executable
        :param tmp_dir: temporary directory
        :param output_dir: output directory
        :param results_file: results_file (use this if you want to load a previously computed data frame and perform
            downstream analysis such as mixture modeling/KDE analysis)
        """
        self.species = species
        self.gene_pattern = gene_pattern
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
        self.best_model_ = None
        self.aic_ = None
        if results_file:
            self.results = pd.read_csv(results_file, index_col=0)

    def __str__(self):
        return "Ks distribution analysis for: " + self.species

    def ks_analysis_pipeline(self, verbose=True, preserve=False, prompt=True, times=1):
        """
        Full analysis pipeline wrapper function, non-parallel implementation. (DEPRECATED)

        :return: Ks, Ka and w distributions for the whole paranome as csv files and histograms.
        """
        self.nucleotide = read_fasta(self.nucleotide_file, split_on_pipe=True)
        self.gf_dict = process_gene_families_ortho_mcl(self.gene_families)
        self.paralogs = _get_paralogs(self.gf_dict, self.gene_pattern)

        if not self.protein_file:
            self.protein = translate_cds(self.nucleotide)
        else:
            self.protein = get_sequences(self.paralogs, self.protein_file)

        # Initialize data frames
        output_dict = {
            'Ks': pd.DataFrame(columns=['gene_1', 'gene_2', 'Ks', 'weight']),
            'Kn': pd.DataFrame(columns=['gene_1', 'gene_2', 'Kn', 'weight']),
            'w': pd.DataFrame(columns=['gene_1', 'gene_2', 'w', 'weight'])
            }

        results_frame = pd.DataFrame(columns=['gene_1', 'gene_2', 'weight', 'Ks', 'Kn', 'w'])

        check_dirs(self.tmp_dir, self.output_dir, prompt, preserve)

        # Start Ks calculation pipeline

        counts = 0

        codeml = Codeml(tmp=self.tmp_dir)
        muscle = Muscle(muscle=self.muscle_path, tmp=self.tmp_dir)

        for family_ID in sorted(self.protein.keys()):
            if verbose:
                print('{:*^38}'.format(family_ID))

            # Make the multiple sequence alignment
            print("Performing multiple sequence alignment")
            family = self.protein[family_ID]
            msa_path = muscle.run_muscle(family, file=family_ID)
            msa_path = multiple_sequence_aligment_nucleotide(msa_path, self.nucleotide)

            if msa_path is None:
                print("{:>38}".format('Stripped MSA too short'))
                print("{:>38}".format('Did not add gene family {}'.format(family_ID)))
                continue

            print("Performing codeml analysis")

            # Calculate Ks values (codeml)
            results_dict = codeml.run_codeml(msa_path)

            # Average linkage clustering based on Ks
            clustering = _average_linkage_clustering(results_dict)

            if clustering is not None:

                out = _calculate_weighted_ks(clustering, results_dict)

                # Add to data frame
                if out is not None:
                    results_frame = pd.concat([results_frame, out])
                    counts += 1
                    if verbose:
                        print("{:>38}".format("Added gene family {}".format(family_ID)))
                else:
                    if verbose:
                        print("{:>38}".format("Did not add gene family {}".format(family_ID)))
            else:
                if verbose:
                    print("{:>38}".format("Did not add gene family {}".format(family_ID)))

            # If preserve, move tmp data to output directory
            if preserve:
                os.system('mv {0} {1}'.format(msa_path, os.path.join(self.output_dir, 'msa')))

            # Else, remove data from tmp
            else:
                if msa_path is not None:
                    os.remove(msa_path)

        # save the full results as csv
        results_frame.to_csv(os.path.join(self.output_dir, 'full_analysis.csv'))

        # Process the Ks distribution
        for key in ['Ks', 'Kn', 'w']:
            distribution = results_frame[['gene_1', 'gene_2', 'weight', key]]
            distribution.to_csv(os.path.join(self.output_dir, '{}_distribution.csv'.format(key)))
            distribution = distribution[distribution[key] <= 5]

            metric = key
            if metric == 'w':
                metric = '\omega'
            elif metric == 'Ks':
                metric = 'K_S'
            else:
                metric = 'K_N'

            # Make Ks plot, omitting Ks values <= 0.1 to avoid the incorporation of allelic and/or splice variants
            plt.figure()
            plt.title("{0} ${1}$ distribution".format(self.species, metric))
            plt.xlabel("Binned ${}$".format(metric))
            plt.hist(distribution[distribution[key] >= 0.1][key], bins=100,
                     weights=distribution[distribution[key] >= 0.1]['weight'],
                     histtype='stepfilled', color='#82c982')
            plt.xlim(0.1, max(distribution[key]))
            plt.savefig(os.path.join(self.output_dir, '{}_distribution.png'.format(key)), dpi=400)

        self.results = results_frame
        self.counts = counts

        return self

    def ks_analysis_parallel(self, check=True, preserve=True, times=1):
        """
        Parallel implementation for full analysis pipeline. Uses the asyncio library.

        :return: Ks, Ka and w distributions for the whole paranome as csv files and histograms.
        """
        logging.debug('Pre-processing sequence and gene family data')
        self.nucleotide = read_fasta(self.nucleotide_file, split_on_pipe=True)
        self.gf_dict = process_gene_families_ortho_mcl(self.gene_families)
        self.paralogs = _get_paralogs(self.gf_dict, self.gene_pattern)

        if not self.protein_file:
            logging.info('Translating CDS to protein (standard genetic code)')
            protein = translate_cds(self.nucleotide)
            self.protein = get_sequences(self.paralogs, protein)
        else:
            self.protein = get_sequences(self.paralogs, self.protein_file)

        if check:
            logging.debug('Checking directories (tmp, output)')
            check_dirs(self.tmp_dir, self.output_dir, prompt=True, preserve=preserve)

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
        results_frame = pd.DataFrame(columns=['Paralog1', 'Paralog2', 'Family',
                                              'WeightOutliersIncluded', 'Ks', 'Ka', 'Omega'])

        counts = 0
        for f in os.listdir(self.tmp_dir):
            if f[-3:] == '.Ks':
                counts += 1
                df = pd.read_csv(os.path.join(self.tmp_dir, f), index_col=0)
                results_frame = pd.concat([results_frame, df])

        results_frame.index = list(range(len(results_frame.index)))

        logging.info('Removing tmp directory')
        os.system('rm -r {}'.format(self.tmp_dir))
        results_frame.to_csv(os.path.join(self.output_dir, 'all.csv'))

        # Process the Ks distribution
        for key in ['Ks', 'Ka', 'Omega']:
            logging.info('Processing {} distribution'.format(key))

            distribution = results_frame[['Paralog1', 'Paralog2', 'WeightOutliersIncluded', key]]
            distribution.to_csv(os.path.join(self.output_dir, '{}_distribution.csv'.format(key)))
            distribution = distribution[distribution[key] <= 5]

            metric = key
            if metric == 'Omega':
                metric = 'ln(\omega)'
            elif metric == 'Ks':
                metric = 'K_S'
            else:
                metric = 'K_A'

            # Make Ks plot, omitting Ks values <= 0.1 to avoid the incorporation of allelic and/or splice variants
            logging.info('Generating png for {} distribution'.format(key))
            plt.figure()
            plt.title("{0} ${1}$ distribution".format(self.species, metric))
            plt.xlabel("Binned ${}$".format(metric))
            if key == 'Omega':
                plt.hist(np.log(distribution[key]), bins=100,
                         weights=distribution['WeightOutliersIncluded'],
                         histtype='stepfilled', color='#82c982')
            else:
                plt.hist(distribution[distribution[key] >= 0.1][key], bins=100,
                         weights=distribution[distribution[key] >= 0.1]['WeightOutliersIncluded'],
                         histtype='stepfilled', color='#82c982')
                plt.xlim(0.1, max(distribution[key]))
            plt.savefig(os.path.join(self.output_dir, '{}_distribution.png'.format(key)), dpi=400)

        self.results = results_frame
        self.counts = counts

        return self

    def mixture_modeling(self, method='bgmm', **kwargs):
        """
        Mixture modeling method. Uses by default a Bayesian Gaussian mixture modeling approach.

        :param kwargs: Keyword arguments for :py:func:`wgd.modeling.mixture_model_bgmm`
        :param method: Mixture modeling method (bayesian/gaussian)
        :return: sklearn.mixture models for different numbers of components
        """
        if self.results is None:
            raise AttributeError('No results attribute found, please run the analysis or load a data frame first.')

        else:
            if method == 'bgmm':
                self.mixture_models = mixture_model_bgmm(
                    self.results[['Paralog1', 'Paralog2', 'WeightOutliersIncluded', 'Ks']], output_dir=self.output_dir,
                    **kwargs
                )

            elif method == 'gmm':
                self.mixture_models, self.best_model_, self.aic_ = mixture_model_gmm(
                    self.results[['Paralog1', 'Paralog2', 'WeightOutliersIncluded', 'Ks']], output_dir=self.output_dir,
                    **kwargs
                )

            else:
                raise ValueError("No such mixture modeling method: {}".format(method))

        return self

    def kernel_density(self, **kwargs):
        """
        Kernel density estimation method.

        :param kwargs: Keyword arguments for :py:func:`wgd.modeling.kernel_density_estimation`
        :return: sklearn.KernelDensity objects for different bandwidths
        """
        if self.results is None:
            raise AttributeError('No results attribute found, please run the analysis or load a data frame first.')

        else:
            self.kde = kernel_density_estimation(
                self.results[['Paralog1', 'Paralog2', 'WeightOutliersIncluded', 'Ks']], output_dir=self.output_dir,
                **kwargs
            )

        return self

    def write_report(self, mixture, kde):
        """
        Method for writing a html report.

        :return: writes a `html` report of the performed analysis
        """
        write_report_ks(self.output_dir, self.species, self.gene_families,
                        self.nucleotide_file, self.protein_file, self.counts, mixture, kde)

        return self

    def load_distribution(self, csv_file, species='species'):
        """
        Method for loading a data frame with results.
        For visualization and modeling of precomputed results.

        :param csv_file: CSV file with at least following columns:
            ``'Paralog1','Paralog2','WeightOutliersIncluded','Ks'``
        :return: ``self``
        """
        data_frame = pd.read_csv(csv_file, index_col=0)
        self.results = data_frame
        self.species = species

    def histogram(self, out_file=None, bins=50, key='Ks', ks_range=(0.1, 5), color='#82c982', alpha=1,
                  histtype='stepfilled', rwidth=0.85, fig_size=(12,9)):
        """
        Plot histogram of a distribution in the results frame (Ks, Ka or omega).

        :param key: distribution to plot, usually one of ``'Ks'``, ``'Ka'`` or ``'Omega'``
        :return: a histogram of the distributions
        """
        if type(self.results) != pd.core.frame.DataFrame:
            raise AttributeError("No results data frame found, please run the analysis or load a data frame first.")

        for col in ['Paralog1', 'Paralog2', 'WeightOutliersIncluded', key]:
            if col not in self.results.columns:
                raise AttributeError("{} not found in columns of the results data frame (self.results)".format(col))

        distribution = self.results[['Paralog1', 'Paralog2', 'WeightOutliersIncluded', key]]
        distribution = distribution[distribution[key] <= ks_range[1]]
        distribution = distribution[distribution[key] >= ks_range[0]]

        metric = key
        if metric == 'Omega':
            metric = '\omega'
        elif metric == 'Ks':
            metric = 'K_S'
        else:
            metric = 'K_A'

        plt.figure(figsize=fig_size)
        plt.title("\\textit{{{0}}} ${1}$ distribution".format(self.species, metric))
        plt.xlabel("Binned ${}$".format(metric))

        plt.hist(distribution[key], bins=bins,
                 weights=distribution[distribution[key] >= 0.1]['WeightOutliersIncluded'],
                 histtype=histtype, color=color, alpha=alpha, rwidth=rwidth)
        plt.xlim(ks_range[0], ks_range[1])

        if out_file:
            plt.savefig(out_file, dpi=400)

        else:
            plt.show()
