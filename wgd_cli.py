#!/usr/bin/python3.5
"""
Arthur Zwaenepoel
"""
# TODO: separate subcommand for mixture modeling of Ks distributions
# TODO: this can than also include a mixture model + peak based paralog extraction tool?

import click
import logging
import sys
import os
import datetime
import pandas as pd
import uuid
import coloredlogs
from wgd.utils import translate_cds, read_fasta, write_fasta, get_one_v_one_orthologs_rbh, Genome


# CLI ENTRYPOINT -------------------------------------------------------------------------------------------------------
@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--verbose', type=click.Choice(['silent', 'info', 'debug']),
              default='info', help="Verbosity level, default = info.")
def cli(verbose):
    """
    Welcome to the wgd command line interface!

    \b
                           _______
                           \\  ___ `'.
           _     _ .--./)   ' |--.\\  \\
     /\\    \\\\   ///.''\\\\    | |    \\  '
     `\\\\  //\\\\ //| |  | |   | |     |  '
       \\`//  \\'/  \\`-' /    | |     |  |
        \\|   |/   /("'`     | |     ' .'
         '        \\ '---.   | |___.' /'
                   /'""'.\\ /_______.'/
                  ||     ||\\_______|/
                  \\'. __//
                   `'---'
    \b
    Arthur Zwaenepoel - 2017
    """
    coloredlogs.install(fmt='%(asctime)s: %(levelname)s\t%(message)s', level=verbose.upper(), stream=sys.stdout)
    pass


# BLAST AND MCL --------------------------------------------------------------------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--cds', is_flag=True, help='Sequences are CDS.')
@click.option('--mcl', is_flag=True, help='Perform MCL clustering.')
@click.option('--one_v_one', is_flag=True, help='Get one vs. one orthologs')
@click.option('--sequences','-s', default=None,
              help='Input fasta files, as a comma separated string (e.g. x.fasta,y.fasta,z.fasta).')
@click.option('--species_ids','-id', default=None,
              help='Species identifiers for respective input sequence files, as a comma separated '
                   'string (e.g. x,y,z). (optional)')
@click.option('--blast_results','-b', default=None,
              help='Input precomputed tab separated blast results.')
@click.option('--inflation_factor', '-I', default=2.0,
              help="Inflation factor for MCL clustering. (Default = 2)")
@click.option('--eval_cutoff', '-e', default=1e-10,
              help="E-value cut-off for Blast results (Default = 1e-10)")
@click.option('--output_dir','-o', default='wgd.blast.out',
              help='Output directory.')
def blast(cds, mcl, one_v_one, sequences, species_ids, blast_results, inflation_factor, eval_cutoff, output_dir):
    """
    Perform all-vs.-all Blastp (+ MCL) analysis.

    Example 1 - whole paranome delineation:

        wgd blast --cds --mcl -s thorny_devil.fasta -o thorny_devil_blast_out

    Example 2 - one vs. one ortholog delineation:

        wgd blast --cds --one_v_one -s equus_ferus.fasta,ursus_arctos.fasta -id horse,bear -e 1e-8 -o bear_horse_out
    """
    # lazy imports
    from wgd.blast_mcl import all_v_all_blast, run_mcl_ava_2, ava_blast_to_abc_2

    if not sequences and not blast_results:
        logging.error('No sequences nor blast results provided! Please use the --help flag for usage instructions.')
        return

    if not os.path.exists(output_dir):
        logging.info('Output directory: {} does not exist, will make it.'.format(output_dir))
        os.mkdir(output_dir)

    if not blast_results:
        sequence_files = sequences.strip().split(',')

        if len(sequence_files) != 1 and not one_v_one:
            logging.error('Please provide only one fasta file for whole paranome all-vs-all blast')
            return

        if len(sequence_files) != 2 and one_v_one:
            logging.error('Please provide (only) two fasta files for one-vs-one ortholog finding')
            return

        if species_ids:
            ids = species_ids.strip().split(',')
            if len(ids) != len(sequence_files):
                logging.error('Number of species identifiers ({0}) does not match number of provided sequence '
                              'files ({1}).'.format(len(ids), len(sequence_files)))
        else:
            ids = [''] * len(sequence_files)

        if cds:
            logging.info("CDS sequences provided, will first translate.")

        protein_sequences = []
        for i in range(len(sequence_files)):
            if cds:
                protein_seqs = translate_cds(read_fasta(sequence_files[i], prefix=ids[i]))
                protein_sequences.append(protein_seqs)
            else:
                protein_sequences.append(read_fasta(sequence_files[i], prefix=ids[i]))

        logging.info('Writing blastdb sequences to seqs.fasta.')
        db = os.path.join(output_dir,'seqs.fasta')
        write_fasta(protein_sequences[0], db)
        query = db

        if one_v_one:
            query = os.path.join(output_dir, 'query.fasta')
            logging.info('Writing query sequences to query.fasta.')
            write_fasta(protein_sequences[1], query)

        logging.info('Performing all_v_all_blastp (this might take a while)')
        blast_results = all_v_all_blast(query, db, output_dir, eval_cutoff=eval_cutoff)

    if one_v_one:
        logging.info('Retrieving one vs. one orthologs')
        get_one_v_one_orthologs_rbh(blast_results, output_dir)

    if mcl:
        logging.info('Performing MCL clustering (inflation factor = {0})'.format(inflation_factor))
        ava_graph = ava_blast_to_abc_2(blast_results)
        mcl_out = run_mcl_ava_2(ava_graph, output_dir=output_dir, output_file='{}.paranome.mcl'.format(sequences),
                                inflation=inflation_factor)

    logging.info('Done')
    pass


# Ks ANALYSIS USING JOBLIB/ASYNC  --------------------------------------------------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--gene_families', '-gf', default=None,
              help='Gene families (paralogs or one-to-one orthologs). Every '
                   'family should be provided as a tab separated line of gene IDs. See also `wgd blast`.')
@click.option('--sequences', '-s', default=None,
              help='CDS sequences file in fasta format. Multiple files should be specified separated by commas.')
@click.option('--output_directory', '-o', default='ks.out',
              help='Output directory. (Default = ks.out)')
@click.option('--protein_sequences', '-ps', default=None,
              help="Protein sequences fasta file. Optional since by default the CDS file will be translated (universal"
                   "genetic code, with no ambiguous nucleotides). Multiple files should be specified separated by "
                   "commas.")
@click.option('--tmp_dir', '-tmp', default=None,
              help="Path to store temporary files. (Default = automatically generated unique id)")
@click.option('--muscle', '-m', default='muscle',
              help="Absolute path to muscle executable, not necessary if in PATH environment variable.")
@click.option('--codeml', '-c', default='codeml',
              help="Absolute path to codeml executable, not necessary if in PATH environment variable.")
@click.option('--times', '-t', default=1,
              help="Number of times to perform ML estimation (for more stable estimates). (Default = 1)")
@click.option('--min_msa_length', '-mml', default=100,
              help="Minimum MSA length for Ks analysis. (Default = 100)")
@click.option('--ignore_prefixes', is_flag=True,
              help="Ignore gene ID prefixes (defined by the '|' symbol) in the gene families file.")
@click.option('--one_v_one', is_flag=True,
              help="One vs one ortholog distribution.")
@click.option('--preserve', is_flag=True,
              help="Keep multiple sequence alignment and codeml output. ")
@click.option('--async', is_flag=True, default=False,
              help="Use asyncio module for parallelization. (Default uses joblib)")
@click.option('--n_cores','-n', default=4,
              help="Number of CPU cores to use.")
def ks(gene_families, sequences, output_directory, protein_sequences,
        tmp_dir, muscle, codeml, times, min_msa_length, ignore_prefixes, one_v_one, preserve, async, n_cores):
    """
    Construct a Ks distribution.

    Ks distribution construction for a set of paralogs or one-to-one orthologs.
    This implementation uses either the joblib or the asyncio library for parallellization.

    Example 1 - whole paranome Ks distribution:

        wgd ks -gf fringilla_coelebs.mcl -s fringilla_coelebs.cds.fasta -o finch_ks_out --n_cores 8

    Example 2 - one vs. one ortholog Ks distribution:

        wgd ks -gf beaver_eagle -s castor_fiber.cds.fasta,aquila_chrysaetos.cds.fasta -o beaver_eagle_ks_out
    """
    # lazy imports
    from wgd.ks_distribution import ks_analysis_paranome, ks_analysis_one_vs_one
    from wgd.viz import plot_selection, syntenic_dotplot, syntenic_dotplot_ks_colored

    # input check
    if not (gene_families and sequences):
        logging.error('No gene families or no sequences provided.')
        return 1

    if not tmp_dir:
        tmp_dir = os.path.join('.', str(uuid.uuid4()))

    output_directory = os.path.abspath(output_directory)
    tmp_dir = os.path.abspath(tmp_dir)
    gene_families = os.path.abspath(gene_families)

    if os.path.exists(output_directory):
        logging.warning(
            'Output directory already exists, will possibly overwrite')
    else:
        os.mkdir(output_directory)

    if os.path.exists(tmp_dir):
        logging.warning(
            'tmp directory already exists, be sure not to run two analyses simultaneously '
            'in the same tmp directory as this will mess up the results!')
    else:
        os.mkdir(tmp_dir)

    logging.debug('Reading CDS sequences')
    seq_list = [os.path.abspath(x) for x in sequences.strip().split(',')]
    cds_seqs = {}
    for seq_file in seq_list:
        cds_seqs.update(read_fasta(seq_file))

    # translate CDS file(s)
    if not protein_sequences:
        logging.info('Translating CDS file')
        protein_seqs = translate_cds(cds_seqs)

    else:
        logging.debug('Reading protein sequences')
        seq_list = [os.path.abspath(x) for x in protein_sequences.strip().split(',')]
        protein_seqs = {}
        for seq_file in seq_list:
            protein_seqs.update(read_fasta(seq_file))

    # one-vs-one ortholog input
    if one_v_one:
        os.chdir(tmp_dir)
        logging.info('Started one-vs-one ortholog Ks analysis')
        results = ks_analysis_one_vs_one(cds_seqs, protein_seqs, gene_families, tmp_dir, output_directory,
                                         muscle, codeml, async=async, n_cores=n_cores, preserve=preserve, check=False,
                                         times=times, min_length=min_msa_length)
        logging.info('Generating plots')
        plot_selection(results, output_file=os.path.join(output_directory, '{}.ks.png'.format(os.path.basename(
            gene_families))), title=os.path.basename(gene_families))

    # whole paranome ks analysis
    else:
        os.chdir(tmp_dir)  # change directory to the tmp dir, as codeml writes non-unique file names to the working dir
        logging.info('Started whole paranome Ks analysis')
        results = ks_analysis_paranome(cds_seqs, protein_seqs, gene_families, tmp_dir, output_directory,
                                       muscle, codeml, preserve=preserve, check=False, times=times,
                                       ignore_prefixes=ignore_prefixes, async=async, n_cores=n_cores,
                                       min_length=min_msa_length)
        logging.info('Generating plots')
        plot_selection(results, output_file=os.path.join(output_directory, '{}.ks.png'.format(os.path.basename(
            gene_families))), title=os.path.basename(gene_families))

    logging.info('Done')


# CO-LINEARITY ---------------------------------------------------------------------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--gff_file', '-gff', default=None, help='Annotation in gff3 format')
@click.option('--families', '-gf', default=None, help='Gene families as outline per MCL analysis')
@click.option('--output_dir', '-o', default='./coll_out', help='Output directory')
@click.option('--ks_distribution', '-ks', default=None,
              help="Ks distribution for the whole paranome of the species of interest, "
                   "csv file as generated using `wgd ks`.")
@click.option('--keyword', '-kw', default='mRNA',
              help="Keyword for parsing the genes from the GFF file (column 3). (Default = 'mRNA').")
@click.option('--id_string', '-id', default='ID',
              help="Keyword for parsing the gene IDs from the GFF file (column 9). (Default = 'ID').")
def syn(gff_file, families, output_dir, ks_distribution, keyword, id_string):
    """
    Co-linearity analyses.
    Requires I-ADHoRe
    """
    # lazy imports
    from wgd.collinearity import write_families_file, write_gene_lists, write_config_adhore, run_adhore
    from wgd.collinearity import get_anchor_pairs
    from wgd.viz import plot_selection, syntenic_dotplot, syntenic_dotplot_ks_colored

    # input check
    if not gff_file:
        logging.error('No gff file provided! Run `wgd syn --help` for usage instructions.')
        return 1

    if not families:
        logging.error('No gene families provided! Run `wgd syn --help` for usage instructions.')
        return 1

    if os.path.exists(output_dir):
        logging.warning(
            'Output directory already exists, will possibly overwrite')

    else:
        os.mkdir(output_dir)
        logging.info('Made output directory {0}'.format(output_dir))

    logging.info("Parsing GFF file")
    genome = Genome()
    genome.parse_plaza_gff(gff_file, keyword=keyword, id_string=id_string)

    logging.info("Writing gene lists")
    all_genes = write_gene_lists(
        genome, os.path.join(output_dir, 'gene_lists'))

    logging.info("Writing families file")
    write_families_file(families, all_genes,
                        os.path.join(output_dir, 'families.tsv'))

    logging.info("Writing configuration file")
    write_config_adhore(os.path.join(output_dir, 'gene_lists'), os.path.join(output_dir, 'families.tsv'),
                        config_file_name=os.path.join(output_dir, 'adhore.conf'),
                        output_path=os.path.join(output_dir, 'i-adhore-out'))

    logging.info("Running I-ADHoRe")
    run_adhore(os.path.join(output_dir, 'adhore.conf'))

    logging.info('Drawing co-linearity dotplot')
    syntenic_dotplot(pd.read_csv(os.path.join(output_dir, 'i-adhore-out', 'multiplicons.txt'), sep='\t'),
                     output_file=os.path.join(output_dir, 'dotplot.png'))

    if ks_distribution:
        logging.info("Constructing Ks distribution for anchors")
        ks, anchors = get_anchor_pairs(os.path.join(output_dir, 'i-adhore-out', 'anchorpoints.txt'), ks_distribution,
                                       out_file=os.path.join(output_dir, 'ks_anchors.csv'))

        syntenic_dotplot_ks_colored(pd.read_csv(os.path.join(output_dir, 'i-adhore-out', 'multiplicons.txt'), sep='\t'),
                                    pd.read_csv(os.path.join(output_dir, 'i-adhore-out', 'anchorpoints.txt'), sep='\t'),
                                    anchors, output_file=os.path.join(output_dir, 'dotplot.ks.png'))

        logging.info("Generating histogram")
        plot_selection([ks, anchors], alphas=[0.2,0.7], output_file=os.path.join(output_dir, '{}.ks.png'.format(
            os.path.basename(families))), title=os.path.basename(families))

    logging.info("Done")


# MIXTURE MODELING -----------------------------------------------------------------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--ks_distribution', '-ks', default=None, help="Ks distribution csv file, as generated with `wgd ks`.")
@click.option('--method', type=click.Choice(['bgmm', 'gmm', 'both']), default='bgmm',
              help="Mixture modeling method, default is `bgmm` (Bayesian Gaussian mixture model).")
@click.option('--n_range', '-n', default='1,4',
              help='Range of number of components to fit. Default = 1,4')
@click.option('--ks_range', '-r', default='0.1,3',
              help='Ks range to use for modeling. Default = 0.1,3')
@click.option('--output_dir', '-o', default=None, help='Output directory')
@click.option('--gamma', '-g', default=1,
              help='Gamma parameter (inverse of regularization strength) for bgmm models. Default = 1')
@click.option('--sequences', '-s', default=None,
              help='Corresponding sequence files, if provided then the paralogs corresponding to each component '
                   'will be in the output.')
def mix(ks_distribution, method, n_range, ks_range, output_dir, gamma, sequences):
    """
    Mixture modeling of Ks distributions
    """
    # lazy imports
    from wgd.modeling import mixture_model_bgmm, mixture_model_gmm

    if not ks_distribution:
        logging.error('No Ks distribution provided! Run `wgd mix --help` for usage instructions.')
        return 1

    if not output_dir:
        output_dir = './{0}_mixture_{1}'.format(method, datetime.datetime.now().strftime("%d%m%y_%H%M%S"))
        logging.info('Output will be in {}'.format(output_dir))
        os.mkdir(output_dir)

    logging.info("Reading Ks distribution")
    df = pd.read_csv(ks_distribution, index_col=0, sep='\t')
    ks_range = [float(x) for x in ks_range.strip().split(',')]
    n_range = [int(x) for x in n_range.strip().split(',')]

    df = df[df['Ks'] < ks_range[1]]
    df = df[df['Ks'] > ks_range[0]]
    df = df.drop_duplicates(keep='first')
    df = df.dropna()

    if method == 'bgmm' or method == 'both':
        logging.info('Started Bayesian Gaussian Mixture modeling (BGMM)')
        models_bgmm = mixture_model_bgmm(df, n_range=n_range, plot_save=True, output_dir=output_dir,
                                         output_file='bgmm.mixture.png', Ks_range=ks_range, gamma=gamma)

    if method == 'gmm' or method == 'both':
        logging.info('Started Gaussian Mixture modeling (GMM)')
        models_gmm, bic, aic, best = mixture_model_gmm(df, Ks_range=ks_range, n=n_range[1], output_dir=output_dir,
                                                       output_file='gmm.mixture.png')

    # TODO Add both method with comparison plots (3 panels ~ cedalion)
    # TODO Get paralogs method (see lore notebook) and finetune plots


@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--ks_distributions', '-ks', default=None,
              help="Comma-separated Ks distribution csv files, as generated with `wgd ks`.")
@click.option('--alpha_values', '-a', default=None,
              help="Comma-separated alpha values, optional.")
@click.option('--colors', '-c', default=None,
              help="Comma-separated colors, optional.")
@click.option('--labels', '-l', default=None,
              help="Comma-separated labels (for legend), optional.")
@click.option('--hist_type', '-ht', default='barstacked', type=click.Choice(['barstacked', 'step', 'stepfilled']),
              help="Histogram type.")
@click.option('--title', '-t', default='WGD histogram',
              help="Plot title.")
@click.option('--output_file', '-o', default='wgd_hist.png',
              help="Output file, default='wgd_hist.png'.")
def hist(ks_distributions, alpha_values, colors, labels, hist_type, title, output_file):
    """ Plot (stacked) histograms """
    dists = [pd.read_csv(x, sep='\t') for x in ks_distributions.split(',')]
    if alpha_values:
        alpha_values = [float(x) for x in alpha_values.split(',')]
    if colors:
        colors = [x for x in colors.split(',')]
    if not labels:
        labels = [x for x in ks_distributions.split(',')]
    else:
        labels = labels.split(',')

    plot_selection(dists, alphas=alpha_values, colors=colors, labels=labels, output_file=output_file,
                   title=title, histtype=hist_type)
    pass


@cli.command(context_settings={'help_option_names': ['-h', '--help']})
def pipeline_1(sequences, gff_file, output_dir):
    """ Standard workflow whole paranome Ks. """
    # TODO
    pass


@cli.command(context_settings={'help_option_names': ['-h', '--help']})
def pipeline_2(sequences, output_dir):
    """ Standard workflow one-vs-one ortholog Ks."""
    # TODO
    pass


if __name__ == '__main__':
    cli()
