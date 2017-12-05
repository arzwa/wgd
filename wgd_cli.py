#!/usr/bin/python3.5
"""
Arthur Zwaenepoel - 2017

The idea behind the wgd CLI is to provide script-like programs bundled under one command (`wgd`)
for common WGD analyses for which the wgd package provides the underlying modular functionalities.
"""

# keep these imports to a minimum tos peed up initial CLI loading
import click
import logging
import sys
import os
import datetime
import pandas as pd
import uuid
import coloredlogs
from wgd.utils import translate_cds, read_fasta, write_fasta, Genome, can_i_run_software


# CLI ENTRY POINT ------------------------------------------------------------------------------------------------------
@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--verbosity','-v', type=click.Choice(['info', 'debug']),
              default='info', help="Verbosity level, default = info.")
def cli(verbosity):
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
    (arzwa@psb.vib-ugent.be)
    """
    coloredlogs.install(fmt='%(asctime)s: %(levelname)s\t%(message)s', level=verbosity.upper(), stream=sys.stdout)
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
@click.option('--output_dir','-o', default='wgd_blast',
              help='Output directory.')
@click.option('--n_threads','-n', default=4,
              help='Number of threads used by blastp.')
def blast(cds, mcl, one_v_one, sequences, species_ids, blast_results, inflation_factor, eval_cutoff, output_dir,
          n_threads):
    """
    All-vs.-all blastp (+ MCL) analysis.

    Requires ``blastp``, ``makeblastdb`` (ncbi-blast+ suite) and ``mcl``.

    Example 1 - whole paranome delineation:

        wgd blast --cds --mcl -s thorny_devil.fasta -o thorny_devil_blast_out

    Example 2 - one vs. one ortholog delineation:

        wgd blast --cds --one_v_one -s equus_ferus.fasta,ursus_arctos.fasta -id horse,bear -e 1e-8 -o bear_horse_out
    """
    # lazy imports
    blast_(cds, mcl, one_v_one, sequences, species_ids, blast_results, inflation_factor, eval_cutoff, output_dir,
           n_threads)


def blast_(cds=True, mcl=True, one_v_one=False, sequences=None, species_ids=None, blast_results=None,
           inflation_factor=2.0, eval_cutoff=1e-10, output_dir='wgd_blast', n_threads=4):
    """
    All vs. all Blast pipeline
    For usage in the ``wgd`` CLI.

    :param cds: boolean, provided sequences are CDS
    :param mcl: boolean, perform MCL clustering
    :param one_v_one: boolean, identify one vs. one orthologs (reciprocal best hits)
    :param sequences: CDS fasta files, if multiple (for one vs. one ortholog identification), then as a comma-separated
        string e.g. ath.fasta,aly.fasta
    :param species_ids: comma-separated species ids, optional for one-vs-one ortholog delineation
    :param blast_results: precomputed blast results (tab separated blast output style)
    :param inflation_factor: inflation factor for MCL clustering
    :param eval_cutoff: e-value cut off for blastp analysis
    :param output_dir: output directory
    :param n_threads: number of CPU threads to use
    :return: output file name
    """
    # lazy imports
    from wgd.blast_mcl import all_v_all_blast, run_mcl_ava_2, ava_blast_to_abc, get_one_v_one_orthologs_rbh

    # software checks
    software = []
    if not blast_results:
        software += ['makeblastdb', 'blastp']
    if mcl:
        software.append('mcl')
    if can_i_run_software(software) == 1:
        logging.error('Could not run all required software, will exit here.')
        return

    # input checks
    if not sequences and not blast_results:
        logging.error('No sequences nor blast results provided! Please use the --help flag for usage instructions.')
        return

    if not os.path.exists(output_dir):
        logging.info('Output directory: {} does not exist, will make it.'.format(output_dir))
        os.mkdir(output_dir)

    # all vs. all blast
    if not blast_results:
        sequence_files = sequences.strip().split(',')

        # input checks
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
                return
        elif one_v_one:
            ids = [os.path.basename(x) for x in sequence_files]
        else:
            ids = [''] * len(sequence_files)

        # get protein sequences
        if cds:
            logging.info("CDS sequences provided, will first translate.")

        protein_sequences = []
        for i in range(len(sequence_files)):
            if cds:
                protein_seqs = translate_cds(read_fasta(sequence_files[i], prefix=ids[i]))
                protein_sequences.append(protein_seqs)
            else:
                protein_sequences.append(read_fasta(sequence_files[i], prefix=ids[i]))

        # blast
        logging.info('Writing blastdb sequences to db.fasta.')
        db = os.path.join(output_dir, str(uuid.uuid4()) + '.db.fasta')
        write_fasta(protein_sequences[0], db)
        query = db

        if one_v_one:
            query = os.path.join(output_dir, str(uuid.uuid4()) + '.query.fasta')
            logging.info('Writing query sequences to query.fasta.')
            write_fasta(protein_sequences[1], query)

        logging.info('Performing all-vs.-all Blastp (this might take a while)')
        blast_results = all_v_all_blast(query, db, output_dir, output_file='{}.blast.tsv'.format(
            '_'.join([os.path.basename(x) for x in sequence_files])), eval_cutoff=eval_cutoff, n_threads=n_threads)
        logging.info('Blast done')

        # remove redundant files
        os.system('rm {}'.format(db))
        if db != query:
            os.system('rm {}'.format(query))

    # get one-vs-one orthologs (RBHs)
    if one_v_one:
        logging.info('Retrieving one vs. one orthologs')
        one_v_one_out = get_one_v_one_orthologs_rbh(blast_results, output_dir)
        logging.info('Done')
        return one_v_one_out

    # get paranome (MCL)
    if mcl:
        logging.info('Performing MCL clustering (inflation factor = {0})'.format(inflation_factor))
        ava_graph = ava_blast_to_abc(blast_results)
        mcl_out = run_mcl_ava_2(ava_graph, output_dir=output_dir, output_file='{}.mcl'.format(
            os.path.basename(blast_results)), inflation=inflation_factor)
        logging.info('Done')
        return mcl_out

    return blast_results


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
@click.option('--aligner', '-a', default='muscle', type=click.Choice(['muscle', 'prank']),
              help='Aligner program to use, from fast to slow: muscle, prank. (Default = muscle)')
@click.option('--times', '-t', default=1,
              help="Number of times to perform ML estimation (for more stable estimates). (Default = 1)")
@click.option('--min_msa_length', '-mml', default=100,
              help="Minimum MSA length for Ks analysis. (Default = 100)")
@click.option('--n_threads','-n', default=4,
              help="Number of threads to use.")
@click.option('--wm', '-w', type=click.Choice(['alc', 'fasttree', 'phyml']), default='fasttree',
              help="Node weighting method, from fast to slow: alc, fasttree, phyml")
@click.option('--ignore_prefixes', is_flag=True,
              help="Ignore gene ID prefixes (defined by the '|' symbol) in the gene families file.")
@click.option('--one_v_one', is_flag=True,
              help="One vs one ortholog distribution.")
@click.option('--preserve', is_flag=True,
              help="Keep multiple sequence alignment and codeml output. ")
@click.option('--async', is_flag=True, default=False,
              help="Use asyncio module for parallelization. (Default uses joblib)")
def ks(gene_families, sequences, output_directory, protein_sequences, tmp_dir, aligner, times, min_msa_length,
       n_threads, wm, ignore_prefixes, one_v_one, preserve, async):
    """
    Ks distribution construction.

    Ks distribution construction for a set of paralogs or one-to-one orthologs.
    This implementation uses either the joblib or the asyncio library for parallellization.
    Requires both ``codeml`` and ``muscle``. Depending on the weighting method chosen (``wm`` option)
    ``phyml`` or ``FastTree`` might also be required.

    Example 1 - whole paranome Ks distribution:

        wgd ks -gf fringilla_coelebs.mcl -s fringilla_coelebs.cds.fasta -o finch_ks_out --n_cores 8

    Example 2 - one vs. one ortholog Ks distribution:

        wgd ks -gf beaver_eagle -s castor_fiber.cds.fasta,aquila_chrysaetos.cds.fasta -o beaver_eagle_ks_out
    """
    ks_(gene_families, sequences, output_directory, protein_sequences, tmp_dir, aligner, muscle='muscle', codeml='codeml',
        times=times, min_msa_length=min_msa_length, ignore_prefixes=ignore_prefixes, one_v_one=one_v_one,
        preserve=preserve, async=async, n_threads=n_threads, weighting_method=wm)


def ks_(gene_families, sequences, output_directory, protein_sequences=None, tmp_dir=None, aligner='muscle',
        muscle='muscle', codeml='codeml', times=1, min_msa_length=100, ignore_prefixes=False, one_v_one=False,
        preserve=False, async=False, n_threads=4, weighting_method='fasttree'):
    """
    Ks distribution construction pipeline. For usage in the ``wgd`` CLI.

    :param gene_families: gene families, i.e. tab separated paralogs or one-vs-one orthologs (see :py:func:`blast_`)
    :param sequences: CDS fasta files, if multiple (one-vs.-one ortholog distribution) then as a comma separated string
    :param output_directory: output directory
    :param protein_sequences: potein sequences (optional)
    :param tmp_dir: tmp directory name (optional)
    :param muscle: path to MUSCLE executable
    :param codeml: path to codeml executable
    :param times: number of times to iteratively perform ML estimation of Ks, Ka and omega values.
    :param min_msa_length: minimum multiple sequence alignment length
    :param ignore_prefixes: ignore prefixes defined by '|' in gene IDs
    :param one_v_one: boolean, one-vs.-one ortholog analysis
    :param preserve: boolean, preserve codeml output files and multiple sequence alignments?
    :param async: use the async library for parallelization
    :param n_threads: number of CPU threads to use
    :param weighting_method: weighting method (fasttree, phyml or alc)
    :return: output file name
    """
    # lazy imports
    from wgd.ks_distribution import ks_analysis_paranome, ks_analysis_one_vs_one
    from wgd.viz import plot_selection

    # software check
    software = [codeml, aligner]
    if weighting_method == 'fasttree':
        software.append('FastTree')
    elif weighting_method == 'phyml':
        software.append('phyml')
    if can_i_run_software(software) == 1:
        logging.error('Could not run all required external software, exit here.')
        return 1

    # input check
    if not (gene_families and sequences):
        logging.error('No gene families or no sequences provided.')
        return 1

    if not tmp_dir:
        tmp_dir = os.path.join('.', 'ks_tmp.' + str(uuid.uuid4()))  # unique tmp directory

    # get absolute paths before changing dir
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

    # define a base name for the output files
    base = '_'.join([os.path.basename(x) for x in seq_list])

    # one-vs-one ortholog input
    if one_v_one:
        os.chdir(tmp_dir)  # chdir is necessary because codeml produces these rub, rst and rst1 files
        logging.info('Started one-vs-one ortholog Ks analysis')
        results = ks_analysis_one_vs_one(cds_seqs, protein_seqs, gene_families, tmp_dir, output_directory,
                                         muscle, codeml, async=async, n_cores=n_threads, preserve=preserve,
                                         times=times, min_length=min_msa_length, aligner=aligner)
        results.round(5).to_csv(os.path.join(output_directory, '{}.ks.tsv'.format(base)), sep='\t')

        logging.info('Generating plots')
        plot_selection(results, output_file=os.path.join(output_directory, '{}.ks.png'.format(base)),
                       title=os.path.basename(gene_families))

        logging.info('Done')
        return os.path.join(output_directory, '{}.ks.tsv'.format(base))

    # whole paranome ks analysis
    else:
        os.chdir(tmp_dir)  # change directory to the tmp dir, as codeml writes non-unique file names to the working dir
        logging.info('Started whole paranome Ks analysis')
        results = ks_analysis_paranome(cds_seqs, protein_seqs, gene_families, tmp_dir, output_directory,
                                       muscle, codeml, preserve=preserve, times=times, aligner=aligner,
                                       ignore_prefixes=ignore_prefixes, async=async, n_cores=n_threads,
                                       min_length=min_msa_length, method=weighting_method)
        results.round(5).to_csv(os.path.join(output_directory, '{}.ks.tsv'.format(base)), sep='\t')

        logging.info('Generating plots')
        plot_selection(results, output_file=os.path.join(output_directory, '{}.ks.png'.format(base)),
                       title=os.path.basename(gene_families))

        logging.info('Done')
        return os.path.join(output_directory, '{}.ks.tsv'.format(base))


# CO-LINEARITY ---------------------------------------------------------------------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--gff_file', '-gff', default=None, help='Annotation in gff3 format')
@click.option('--gene_families', '-gf', default=None, help='Gene families as outline per MCL analysis')
@click.option('--output_dir', '-o', default='./coll_out', help='Output directory')
@click.option('--ks_distribution', '-ks', default=None,
              help="Ks distribution for the whole paranome of the species of interest, "
                   "csv file as generated using `wgd ks`.")
@click.option('--keyword', '-kw', default='mRNA',
              help="Keyword for parsing the genes from the GFF file (column 3). (Default = 'mRNA').")
@click.option('--id_string', '-id', default='Parent',
              help="Keyword for parsing the gene IDs from the GFF file (column 9). (Default = 'ID').")
def syn(gff_file, gene_families, output_dir, ks_distribution, keyword, id_string):
    """
    Co-linearity analyses.

    Requires I-ADHoRe 3.0.

    Example:

        wgd syn -gff ailuropoda.gff -gf ailuropoda.paranome.mcl -ks panda.ks -o panda.anchors_out
    """
    syn_(gff_file, gene_families, output_dir, ks_distribution, keyword, id_string)


def syn_(gff_file, families, output_dir, ks_distribution, keyword='mRNA', id_string='Parent'):
    """
    Co-linearity analysis with I-ADHoRe 3.0
    For usage in the ``wgd`` CLI.

    :param gff_file: GFF annotation file
    :param families: gene families as tab seperated gene IDs, see :py:func:`blast_`
    :param output_dir: output directory
    :param ks_distribution: Ks distribution tsv file, see :py:func:`ks_`
    :param keyword: keyword for entities of interest in the GFF file, e.g. 'CDS' or 'mRNA'
    :param id_string: ID string for the gene ID in the GFF (9th column), e.g. 'ID' or 'Parent'
    :return: nothing at all
    """
    # lazy imports
    from wgd.colinearity import write_families_file, write_gene_lists, write_config_adhore, run_adhore
    from wgd.colinearity import get_anchor_pairs
    from wgd.viz import plot_selection, syntenic_dotplot, syntenic_dotplot_ks_colored

    # software check
    if can_i_run_software(['i-adhore']) == 1:
        logging.error('Could not run all software, exit here.')
        return 1

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

    # parse the gff
    logging.info("Parsing GFF file")
    genome = Genome()
    try:
        genome.parse_plaza_gff(gff_file, keyword=keyword, id_string=id_string)
    except:
        logging.error('Invalid GFF file, be sure to have 9 columns, with your features of interest marked by the ')
        logging.error('keyword {0} (set with the -kw option) in the third column and the gene too which it refers ')
        logging.error('identified by {1}=gene in the ; separated string in the last column. ')
        return

    # generate necessary files for i-adhore
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

    # run i-adhore
    logging.info("Running I-ADHoRe 3.0")
    run_adhore(os.path.join(output_dir, 'adhore.conf'))

    # dotplot
    logging.info('Drawing co-linearity dotplot')
    syntenic_dotplot(
        pd.read_csv(os.path.join(output_dir, 'i-adhore-out', 'multiplicons.txt'), sep='\t'),
        output_file=os.path.join(output_dir, '{}.dotplot.png'.format(os.path.basename(families))))

    # Ks distribution for acnhors and Ks colored dotplot
    if ks_distribution:
        logging.info("Constructing Ks distribution for anchors")
        ks, anchors = get_anchor_pairs(
            os.path.join(output_dir, 'i-adhore-out', 'anchorpoints.txt'),
            pd.read_csv(ks_distribution, index_col=0, sep='\t'),
            out_file=os.path.join(output_dir, '{}.ks_anchors.tsv'.format(os.path.basename(families)))
        )

        syntenic_dotplot_ks_colored(
            pd.read_csv(os.path.join(output_dir, 'i-adhore-out', 'multiplicons.txt'), sep='\t'),
            pd.read_csv(os.path.join(output_dir, 'i-adhore-out', 'anchorpoints.txt'), sep='\t'),
            anchors, output_file=os.path.join(output_dir, '{}.dotplot.ks.png'.format(os.path.basename(families)))
        )

        logging.info("Generating histogram")
        plot_selection([ks, anchors], alphas=[0.2,0.7], output_file=os.path.join(output_dir, '{}.ks_anchors.png'.format(
            os.path.basename(families))), title=os.path.basename(families), labels=['Whole paranome', 'Anchors'])

    logging.info("Done")


# MIXTURE MODELING -----------------------------------------------------------------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--ks_distribution', '-ks', default=None, help="Ks distribution csv file, as generated with `wgd ks`.")
@click.option('--method', type=click.Choice(['bgmm', 'gmm', 'both']), default='gmm',
              help="Mixture modeling method, default is `gmm` (Gaussian mixture model).")
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
@click.option('--cut_off', '-c', default=0.95,
              help='Minimum probability to belong to a particular component for selected sequence files. '
                   'Default = 0.95')
def mix(ks_distribution, method, n_range, ks_range, output_dir, gamma, sequences, cut_off):
    """
    Mixture modeling of Ks distributions.

    Example:

        wgd mix -ks lama.ks.tsv --method gmm -n 2,5
    """
    mix_(ks_distribution, method, n_range, ks_range, output_dir, gamma, sequences, cut_off)


def mix_(ks_distribution, method, n_range, ks_range, output_dir, gamma, sequences, cut_off=0.95):
    """
    Mixture modeling using sklearn

    :param ks_distribution: Ks distribution tsv file
    :param method: mixture modeling method, either 'bgmm' or 'gmm'
    :param n_range: range of number of components to fit
    :param ks_range: range of Ks values to fit model to
    :param output_dir: output directory
    :param gamma: gamma parameter for regulaization in BGMM (lower, stronger regularization)
    :param sequences: fasta file (if provided will write component-wise fasta files)
    :return: nada
    """
    # lazy imports
    from wgd.modeling import mixture_model_bgmm, mixture_model_gmm, get_component_probabilities
    from wgd.utils import get_paralogs_fasta

    if not ks_distribution:
        logging.error('No Ks distribution provided! Run `wgd mix --help` for usage instructions.')
        return 1

    if not output_dir:
        output_dir = './{0}_mixture_{1}'.format(method, datetime.datetime.now().strftime("%d%m%y_%H%M%S"))
        logging.info('Output will be in {}'.format(output_dir))
        os.mkdir(output_dir)

    elif not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    logging.info("Reading Ks distribution")
    df = pd.read_csv(ks_distribution, index_col=0, sep='\t')
    ks_range = [float(x) for x in ks_range.strip().split(',')]
    n_range = [int(x) for x in n_range.strip().split(',')]

    df = df[df['Ks'] < ks_range[1]]
    df = df[df['Ks'] > ks_range[0]]
    df = df.drop_duplicates(keep='first')
    df = df.dropna()

    if method == 'bgmm':
        logging.info('Started Bayesian Gaussian Mixture modeling (BGMM)')
        models_bgmm = mixture_model_bgmm(df, n_range=n_range, plot_save=True, output_dir=output_dir,
                                         output_file='bgmm.mixture.png', Ks_range=ks_range, gamma=gamma)
        best = models_bgmm[-1]
        logging.info('Bayesian Gaussian mixture option, will take the highest number of ')
        logging.info('components as best model under assumption of sufficient regularization.')

    else:
        logging.info('Started Gaussian Mixture modeling (GMM)')
        models_gmm, bic, aic, best = mixture_model_gmm(df, Ks_range=ks_range, n=n_range[1], output_dir=output_dir,
                                                       output_file='gmm.mixture.png')

    logging.info('Saving data frame with probabilities for each component to {}'.format(
        os.path.join(output_dir, '{}.{}.tsv'.format(ks_distribution, method))))
    df = get_component_probabilities(df, best)
    df.to_csv(os.path.join(output_dir, '{}.{}.tsv'.format(ks_distribution, method)), sep='\t')

    logging.info('Saving fasta file for each component')
    if sequences:
        for i in range(best.n_components):
            selected = df[df['p_component{}'.format(i+1)] >= cut_off]
            get_paralogs_fasta(sequences, selected, os.path.join(output_dir, 'component{}.fasta'.format(i+1)))

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
@click.option('--interactive','-i', is_flag=True,
              help="Interactive visualization with bokeh.")
def viz(ks_distributions, alpha_values, colors, labels, hist_type, title, output_file, interactive):
    """
    Plot histograms/densities (interactively).

    Requires a running bokeh server for interactive visualization.
    Run a bokeh serve instance (in the background) with `bokeh serve &`.

    more information about bokeh: https://bokeh.pydata.org/en/latest/index.html
    """
    viz_(ks_distributions, alpha_values, colors, labels, hist_type, title, output_file, interactive)


def viz_(ks_distributions, alpha_values, colors, labels, hist_type, title, output_file, interactive=False):
    """
    Plot (stacked) histograms (interactively)

    :param ks_distributions: a directory with ks distributions (other files are ignored)
        or a comma-separated string of file names
    :param alpha_values: alpha values for the different distributions (in the same order).
        Only relevant for non-interactive visualization.
    :param colors: as in ``alpha_values`` but for colors
    :param labels: as in ``alpha_values`` but for legend labels (by default the file names are used),
        this is also relevant for the interactive bokeh visualization (as opposed to ``alpha_values`` and ``colors``.
    :param hist_type: histogram type (matplotlib), either 'barstacked', 'step' or 'stepfilled'.
    :param title: plot title
    :param output_file: output file name
    :param interactive: render an interactive bokeh plot. This makes some of the above arguments redundant
    :return: nada
    """
    from wgd.viz import plot_selection

    # basic input check
    if not ks_distributions:
        logging.error('No Ks distributions provided, run `wgd viz -h` for usage instructions')
        logging.error('You have to provide one or more computed Ks distributions! See for example `wgd ks -h`')
        return 1

    # directory provided
    if os.path.isdir(ks_distributions):
        ks_distributions = [os.path.join(ks_distributions, x) for x in os.listdir(ks_distributions)]

    # separate distributions provided
    else:
        ks_distributions = ks_distributions.split(',')

    # get distributions
    dists_files = []
    dists = []
    for i in ks_distributions:
        # check if valid distributions
        try:
            d = pd.read_csv(i, sep='\t', index_col=0)
            _ = d[['Ks', 'WeightOutliersExcluded']]
            dists_files.append(i)
            dists.append(d)
        except:
            logging.info('Not a Ks distribution: {}'.format(i))

    # interactive bokeh visualizatoin
    if interactive:
        from wgd.viz import histogram_bokeh
        histogram_bokeh(dists_files, labels)
        return

    # normal matplotlib plots
    # features
    if alpha_values:
        alpha_values = [float(x) for x in alpha_values.split(',')]
        if len(alpha_values) != len(dists):
            logging.error('Please provide as much alpha values as there are distributions')
    if colors:
        colors = [x for x in colors.split(',')]
        if len(colors) != len(dists):
            logging.error('Please provide as much colors as there are distributions')
    if not labels:
        labels = dists_files
        if len(labels) != len(dists):
            logging.error('Please provide as much labels as there are distributions')
    else:
        labels = labels.split(',')

    # make plot
    logging.info('Plotting Ks distributions overlay')
    plot_selection(dists, alphas=alpha_values, colors=colors, labels=labels, output_file=output_file,
                   title=title, histtype=hist_type)
    return


@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('sequences', default=None)
@click.argument('output_dir', default=None)
@click.option('--gff_file', '-gff', default=None, help="GFF file for co-linearity analysis.")
@click.option('--n_threads','-n', default=4,
              help="Number of threads to use.")
def pipeline_1(sequences, output_dir, gff_file, n_threads):
    """
    Standard workflow whole paranome Ks.

    Example:

        wgd pipeline_1 -gff snail.gff -n 16 snail.fasta snail_wgd_out
    """
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # wgd blast
    blast_dir = os.path.join(output_dir, 'wgd_blast')
    mcl_out = blast_(cds=True, mcl=True, sequences=sequences, n_threads=n_threads, output_dir=blast_dir)

    # wgd ks
    ks_dir = os.path.join(output_dir, 'wgd_ks')
    ks_results = ks_(gene_families=mcl_out, sequences=sequences, n_threads=n_threads, output_directory=ks_dir)

    # wgd syn
    if gff_file:
        syn_dir = os.path.join(output_dir, 'wgd_syn')
        syn_(gff_file=gff_file, families=mcl_out, output_dir=syn_dir, ks_distribution=ks_results)

    return


@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('sequences', default=None)
@click.argument('output_dir', default=None)
@click.option('--n_threads','-n', default=4,
              help="Number of threads to use.")
def pipeline_2(sequences, output_dir, n_threads):
    """
    Standard workflow one-vs-one ortholog Ks.

    Example:

        wgd pipeline_2 -n 8 snail.fasta,whale.fasta snail_vs_whale_ks_out
    """
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # wgd blast
    blast_dir = os.path.join(output_dir, 'wgd_blast')
    ovo_out = blast_(cds=True, one_v_one=True, mcl=False, sequences=sequences, n_threads=n_threads,
                     output_dir=blast_dir)

    # wgd ks
    ks_dir = os.path.join(output_dir, 'wgd_ks')
    ks_(gene_families=ovo_out, sequences=sequences, n_threads=n_threads, output_directory=output_dir, one_v_one=True)

    return


if __name__ == '__main__':
    cli()
