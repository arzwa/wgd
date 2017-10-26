#!/usr/bin/python3.5
"""
Arthur Zwaenepoel
"""
# TODO: separate subcommand for mixture modeling of Ks distributions
# TODO: this can than also include a mixture model + peak based paralog extraction tool?

import click
import coloredlogs
import logging
import sys
import os
import re
import datetime
import pandas as pd
from wgd.modeling import mixture_model_bgmm, mixture_model_gmm
from wgd.ks_distribution_ import ks_analysis_paranome, ks_analysis_one_vs_one, ks_analysis_paranome_2
from wgd.positive_selection import PositiveSelection
from wgd.mcl import run_mcl_ava, all_v_all_blast, run_mcl_ava_2, ava_blast_to_abc_2, family_stats
from wgd.utils import check_dirs, translate_cds, read_fasta, write_fasta, prefix_fasta, prefix_multi_fasta, prefix_mcl
from wgd.collinearity import write_families_file, write_gene_lists, write_config_adhore, run_adhore
from wgd.collinearity import segments_to_chords_table, visualize, get_anchor_pairs, stacked_histogram
from wgd.gff_parser import Genome


# CLI ENTRYPOINT -------------------------------------------------------------------------------------------------------
@click.group()
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
    coloredlogs.install()

    if verbose == 'debug':
        logging.basicConfig(format='%(asctime)s: %(levelname)s\t%(message)s', level=logging.DEBUG,
                            stream=sys.stdout)

    elif verbose == 'info':
        logging.basicConfig(format='%(asctime)s: %(levelname)s\t%(message)s', level=logging.INFO,
                            stream=sys.stdout)
    pass


# BLAST AND MCL --------------------------------------------------------------------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--cds', is_flag=True, help='Sequences are CDS.')
@click.option('--mcl', is_flag=True, help='Perform MCL clustering.')
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
def blast(cds, mcl, sequences, species_ids, blast_results, inflation_factor, eval_cutoff, output_dir):
    """
    Perform all-vs.-all Blastp (+ MCL) analysis.
    """

    if not os.path.exists(output_dir):
        logging.info('Output directory: {} does not exist, will make it.'.format(output_dir))
        os.mkdir(output_dir)

    if not blast_results:
        if not sequences:
            logging.error('No sequences provided (use the -s flag)')

        sequence_files = sequences.strip().split(',')

        if species_ids:
            ids = species_ids.strip().split(',')
            if len(ids) != len(sequence_files):
                logging.error('Number of species identifiers ({0}) does not match number of provided sequence '
                              'files ({1}).'.format(len(ids), len(sequence_files)))
        else:
            ids = [''] * len(sequence_files)

        if cds:
            logging.info("CDS sequences provided, will first translate.")

        sequences_dict = {}
        for i in range(len(sequence_files)):
            if cds:
                protein_seqs = translate_cds(read_fasta(sequence_files[i], prefix=ids[i]))
                sequences_dict.update(protein_seqs)
            else:
                sequences_dict.update(read_fasta(sequence_files[i], prefix=ids[i]))

        logging.info('Writing merged sequences file to seqs.fasta.')
        write_fasta(sequences_dict, os.path.join(output_dir,'seqs.fasta'))

        logging.info('Performing all_v_all_blastp (this might take a while)')
        blast_results = all_v_all_blast(os.path.join(output_dir,'seqs.fasta'), output_dir, eval_cutoff=eval_cutoff)

    if mcl:
        logging.info('Performing MCL clustering (inflation factor = {0})'.format(inflation_factor))
        ava_graph = ava_blast_to_abc_2(blast_results)
        mcl_out = run_mcl_ava_2(ava_graph, output_dir=output_dir, output_file='out.mcl',
                                inflation=inflation_factor)

        # family_stats(mcl_out) # TODO add some statistics from the mcl output

    logging.info('Done')
    pass


# PREFIX ---------------------------------------------------------------------------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('input_file')
@click.option('--fasta', is_flag=True, help='Fasta file input')
@click.option('--mcl', is_flag=True, help='MCL output file input')
@click.option('--prefixes','-p', default='soi',
              help="Prefix(es), multiple prefixes should be provided as a comma separated string "
                   "(don't forget the regexes).")
@click.option('--regexes','-r', default='.+',
              help="Regular expressions, comma separated, in the same order as the prefixes.")
@click.option('--output_file','-o', default=None,
              help='Output file name')
def prefix(input_file, fasta, prefixes, regexes, mcl, output_file):
    """
    Add prefixes to fasta or gene family data
    """
    if not output_file:
        output_file = input_file + '.prefixed'

    multi_fasta = False
    if len(prefixes.split(',')) > 1:
        multi_fasta = True

    prefix_list = prefixes.split(',')
    regex_list = [re.compile(x) for x in regexes.split(',')]
    if len(prefix_list) != len(regex_list):
        logging.error('Number of prefixes is different from number of regular expressions')
    prefix_dict = {prefix_list[i]: regex_list[i] for i in range(len(regex_list))}

    if fasta and not multi_fasta:
        prefix_fasta(prefixes, input_file, output_file)

    if fasta and multi_fasta:
        prefix_multi_fasta(prefix_dict, input_file, output_file)

    if mcl:
        prefix_mcl(prefix_dict, input_file, output_file)


# Ks ANALYSIS USING JOBLIB/ASYNC  --------------------------------------------------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--gene_families', '-gf', default=None,
              help='Gene families (paralogs or one-to-one orthologs). Every '
                   'family should be provided as a tab separated line of gene IDs.')
@click.option('--sequences', '-s', default=None,
              help='CDS sequences file in fasta format.')
@click.option('--species_prefixes', '-p', default=None,
              help='Comma-separated species prefixes for one-to-one ortholog '
                   'distributions (not necessary for whole paranome distributions).')
@click.option('--output_directory', '-o', default='ks.out',
              help='Output directory (should not yet exist). (Default = ks.out)')
@click.option('--protein_sequences', '-ps', default=None,
              help="Protein sequences fasta file. Optional since by default the CDS file will be translated.")
@click.option('--tmp_dir', '-tmp', default='./',
              help="Path to store temporary files. (Default = ./)")
@click.option('--muscle', '-m', default='muscle',
              help="Absolute path to muscle executable, not necessary if in PATH environment variable.")
@click.option('--codeml', '-c', default='codeml',
              help="Absolute path to codeml executable, not necessary if in PATH environment variable.")
@click.option('--times', '-t', default=1,
              help="Number of times to perform ML estimation (for more stable estimates). (Default = 1)")
@click.option('--ignore_prefixes', is_flag=True,
              help="Ignore gene ID prefixes (defined by the '|' symbol) in the gene families file.")
@click.option('--preserve', is_flag=True,
              help="Keep multiple sequence alignment and codeml output. ")
@click.option('--no_prompt', is_flag=True,
              help="Disable prompt for directory clearing.")
@click.option('--async', is_flag=True, default=True,
              help="Use asyncio module for parallelization. (Default)")
@click.option('--joblib', is_flag=True, default=True,
              help="Use joblib module for parallelization, you might want to set the --n_cores parameter. ")
@click.option('--n_cores','-n', default=4,
              help="Number of CPU cores to use.")
def ks(gene_families, sequences, species_prefixes, output_directory, protein_sequences,
        tmp_dir, muscle, codeml, times, ignore_prefixes, preserve, no_prompt, async, joblib, n_cores):
    """
    Construct a Ks distribution.

    Ks distribution construction for a set of paralogs or one-to-one orthologs.
    This implementation uses the joblib library for parallellization.
    """
    if not (gene_families and sequences):
        logging.error('No gene families or no sequences provided.')

    tmp_dir = os.path.join(tmp_dir, output_directory + '.tmp')
    check_dirs(tmp_dir, output_directory, prompt=not no_prompt, preserve=preserve)

    logging.debug("Constructing KsDistribution object")

    # translate CDS file(s)
    if not protein_sequences:
        logging.info('Translating CDS file')
        protein_seqs = translate_cds(read_fasta(sequences))
        protein_sequences = os.path.join(sequences + '.tfa')
        write_fasta(protein_seqs, protein_sequences)

    # one-vs-one ortholog input
    if species_prefixes:
        sp = species_prefixes.strip().split(',')
        if len(sp) != 2:
            logging.error('Number of species prefixes provided differs from 2.')
        s1, s2 = sp[0], sp[1]
        logging.info('Species prefixes provided, will calculate a one-vs-one ortholog Ks distribution for'
                     '{0} vs {1}'.format(s1, s2))
        ks_analysis_one_vs_one(sequences, protein_sequences, gene_families, s1, s2, tmp_dir,
                               output_directory, muscle, codeml, preserve=preserve, check=False, times=times)

    # whole paranome ks analysis
    else:
        os.chdir(tmp_dir)  # change directory to the tmp dir, as codeml writes non-unique file names to the working dir
        if joblib:
            logging.info('Started whole paranome Ks analysis')
            ks_analysis_paranome_2(sequences, protein_sequences, gene_families, tmp_dir, output_directory,
                                 muscle, codeml, preserve=preserve, check=False, times=times,
                                 ignore_prefixes=ignore_prefixes, n_cores=n_cores)
        else:
            logging.info('Started whole paranome Ks analysis')
            ks_analysis_paranome(sequences, protein_sequences, gene_families, tmp_dir, output_directory,
                                 muscle, codeml, preserve=preserve, check=False, times=times,
                                 ignore_prefixes=ignore_prefixes)

    logging.info('Done')


# COLLINEARITY ---------------------------------------------------------------------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('gff_file')
@click.argument('families')
@click.argument('output_dir')
@click.option('--ks_distribution', '-ks', default=None,
              help="Ks distribution for the whole paranome of the species of interest, "
                   "csv file as generated using `wgd ks`.")
@click.option('--keyword', '-kw', default='mRNA',
              help="Keyword for parsing the genes from the GFF file (column 3). (Default = 'mRNA').")
@click.option('--id_string', '-id', default='ID',
              help="Keyword for parsing the gene IDs from the GFF file (column 9). (Default = 'ID').")
def coll(gff_file, families, output_dir, ks_distribution, keyword, id_string):
    """
    Collinearity analyses.
    Requires I-ADHoRe
    """
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

    logging.info("Generating genome.json")
    genome.karyotype_json(out_file=os.path.join(output_dir, 'genome.json'))

    logging.info("Generating visualization")
    segments_to_chords_table(os.path.join(output_dir, 'i-adhore-out', 'segments.txt'),
                             genome, output_file=os.path.join(output_dir, 'chords.tsv'))
    visualize(output_dir)

    if ks_distribution:
        logging.info("Constructing Ks distribution for anchors")
        ks, anchors = get_anchor_pairs(os.path.join(output_dir, 'i-adhore-out', 'anchorpoints.txt'), ks_distribution,
                                   out_file=os.path.join(output_dir, 'ks_anchors.csv'))

        logging.info("Generating histogram")
        stacked_histogram(ks_dist=ks, anchors=anchors,
                          out_file=os.path.join(output_dir, 'histogram.png'))

    logging.info("Done")


# TRANSLATE CDS SEQUENCES ----------------------------------------------------------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--fasta_dir', '-d', default=None, help="Directory with fasta files to translate.")
@click.option('--fasta_file', '-f', default=None, help="Fasta file to translate.")
def trans(fasta_dir, fasta_file):
    """
    Translate a CDS fasta file
    """
    fastas = []

    if fasta_dir:
        fastas = [os.path.join(fasta_dir, x) for x in os.listdir(fasta_dir)]

    elif fasta_file:
        fastas.append(fasta_file)

    else:
        logging.error('Neither fasta file nor directory provided!')

    for ffile in fastas:
        logging.info('Translating {}'.format(ffile))
        protein_seqs = translate_cds(read_fasta(ffile))
        protein_sequences = os.path.join(ffile + '.tfa')
        write_fasta(protein_seqs, protein_sequences)

    logging.info('DONE')


# MIXTURE MODELING -----------------------------------------------------------------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--ks_distribution', '-ks', default=None, help="Ks distribution csv file, as generated with `wgd ks`.")
@click.option('--method', type=click.Choice(['bgmm', 'gmm', 'both']), default='bgmm',
              help="Mixture modeling method, default is `bgmm` (Bayesian Gaussian mixture model).")
@click.option('--n_range', '-n', default='1,4',
              help='Range of number of components to fit. Default = 1,4')
@click.option('--ks_range', '-r', default='0.1,3',
              help='Ks range to use for modeling. Default = 0.1,3')
@click.option('--output_dir', '-o', default='./mixtures_[TIMESTAMP]', help='Output directory')
@click.option('--gamma', '-g', default=1,
              help='Gamma parameter (inverse of regularization strength) for bgmm models. Default = 1')
@click.option('--sequences', '-s', default=None,
              help='Corresponding sequence files, if provided then the paralogs corresponding to each component '
                   'will be in the output.')
def mix(ks_distribution, method, n_range, ks_range, output_dir, gamma, sequences):
    if output_dir == 'mixtures_[TIMESTAMP]':
        output_dir = './mixture' + datetime.datetime.now().strftime("%d%m%y_%H%M%S")
        logging.info('Output will be in {}'.format(output_dir))
        os.mkdir(output_dir)

    logging.info("Reading Ks distribution")
    df = pd.read_csv(ks_distribution, index_col=0)
    ks_range = ks_range.split(',')
    n_range = n_range.split(',')

    df = df[df['Ks'] < ks_range[1]]
    df = df[df['Ks'] > ks_range[0]]
    df = df.drop_duplicates(keep='first')
    df = df.dropna()

    if method == 'bgmm' or method == 'both':
        models_bgmm = mixture_model_bgmm(df, n_range=n_range, plot_save=True, output_dir=output_dir,
                                         output_file='bgmm.mixture.png', Ks_range=ks_range)

    if method == 'gmm' or method == 'both':
        models_gmm = mixture_model_gmm(df, Ks_range=ks_range, n=n_range[1], output_dir=output_dir,
                                       output_file='gmm.mixture.png')

    # TODO Add both method with comparison plots (3 panels ~ cedalion)
    # TODO Get paralogs method (see lore notebook) and finetune plots


if __name__ == '__main__':
    cli()
