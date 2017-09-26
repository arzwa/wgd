#!/usr/bin/python3.5
"""
Arthur Zwaenepoel
"""

import click
import coloredlogs
import logging
import sys
import os
import re
from wgd.ks_distribution import KsDistribution
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
    if verbose == 'debug':
        coloredlogs.install()
        logging.basicConfig(format='%(asctime)s: %(levelname)s\t%(message)s', level=logging.DEBUG,
                            stream=sys.stdout)

    elif verbose == 'info':
        coloredlogs.install()
        logging.basicConfig(format='%(asctime)s: %(levelname)s\t%(message)s', level=logging.INFO,
                            stream=sys.stdout)

    # click.echo(click.style('WGD v1.0', fg='green'))
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
@click.option('--inflation_factor', '-I', default=2.0,
              help="Inflation factor for MCL clustering, when blast results provided. (Default = 2)")
@click.option('--output_dir','-o', default='wgd.blast.out',
              help='Output directory.')
def blast(cds, mcl, sequences, species_ids, inflation_factor, output_dir):
    """
    Perform all-vs.-all Blastp (+ MCL) analysis.
    """
    if not sequences:
        logging.error('No sequences provided (use the -s flag)')

    logging.info('Output directory: {} does not exist, will make it.'.format(output_dir))
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

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
    all_vs_all = all_v_all_blast(os.path.join(output_dir,'seqs.fasta'), output_dir)

    if mcl:
        logging.info('Performing MCL clustering (inflation factor = {0})'.format(inflation_factor))
        ava_graph = ava_blast_to_abc_2(all_vs_all)
        mcl_out = run_mcl_ava_2(ava_graph, output_dir=output_dir, output_file='out.mcl')
        family_stats(mcl_out)

    pass


# PHYLOGENETIC PROFILE -------------------------------------------------------------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('input_file')
@click.option('--output_file','-o', default=None,
              help='Output file name, if none given than only summary statistics are shown.')
def profile(input_file, output_file):
    """
    Generate a phylogenetic profile.

    The input file should be a typical (Ortho)MCL output file,
    with gene IDs prefixed with their respective species ID.
    """
    family_stats(input_file, pro_file=output_file)


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
    Generate a phylogenetic profile.

    The input file should be a typical (Ortho)MCL output file,
    with gene IDs prefixed with their respective species ID.
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


# Ks ANALYSIS ----------------------------------------------------------------------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('cds_fasta_file', type=click.Path(exists=True))
@click.argument('output_directory', type=click.Path(exists=False))
@click.option('--all_vs_all', '-ava', default=None,
              help="File with all-versus-all Blast results. Will be used for MCL clustering. "
              "Please format your file as `gene 1 [TAB] gene 2 [TAB] e-value`. "
              "Note that you can also provide precomputed gene families using the `--gene_families` flag.")
@click.option('--gene_families', '-gf', default=None,
              help="File with precomputed gene families, e.g. from MCL clustering of all-vs-all Blast results,"
              "OrthoMCL, INPARANOID, OrthoFinder etc. Please format your file with one tab separated gene "
              "family per line. Note that you can also provide raw all-vs-all Blast results using the "
              "`--all_vs_all` flag.")
@click.option('--prefix', '-p', default=None,
              help="Prefix for parsing out the genes of interest, for example 'orysa' for genes named like "
              "`orysa|OS1G11100`. Note that you can also provide a regex pattern (`--regex`). "
              "Note that if neither a prefix nor a regex pattern is provided, all genes are assumed "
              "to be relevant. IMPORTANT, setting a prefix assumes that the sequence IDs in the provided "
              "fasta files are NOT prefixed. (optional)")
@click.option('--regex', '-r', default=None,
              help="Regular expression pattern (Python/Perl type) for parsing out genes of interest. "
              "Especially useful when multi-species gene families are provided. Note that if neither a prefix "
              "nor a regex pattern is provided, all genes are assumed to be relevant. IMPORTANT, setting a "
              "regex assumes that both the genes in the all-vs-all/gene families and sequence input are to be "
              "matched with the regex. (optional)")
@click.option('--inflation_factor', '-I', default=2.0,
              help="Inflation factor for MCL clustering, when blast results provided. (Default = 2)")
@click.option('--protein_sequences', '-ps', default=None,
              help="Protein sequences fasta file. Optional since by default the CDS file will be translated.")
@click.option('--tmp_dir', '-tmp', default='./',
              help="Path to store temporary files. (Default = ./)")
@click.option('--muscle', '-m', default='muscle',
              help="Path to muscle executable, not necessary if in PATH environment variable.")
@click.option('--codeml', '-c', default='codeml',
              help="Path to codeml executable, not necessary if in PATH environment variable.")
@click.option('--times', '-t', default=1,
              help="Number of times to perform ML estimation (for more stable estimates). (Default = 1)")
@click.option('--preserve/--no_preserve', default=False,
              help="Do you want to keep the multiple sequence alignment and codeml output? (Default = False) ")
@click.option('--prompt/--no_prompt', default=True,
              help="Prompt for directory clearing? (Default = True) ")
@click.option('--mixture', '-mix', type=click.Choice(['no', 'bayesian', 'gaussian']), default='bayesian',
              help="Mixture modeling method (requires sklearn.mixture). Set to no if not desired.")
@click.option('--kde/--no_kde', default=True,
              help="Perform mixture modeling (requires sklearn.neighbors.KernelDensity)? (Default = True) ")
def ks(cds_fasta_file, output_directory, all_vs_all, gene_families, prefix, regex, inflation_factor,
       protein_sequences, tmp_dir, muscle, codeml, times, preserve, prompt, mixture, kde):
    """
    Construct a Ks distribution.
    Enables running full analysis pipeline for a species of interest.\n

    \b
    ARGUMENTS (MINIMAL INPUT):
        - CDS_FASTA_FILE:
            File with CDS sequences in fasta format.
        - OUTPUT_DIRECTORY:
            Output directory name, species name is recommended (for plot titles etc.)

    \b
    OUTPUT:
        - csv files and histograms with Ks, Ka and w analysis output
        - modeling plots (Bayesian Gaussian mixture/kde)
        - html report
    """
    # define important directories
    tmp_dir = os.path.join(tmp_dir, output_directory + '.tmp')
    check_dirs(tmp_dir, output_directory, prompt=prompt, preserve=preserve)

    # Basic input checks
    # No families: perform all-v-all first
    if not (gene_families or all_vs_all):
        logging.info("Neither gene families nor all-vs-all Blast results provided, will perform"
                     "all vs. all Blastp first (output in {}/{}.out.blast)".format(output_directory, cds_fasta_file))

        logging.info('Translating CDS file')
        if not protein_sequences:
            protein_seqs = translate_cds(read_fasta(cds_fasta_file))
            protein_sequences = os.path.join(cds_fasta_file + '.tfa')

            with open(os.path.join(cds_fasta_file + '.tfa'), 'w') as o:
                for key, val in protein_seqs.items():
                    o.write('>' + key + '\n')
                    o.write(val + '\n')

        all_vs_all = all_v_all_blast(protein_sequences, output_directory)

    # gene families input
    if gene_families and not (prefix or regex):
        logging.warning("Gene families ({}) provided but no gene pattern, assuming single species gene "
                        "families".format(gene_families))

    if gene_families and (prefix and not regex):
        regex = prefix + '.+'

    logging.info("STARTED Ks ANALYSIS PIPELINE")

    # if blast results provided, perform MCL
    if all_vs_all:
        if os.path.exists(all_vs_all):
            logging.debug("Found all-vs-all Blast results")
        else:
            raise FileNotFoundError(
                "Unable to find blast results file {}".format(all_vs_all))

        gene_families = run_mcl_ava(all_vs_all, regex=regex, prefix=prefix, tmp_dir=tmp_dir,
                                    output_file=os.path.join(
                                        output_directory, 'gene_families.mcl'),
                                    inflation=inflation_factor)

    # perform Ks analysis
    logging.debug("Constructing KsDistribution object")
    ks_dist = KsDistribution(species=output_directory, gene_pattern=regex, gene_families=gene_families,
                             nucleotide=cds_fasta_file, protein=protein_sequences,
                             tmp_dir=tmp_dir, muscle_path=muscle, codeml_path=codeml, output_dir=output_directory)

    logging.info('Started Ks analysis')
    ks_dist.ks_analysis_parallel(preserve=preserve, check=False, times=times)

    if mixture != 'no':
        logging.info('Started {} mixture modeling'.format(mixture))
        ks_dist = ks_dist.mixture_modeling(method=mixture)

    if kde:
        logging.info("Started kernel density estimation")
        ks_dist = ks_dist.kernel_density()

    logging.info("Writing html report")
    ks_dist.write_report(mixture=mixture, kde=kde)


# COLLINEARITY ---------------------------------------------------------------------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('gff_file')
@click.argument('families')
@click.argument('output_dir')
@click.option('--ks_distribution', '-ks', default=None,
              help="Ks distribution for the whole paranome of the species of interest, "
                   "csv file as generated using `wgd ks`.")
def coll(gff_file, families, output_dir, ks_distribution):
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
    genome.parse_plaza_gff(gff_file)

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


# POSITIVE SELECTION SCREENING -----------------------------------------------------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('gene_pattern')
@click.argument('gene_families')
@click.argument('nucleotide_sequences')
@click.argument('output_dir')
@click.option('--job_id', '-i', default='pos',
              help="ID for directory names etc. (e.g. species name). Optional (default = 'pos')")
@click.option('--gene_descriptions', '-gd', default=None,
              help="Gene descriptions to include in data frame. Optional.")
@click.option('--protein_sequences', '-ps', default=None,
              help="Protein sequences fasta file. Optional (will by default translate the CDS file)")
@click.option('--tmp_dir', '-tmp', default='./',
              help="Path to store temporary files. (Default = ./)")
@click.option('--muscle', '-m', default='muscle',
              help="Path to muscle executable, not necessary if in PATH environment variable.")
@click.option('--codeml', '-c', default='codeml',
              help="Path to codeml executable, not necessary if in PATH environment variable.")
@click.option('--times', '-t', default=1,
              help="Number of times to perform ML estimation (for more stable estimates). (Default = 1)")
@click.option('--preserve/--no_preserve', default=True,
              help="Do you want to keep the multiple sequence alignment and codeml output? Default = TRUE. ")
@click.option('--prompt/--no_prompt', default=True,
              help="Prompt for directory clearing? Default = TRUE ")
def pos(gene_pattern, gene_families, protein_sequences, nucleotide_sequences, output_dir,
        job_id, gene_descriptions, tmp_dir, muscle, codeml, times, preserve, prompt):
    """
    Exploratory positive selection analysis.
    Enables running full analysis pipeline for a species of interest.\n

    ARGUMENTS:\n
        - gene pattern: regex pattern for parsing out genes for the species of interest (e.g. ath.+ for genes like
        ath|AT1G00200 or AT.+ for gene IDs like AT1G00200).\n
        - gene families: gene families (e.g. from OrthoMCL), one tab separated family of genes per line\n
        - protein_sequences: fasta file with protein sequences\n
        - nucleotide_sequences: fasta file with corresponding nucleotide sequences\n
        - output_dir: output directory\n

    OUTPUT:\n
        - csv file
    """
    # define important directories
    logging.debug("Defining directory names")
    output_dir = os.path.join(output_dir, '_'.join(job_id.split()))
    tmp_dir = os.path.join(tmp_dir, '_'.join(job_id.split()) + '.tmp')

    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)

    check_dirs(None, output_dir, prompt=prompt, preserve=preserve)

    logging.info('STARTED POSITIVE SELECTION ANALYSIS')
    pos = PositiveSelection(gene_pattern=gene_pattern, gene_families=gene_families, nucleotide=nucleotide_sequences,
                            protein=protein_sequences, species=job_id, output_dir=output_dir, tmp_dir=tmp_dir,
                            muscle_path=muscle, codeml_path=codeml)

    results = pos.positive_selection(
        check=False, preserve=preserve, prompt=True, times=times)

    logging.info('Writing results to csv')
    results.to_csv(os.path.join(output_dir, '{}_out.csv'.format(job_id)))


if __name__ == '__main__':
    cli()
