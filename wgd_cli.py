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

The **command line interface** (CLI) is the main way to interact with the `wgd`
package. The CLI is organized with a Click command that wraps a function with
the same name followed by an underscore (this is chosen mostly so that the
pipeline commands can reuse code from other subcommands).

Upon successful installation you should be able to run::

    $ wgd -h

    Usage: wgd [OPTIONS] COMMAND [ARGS]...

      Welcome to the wgd command line interface!

                             _______
                             \\  ___ `'.
             _     _ .--./)   ' |--.\\  \\
       /\\    \\   ///.''\\    | |    \\  '
       `\\  //\\\\ //| |  | |   | |     |  '
         \\`//  \\'/  \`-' /    | |     |  |
          \\|   |/   /("'`     | |     ' .'
           '        \\ '---.   | |___.' /'
                     /'""'.\ /_______.'/
                    ||     ||\\_______|/
                    \\'. __//
                     `'---'

      wgd  Copyright (C) 2018 Arthur Zwaenepoel
      This program comes with ABSOLUTELY NO WARRANTY;
      This is free software, and you are welcome to redistribute it
      under certain conditions;

      Contact: arzwa@psb.vib-ugent.be

    Options:
      -v, --verbosity [info|debug]  Verbosity level, default = info.
      -l, --logfile TEXT            File to write logs to (optional)
      -h, --help                    Show this message and exit.

    Commands:
      kde  Fit a KDE to a Ks distribution.
      ksd  Ks distribution construction.
      mcl  All-vs.-all blastp + MCL analysis.
      mix  Mixture modeling of Ks distributions.
      syn  Co-linearity analyses.
      viz  Plot histograms/densities (interactively).
      wf1  Standard workflow whole paranome Ks.
      wf2  Standard workflow one-vs-one ortholog Ks.


Note that the verbose flag can be set before the subcommands::

    $ wgd --verbose [silent|info|debug] [COMMAND]

Which sets the verbosity for the logging.

All commands are equipped with usage instructions and documentation of options
and arguments, which can be viewed by using the ``--help`` or ``-h`` flag. These
should be quite self-explanatory, but for further documentation you can refer to
the documentation of the specific functions that are called. These can be found
on this page (e.g. the function called by ``wgd blast`` is
:py:func:`wgd_cli.blast_`).

For more information on the methods used in ``wgd`` to compute |Ks| distributions,
please refer to :ref:`methods`.

Example
=======

Here is a small example on how to use the package through the CLI. This is a
workflow for constructing a |Ks| distribution for a fasta file with CDS
sequences called ``ath.cds.fasta``. File names may be different, but the point
will be clear.

(1) Get the paranome, i.e. perform all-against-all Blastp and MCL clustering,
notice how we specify to use 8 threads::

    $ wgd mcl --cds --mcl -s ath.cds.fasta -o ./ -n 8

(2) Construct a |Ks| distribution, use FastTree for inferring the phylogenetic
trees used in the node weighting procedure::

    $ wgd ksd -o ./ -n 8 ./ath.mcl ath.cds.fasta

(3) Run I-ADHoRe and get an anchor-point |Ks| distribution, as well as dotplots.
here we need a structural annotation in GFF format (see e.g.
https://bioinformatics.psb.ugent.be/plaza/versions/plaza_v4_dicots/download/ for
examples of such files)::

    $ wgd syn ath.gff ath.mcl -ks ath.mcl.ks.tsv -f gene -a ID

(4) Fit some gaussian mixture models with 1 to 5 components (default
parameters)::

    $ wgd mix ath.mcl.ks.tsv -n 1 5

for more information on mixture modeling and some cautionary notes refer to
:ref:`note_on_gmms`

(5) Explore the full and anchor distribution with kernel density estimates
interactively. First run a bokeh server in the background::

    $ bokeh serve &

next, execute the following command::

    $ wgd viz -i -ks ath.mcl.ks.tsv,ath.mcl.ks_anchors.tsv -l full,anchors

a tab in your default browser should appear. See :ref:`viz_info` for more
information on vizualization with ``wgd viz``

--------------------------------------------------------------------------------

Reference
=========
"""
# TODO's & IDEAS:
# - Use biopython for all sequence handling
# - Use tmp files to keep track of the outputs from subprocesses, because using
#   PIPE might result in hangs.
# - Convert input CDS to uppercase (with warning), since wgd blast breaks with
#   lowercase nucleotides

# keep these imports to a minimum to speed up initial CLI loading
import click
import logging
import sys
import os
import warnings
import pandas as pd
import coloredlogs
import subprocess
from wgd.utils import translate_cds, read_fasta, write_fasta, can_i_run_software

warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")


# CLI ENTRY POINT --------------------------------------------------------------
@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.option(
        '--verbosity', '-v', type=click.Choice(['info', 'debug']),
        default='info', help="Verbosity level, default = info."
)
@click.option(
        '--logfile', '-l', default=None,
        help="File to write logs to (optional)"
)
@click.option('--version', is_flag=True, help="Print version number")
def cli(verbosity, logfile, version):
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
    wgd  Copyright (C) 2018 Arthur Zwaenepoel
    This program comes with ABSOLUTELY NO WARRANTY;
    This is free software, and you are welcome to redistribute it
    under certain conditions;

    Contact: arzwa@psb.vib-ugent.be
    """
    if not logfile:
        coloredlogs.install(
                fmt='%(asctime)s: %(levelname)s\t%(message)s',
                level=verbosity.upper(), stream=sys.stdout
        )
    else:
        print('Logs will be written to {}'.format(logfile))
        logging.basicConfig(
                filename=logfile,
                filemode='a',
                format='%(asctime)s: %(levelname)s\t%(message)s',
                datefmt='%H:%M:%S',
                level=verbosity.upper()
        )

    # get round problem with python multiprocessing library that can set all cpu
    # affinities to a single cpu found in OrthoFinder source code
    if sys.platform.startswith("linux"):
        with open(os.devnull, "w") as f:
            subprocess.call("taskset -p 0xffffffffffff %d" % os.getpid(),
                            shell=True, stdout=f)

    if version:
        logging.info("This is wgd v1.1")
    pass


# Check and optionally refactor data -------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('sequences', nargs=-1, type=click.Path(exists=True))
@click.option('--rename', is_flag=True, help='rename sequences')
@click.option('--prefix', default=None, help='prefix strings')
@click.option('--out', default=None, help='output file name prefixes')
def pre(sequences, rename, prefix, out):
    """
    Check and optionally rename CDS files.

    Example usage (renaming)

        wgd pre ath.cds.fasta vvi.cds.fasta --rename --prefix ath,vvi

    Example usage (just checking files)

        wgd pre ath.cds.fasta
    """
    from wgd.pre import check_cds
    if prefix:
        prefixes = prefix.split(",")
    else:
        prefixes = sequences
    if out:
        out = out.split(",")
    else:
        out = [s + ".pre" for s in sequences]
    for i, seq in enumerate(sequences):
        logging.info("({}) Checking {}".format(i, seq))
        check_cds(seq, out[i] + ".good", out[i] + ".bad",
                  rename=rename, prefix=prefixes[i])

# BLAST AND MCL ----------------------------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option(
        '--cds', is_flag=True,
        help='sequences are CDS'
)
@click.option(
        '--mcl', is_flag=True,
        help='perform MCL clustering'
)
@click.option(
        '--one_v_one', is_flag=True,
        help='get one vs. one orthologs'
)
@click.option(
        '--sequences', '-s', default=None,
        help='input fasta files, as a comma separated string (e.g. '
             'x.fasta,y.fasta,z.fasta) or directory name.'
)
@click.option(
        '--species_ids', '-id', default=None,
        help='species identifiers for respective input sequence files, as comma'
             ' separated string (e.g. x,y,z). (optional)'
)
@click.option(
        '--blast_results', '-b', default=None,
        help='input precomputed tab separated blast results'
)
@click.option(
        '--inflation_factor', '-I', default=2.0, show_default=True,
        help="inflation factor for MCL clustering"
)
@click.option(
        '--eval_cutoff', '-e', default=1e-10, show_default=True,
        help="e-value cut-off for blast results"
)
@click.option(
        '--output_dir', '-o', default='wgd_blast', show_default=True,
        help='output directory'
)
@click.option(
        '--n_threads', '-n', default=4, show_default=True,
        help='number of threads used by blastp'
)
def mcl(
        cds, mcl, one_v_one, sequences, species_ids, blast_results,
        inflation_factor, eval_cutoff, output_dir, n_threads
):
    """
    All-vs.-all blastp + MCL clustering.

    Requires blastp, makeblastdb (ncbi-blast+ suite) and mcl. Note the two key
    parameters, being the e-value cut-off and inflation factor. It is advised to
    explore the effects of these on your analysis.

    Example 1 - whole paranome delineation:

        wgd mcl --cds --mcl -s thorny_devil.fasta -o thorny_devil_blast_out

    Example 2 - one vs. one ortholog delineation:

        wgd mcl --cds --one_v_one -s equus_ferus.fasta,ursus_arctos.fasta
        -id horse,bear -e 1e-8 -o bear_horse_out

    wgd  Copyright (C) 2018 Arthur Zwaenepoel
    This program comes with ABSOLUTELY NO WARRANTY;
    """
    # lazy imports
    blast_mcl(cds, mcl, one_v_one, sequences, species_ids, blast_results,
              inflation_factor, eval_cutoff, output_dir, n_threads)


def blast_mcl(
        cds=True, mcl=True, one_v_one=False, sequences=None,
        species_ids=None, blast_results=None, inflation_factor=2.0,
        eval_cutoff=1e-10, output_dir='wgd_blast', n_threads=4
):
    """
    All vs. all Blast + MCL pipeline. For usage in the ``wgd`` CLI. Can be used
    to perform all vs. all Blast, MCL clustering and one vs. one ortholog
    delineation.

    :param cds: boolean, indicates that the provided sequences are CDS
    :param mcl: boolean, perform MCL clustering
    :param one_v_one: boolean, identify whether one vs. one orthologs are to be
        inferred (reciprocal best hits) (True) or a paranome (False).
    :param sequences: CDS fasta files, if multiple (for one vs. one ortholog
        identification), then as a comma-separated string e.g.
        ``ath.fasta,aly.fasta``
    :param species_ids: comma-separated species ids, optional for one-vs-one
        ortholog delineation (will prefix the sequence IDs in that case).
    :param blast_results: precomputed blast results (tab separated blast output
        style)
    :param inflation_factor: inflation factor for MCL clustering
    :param eval_cutoff: e-value cut off for blastp analysis
    :param output_dir: output directory
    :param n_threads: number of threads to use
    :return: output file name
    """
    # lazy imports
    from wgd.blast_mcl import all_v_all_blast, run_mcl_ava, ava_blast_to_abc, \
        get_one_v_one_orthologs_rbh
    from wgd.utils import uniq_id

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
        logging.error('No sequences nor blast results provided! Please use the '
                      '--help flag for usage instructions.')
        return

    if not os.path.exists(output_dir):
        logging.info('Output directory: {} does not exist, will make it.'
                     ''.format(output_dir))
        os.mkdir(output_dir)

    # all vs. all blast
    if not blast_results:
        # get sequences from directory
        if type(sequences) == tuple:
            sequence_files = sequences
        elif os.path.isdir(sequences):
            sequence_files = sorted([os.path.join(sequences, x) for x
                                     in os.listdir(sequences) if
                                     x.endswith(('fa', 'fasta', 'tfa', 'faa'))])
            logging.info('Read the following fasta files from directory {}:'
                         ''.format(sequences))
            logging.info('{}'.format(', '.join(sequence_files)))
        # sequences as comma-separated string
        else:
            sequence_files = sequences.strip().split(',')

        # input checks
        if len(sequence_files) != 1 and not one_v_one:
            logging.info('More then one fasta file provided and one_v_one flag '
                         'not set, will perform all-vs.-all blast.')

        if len(sequence_files) != 2 and one_v_one:
            logging.error('Please two fasta files for one-vs-one ortholog '
                          'finding')
            return

        if species_ids:
            ids = species_ids.strip().split(',')
            if len(ids) != len(sequence_files):
                logging.error('Number of species identifiers ({0}) does not '
                              'match number of provided sequence files ({1}).'
                              ''.format(len(ids), len(sequence_files)))
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
                protein_seqs = translate_cds(read_fasta(
                        sequence_files[i], prefix=ids[i]))
                protein_sequences.append(protein_seqs)
            else:
                protein_sequences.append(read_fasta(
                        sequence_files[i], prefix=ids[i]))

        # blast
        logging.info('Writing blastdb sequences to db.fasta.')
        db = os.path.join(output_dir, uniq_id() + '.db.fasta')

        # all-vs.-all
        d = {}
        for x in protein_sequences:
            d.update(x)
        write_fasta(d, db)
        query = os.path.join(output_dir, uniq_id() + '.query.fasta')
        logging.info('Writing query sequences to query.fasta.')
        write_fasta(d, query)

        # start the blast
        logging.info('Performing all-vs.-all Blastp (this might take a while)')
        blast_results = all_v_all_blast(
                query, db, output_dir,
                output_file='{}.blast.tsv'.format(
                        '_'.join([os.path.basename(x) for x in sequence_files])
                ),
                eval_cutoff=eval_cutoff,
                n_threads=n_threads
        )
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
        logging.info('Performing MCL clustering (inflation factor = {0})'
                     ''.format(inflation_factor))
        ava_graph = ava_blast_to_abc(blast_results)
        mcl_out = run_mcl_ava(
                ava_graph, output_dir=output_dir,
                output_file='{}.mcl'.format(os.path.basename(blast_results)),
                inflation=inflation_factor
        )
        logging.info('Done')
        return mcl_out

    return blast_results


# Diamond + MCL
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('sequences', nargs=-1, type=click.Path(exists=True))
@click.option('--outdir', '-o', default='wgd_dmd', show_default=True,
    help='output directory')
@click.option('--tmpdir', '-t', default=None, show_default=True,
    help='tmp directory')
@click.option('--inflation', '-I', default=2.0, help="inflation factor for MCL")
@click.option('--eval', '-e', default=1e-10, help="e-value cut-off for similarity")
@click.option('--ignorestop', is_flag=True, help="translate through STOP codons")
@click.option('--nostrictcds', is_flag=True, help="do not enforce proper CDS sequences")
def dmd(sequences, outdir, tmpdir, inflation, eval, ignorestop, nostrictcds):
    """
    All-vs.-all diamond blastp + MCL clustering.

    Requires diamond and mcl. Note the two key  parameters, being the e-value
    cut-off and inflation factor. It is advised to explore the effects of these
    on your analysis.

    Example 1 - whole paranome delineation:

        wgd dmd ath.fasta

    Example 2 - one vs. one ortholog delineation:

        wgd dmd ath.fasta vvi.fasta

    Example 3 - one vs. one ortholog delineation for multiple pairs:

        wgd dmd ath.fasta vvi.fasta egr.fasta

    wgd  Copyright (C) 2019 Arthur Zwaenepoel
    This program comes with ABSOLUTELY NO WARRANTY;
    """
    from wgd.diamond import SequenceData
    s = [SequenceData(s, out_path=outdir, tmp_path=tmpdir,
        to_stop=not ignorestop, cds=not nostrictcds) for s in sequences]
    if len(s) == 0:
        logging.error("No sequences provided!")
        return
    elif len(s) == 1:
        logging.info("One CDS file: will compute paranome")
        s[0].get_paranome(inflation=inflation, eval=eval)
        s[0].write_paranome()
    else:
        logging.info("Multiple CDS files: will compute RBH orthologs")
        for i in range(len(s)-1):
            for j in range(i+1, len(s)):
                logging.info("{} vs. {}".format(s[i].prefix, s[j].prefix))
                s[i].get_rbh_orthologs(s[j], eval=eval)
                s[i].write_rbh_orthologs(s[j])
    if tmpdir is None:
        [x.remove_tmp(prompt=False) for x in s]


# Ks ANALYSIS USING JOBLIB/ASYNC  ----------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('gene_families', nargs=1, type=click.Path(exists=True))
@click.argument('sequences', nargs=-1, type=click.Path(exists=True))
@click.option(
        '--output_directory', '-o', default='wgd_ksd', show_default=True,
        help='output directory'
)
@click.option(
        '--protein_sequences', '-p', default=None,
        help="protein sequences fasta file (optional)"
)
@click.option(
        '--tmp_dir', '-tmp', default=None, show_default=True,
        help="path to store temporary files"
)
@click.option(
        '--aligner', '-a', default='mafft', show_default=True,
        type=click.Choice(['muscle', 'prank', 'mafft']),
        help='aligner program to use, from fast to slow: muscle, prank'
)
@click.option(
        '--times', '-t', default=1, show_default=True,
        help="number of times to perform ML estimation"
)
@click.option(
        '--min_msa_length', '-mml', default=100, show_default=True,
        help="minimum MSA length for pairwise mode"
)
@click.option(
        '--n_threads', '-n', default=4, show_default=True,
        help="number of threads to use"
)
@click.option(
        '--wm', '-w', type=click.Choice(['alc', 'fasttree', 'phyml']),
        default='fasttree', show_default=True,
        help="node weighting method, from fast to slow: alc, fasttree, phyml"
)
@click.option(
        '--pairwise', is_flag=True,
        help="perform the analysis using a pairwise approach"
)
@click.option(
        '--max_pairwise', '-mp', default=10000, show_default=True,
        help="maximum number of pairwise combinations a family may have"
)
@click.option(
        '--ignore_prefixes', is_flag=True,
        help="ignore gene ID prefixes (defined by the '|' symbol) in the gene "
             "families file."
)
@click.option(
        '--one_v_one', is_flag=True,
        help="one vs one ortholog distribution"
)
@click.option(
        '--preserve', is_flag=True,
        help="keep multiple sequence alignment, codeml output and trees"
)
@click.option(
        '--codeml_options', default="",
        help="other options for codeml, as a comma separated string, "
             "e.g. getSE=1,CodonFreq=1"
)
def ksd(
        gene_families, sequences, output_directory, protein_sequences, tmp_dir,
        aligner, times, min_msa_length, n_threads, wm, pairwise,
        max_pairwise, ignore_prefixes, one_v_one, preserve, codeml_options
):
    """
    Ks distribution construction.

    Ks distribution construction for a set of paralogs or one-to-one orthologs.
    This implementation uses the joblib library for parallelization.
    Requires both ``codeml`` and ``muscle|mafft|prank`` for alignment. Depending
    on the weighting method  chosen (``wm`` option) ``phyml`` or ``FastTree``
    might also be required.

    Example 1 - whole paranome Ks distribution:

        wgd ksd --n_threads 8 fringilla_coelebs.mcl fringilla_coelebs.cds.fasta

    Example 2 - one vs. one ortholog Ks distribution:

        wgd ksd -o beaver_eagle orthologs.tsv beaver.cds.fasta eagle.cds.fasta

    wgd  Copyright (C) 2018 Arthur Zwaenepoel
    This program comes with ABSOLUTELY NO WARRANTY;
    """
    codeml_opts = {x.split("=")[0].strip(): x.split("=")[1].strip()
                   for x in codeml_options.split(",") if x != ""}
    ksd_(
            gene_families, sequences, output_directory, protein_sequences,
            tmp_dir, aligner, codeml='codeml',
            times=times, min_msa_length=min_msa_length,
            ignore_prefixes=ignore_prefixes, one_v_one=one_v_one,
            preserve=preserve, n_threads=n_threads,
            weighting_method=wm, pairwise=pairwise,
            max_pairwise=max_pairwise, **codeml_opts
    )


def ksd_(
        gene_families, sequences, output_directory, protein_sequences=None,
        tmp_dir=None, aligner='muscle', codeml='codeml', times=1,
        min_msa_length=100, ignore_prefixes=False, one_v_one=False,
        pairwise=False, preserve=False, n_threads=4,
        weighting_method='fasttree', max_pairwise=10000, **kwargs
):
    """
    Ks distribution construction pipeline. For usage in the ``wgd`` CLI.

    :param gene_families: gene families, i.e. tab separated paralogs or
        one-vs-one orthologs (see :py:func:`blast_`)
    :param sequences: CDS fasta files, if multiple (for constructing one-vs.-one
        ortholog distribution) then as a comma separated string
    :param output_directory: output directory
    :param protein_sequences: protein sequences (optional), by default CDS files
        are translated using the standard genetic code.
    :param tmp_dir: tmp directory name (optional)
    :param aligner: aligner to use
    :param codeml: path to codeml executable
    :param times: number of times to iteratively perform ML estimation of Ks,
        Ka and omega values.
    :param min_msa_length: minimum multiple sequence alignment length
    :param ignore_prefixes: ignore prefixes defined by '|' in gene IDs
    :param one_v_one: boolean, one-vs.-one ortholog analysis
    :param pairwise: run in pairwise mode
    :param preserve: boolean, preserve codeml output files, multiple sequence
        alignments and trees?
    :param async: use the async library for parallelization (not recommended)
    :param n_threads: number of threads to use
    :param weighting_method: weighting method (fasttree, phyml or alc)
    :param max_pairwise: maximum number of pairwise combinations a gene family
        may have. This effectively filters out families of size n where
        n*(n-1)/2 exceeds `max_pairwise`.
    :return: output file name
    """
    # lazy imports
    from wgd.ks_distribution import ks_analysis_paranome, ks_analysis_one_vs_one
    from wgd.viz import plot_selection
    from wgd.utils import uniq_id

    # software check
    software = [codeml, aligner]
    if weighting_method == 'fasttree':
        software.append('FastTree')
    elif weighting_method == 'phyml':
        software.append('phyml')
    if can_i_run_software(software) == 1:
        logging.error(
                'Could not run all required external software, exit here.')
        return 1

    # input check
    if not (gene_families and sequences):
        logging.error('No gene families or no sequences provided.')
        return 1

    if not tmp_dir:
        # unique tmp directory
        tmp_dir = os.path.join('.', 'ks_tmp.' + uniq_id())

    # get absolute paths before changing dir
    output_directory = os.path.abspath(output_directory)
    tmp_dir = os.path.abspath(tmp_dir)
    gene_families = os.path.abspath(gene_families)

    # some warnings
    if os.path.exists(output_directory):
        logging.warning('Output directory exists, will possibly overwrite')
    else:
        os.mkdir(output_directory)
    if os.path.exists(tmp_dir):
        logging.info('tmp directory exists, will try to resume analysis.')
    else:
        os.mkdir(tmp_dir)

    logging.debug('Reading CDS sequences')
    seq_list = [os.path.abspath(x) for x in sequences]
    cds_seqs = {}
    for seq_file in seq_list:
        cds_seqs.update(read_fasta(seq_file))

    # translate CDS file(s)
    if not protein_sequences:
        logging.info('Translating CDS file')
        protein_seqs = translate_cds(cds_seqs)

    else:
        logging.debug('Reading protein sequences')
        seq_list = [os.path.abspath(x) for x in
                    protein_sequences.strip().split(',')]
        protein_seqs = {}
        for seq_file in seq_list:
            protein_seqs.update(read_fasta(seq_file))

    # define a base name for the output files
    base = '_'.join([os.path.basename(x) for x in seq_list])

    # one-vs-one ortholog input
    if one_v_one:
        os.chdir(tmp_dir)  # chdir is necessary because codeml produces these
        # rub, rst and rst1 files
        logging.info('Started one-vs-one ortholog Ks analysis')
        results = ks_analysis_one_vs_one(
            cds_seqs, protein_seqs, gene_families,
            tmp_dir, output_directory, codeml,
            n_threads=n_threads, preserve=preserve,
            times=times,
            aligner=aligner,
            **kwargs
        )
        results.round(5).to_csv(os.path.join(
            output_directory, '{}.ks.tsv'.format(base)), sep='\t')

        logging.info('Generating plots')
        plot_selection(
            results, title=os.path.basename(gene_families),
            output_file=os.path.join(output_directory, '{}.ks.svg'.format(base))
        )

        logging.info('Done')
        return os.path.join(output_directory, '{}.ks.tsv'.format(base))

    # whole paranome ks analysis
    else:
        os.chdir(tmp_dir)  # change directory to the tmp dir, as codeml writes
        # non-unique file names to the working dir
        logging.info('Started whole paranome Ks analysis')
        results = ks_analysis_paranome(
            cds_seqs, protein_seqs, gene_families,
            tmp_dir, output_directory,
            codeml, preserve=preserve,
            times=times, aligner=aligner,
            ignore_prefixes=ignore_prefixes,
            n_threads=n_threads,
            min_length=min_msa_length,
            method=weighting_method,
            pairwise=pairwise,
            max_pairwise=max_pairwise,
            **kwargs
        )
        results.round(5).to_csv(os.path.join(
            output_directory, '{}.ks.tsv'.format(base)), sep='\t')

        logging.info('Generating plots')

        plot_selection(
            results, title=os.path.basename(gene_families),
            output_file=os.path.join(output_directory, '{}.ks.svg'.format(base))
        )

        logging.info('Done')
        return os.path.join(output_directory, '{}.ks.tsv'.format(base))


# CO-LINEARITY -----------------------------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('gff_file', default=None, type=click.Path(exists=True))
@click.argument('gene_families', default=None, type=click.Path(exists=True))
@click.option(
        '--ks_distribution', '-ks', default=None,
        help="paranome ks distribution tsv file (optional, see `wgd ks`)"
)
@click.option(
        '--output_dir', '-o', default='./wgd_syn', show_default=True,
        help='output directory'
)
@click.option(
        '--feature', '-f', default='gene', show_default=True,
        help="keyword for parsing the genes from the GFF file (column 3)"
)
@click.option(
        '--gene_attribute', '-a', default='ID', show_default=True,
        help="keyword for parsing the gene IDs from the GFF file (column 9)"
)
@click.option(
        '--min_length', '-l', default=250, show_default=True,
        help="minimum length of a genomic element (in numbers of genes) to be "
             "included in dotplot."
)
@click.option(
        '--ks_range', '-r', nargs=2, default=(0.05, 5), show_default=True,
        type=float, help='Ks range to use for colored dotplot'
)
def syn(
        gff_file, gene_families, ks_distribution, output_dir, feature,
        gene_attribute, min_length, ks_range
):
    """
    Co-linearity analyses.

    Requires I-ADHoRe 3.0 and a structural annotation in GFF format. Use the
    ``--feature`` and ``--gene_attribute`` flags to ensure correspondence
    between your CDS sequences and the right features in the GFF file!

    Example:

        wgd syn -f gene -a id -ks panda.ks panda.gff panda.paranome.mcl

    wgd  Copyright (C) 2018 Arthur Zwaenepoel
    This program comes with ABSOLUTELY NO WARRANTY;
    """
    syn_(
            gff_file, gene_families, output_dir, ks_distribution, feature,
            gene_attribute, min_length, ks_range
    )


def syn_(
        gff_file, families, output_dir, ks_distribution, feature='mRNA',
        gene_attribute='Parent', min_length=250, ks_range=(0.05, 5)
):
    """
    Co-linearity analysis with I-ADHoRe 3.0. For usage in the ``wgd`` CLI.

    :param gff_file: GFF3 annotation file (see the annotation files on PLAZA as
        an example)
    :param families: gene families as tab separated gene IDs, see
        :py:func:`blast_`
    :param output_dir: output directory
    :param ks_distribution: Ks distribution tsv file, see :py:func:`ks_`
    :param feature: keyword for entities of interest in the GFF file, e.g.
        'CDS' or 'mRNA'
    :param gene_attribute: attribute key for the gene ID in the GFF (9th
        column), e.g. 'ID' or 'Parent'
    :return: nothing at all
    """
    # lazy imports
    from wgd.colinearity import write_families_file, write_gene_lists
    from wgd.colinearity import write_config_adhore, run_adhore
    from wgd.colinearity import get_anchor_pairs, gff_parser
    from wgd.viz import plot_selection, syntenic_dotplot
    from wgd.viz import syntenic_dotplot_ks_colored

    # software check
    if can_i_run_software(['i-adhore']) == 1:
        logging.error('Could not run all software, exit here.')
        return 1

    # input check
    if not gff_file:
        logging.error('No gff file provided! Run `wgd syn --help` for usage '
                      'instructions.')
        return 1

    if not families:
        logging.error('No gene families provided! Run `wgd syn --help` for '
                      'usage instructions.')
        return 1

    if os.path.exists(output_dir):
        logging.warning(
                'Output directory already exists, will possibly overwrite')

    else:
        os.mkdir(output_dir)
        logging.info('Made output directory {0}'.format(output_dir))

    # parse the gff
    logging.info("Parsing GFF file")
    try:
        genome, all_genes = gff_parser(
                gff_file, feature=feature, gene_attribute=gene_attribute)
    except IndexError:
        logging.error('Invalid GFF file, number of columns != 9')
        logging.error('Aborted')
        return 1

    # generate necessary files for i-adhore
    logging.info("Writing gene lists")
    write_gene_lists(genome, os.path.join(output_dir, 'gene_lists'))

    logging.info("Writing families file")
    write_families_file(
            families, all_genes, os.path.join(output_dir, 'families.tsv'))

    logging.info("Writing configuration file")
    write_config_adhore(
            os.path.join(output_dir, 'gene_lists'),
            os.path.join(output_dir, 'families.tsv'),
            config_file_name=os.path.join(output_dir, 'adhore.conf'),
            output_path=os.path.join(output_dir, 'i-adhore-out'))

    # run i-adhore
    logging.info("Running I-ADHoRe 3.0")
    run_adhore(os.path.join(output_dir, 'adhore.conf'))

    # dotplot
    multiplicons = pd.read_csv(os.path.join(
            output_dir, 'i-adhore-out', 'multiplicons.txt'), sep='\t')
    logging.info('Drawing co-linearity dotplot')
    syntenic_dotplot(
            multiplicons, min_length=min_length, output_file=os.path.join(
                    output_dir,
                    '{}.dotplot.svg'.format(os.path.basename(families)))
    )

    # Ks distribution for anchors and Ks colored dotplot
    if ks_distribution:
        # input files
        anchor_points = pd.read_csv(os.path.join(
                output_dir, 'i-adhore-out', 'anchorpoints.txt'), sep='\t')
        ks_dist = pd.read_csv(ks_distribution, index_col=0, sep='\t')

        # output file names
        ks_out = os.path.join(output_dir, '{}.ks_anchors.tsv'.format(
                os.path.basename(families)))
        dotplot_out = os.path.join(output_dir, '{}.dotplot.ks.svg'.format(
                os.path.basename(families)))
        hist = os.path.join(output_dir, '{}.ks_anchors.svg'.format(
                os.path.basename(families)))

        # output and plots
        logging.info("Constructing Ks distribution for anchors")
        ks, anchors = get_anchor_pairs(anchor_points, ks_dist, ks_out)

        logging.info("Generating Ks colored (median Ks) dotplot")
        syntenic_dotplot_ks_colored(
                multiplicons, anchor_points, anchors, min_ks=ks_range[0],
                max_ks=ks_range[1], output_file=dotplot_out,
                min_length=min_length
        )

        logging.info("Generating histogram")
        plot_selection([ks, anchors], alphas=[0.2, 0.7], output_file=hist,
                       title=os.path.basename(families), weighted=False,
                       labels=['Whole paranome', 'Anchors'])

    logging.info("Done")


# KDE FITTING ------------------------------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument(
        'ks_distribution', type=click.Path(exists=True), default=None
)
@click.option(
        '--filters', '-f', nargs=3, type=float, default=(0., 300., 0.),
        help="Data frame filters (alignment identity, length and coverage)",
        show_default=True
)
@click.option(
        '--ks_range', '-r', nargs=2, default=(0, 3), show_default=True,
        type=float,
        help='Ks range to use for modeling'
)
@click.option(
        '--bandwidth', '-bw', default=None, show_default=True, type=float,
        help="Bandwidth for Gaussian KDE, by default Scott's rule is used"
)
@click.option(
        '--bins', '-b', default=25, show_default=True, type=int,
        help="Number of histogram bins."
)
@click.option(
        '--output_file', '-o', default="kde.svg", show_default=True,
        help='output file'
)
def kde(
        ks_distribution, filters, ks_range, bandwidth, bins, output_file
):
    """
    Fit a KDE to a Ks distribution.

    This accounts for boundary effects by applying reflection around the minimum
    Ks value, removing spurious peaks in low Ks regions. Please modify the
    bandwidth parameter if the KDE seems to be under- or overfitting.

    Note that `wgd viz` allows interactive plotting of KDEs also and might be
    more convenient for exploratory analysis.

    Note that histogram weighting is done after applying specified filters. Also
    note that mixture models are fitted to node-averaged (not weighted)
    histograms. Not supported for one-vs.-one ortholog Ks distributions (but see
    `wgd viz`).

    wgd  Copyright (C) 2018 Arthur Zwaenepoel
    This program comes with ABSOLUTELY NO WARRANTY;
    """
    kde_(ks_distribution, filters, ks_range, bandwidth, bins, output_file)


def kde_(ks_distribution, filters, ks_range, bandwidth, bins, output_file):
    """
    Fit a KDE to a Ks distribution.

    :param ks_distribution: Ks distribution file
    :param filters: alignment filters
    :param ks_range: Ks range
    :param bandwidth: bandwidth
    :param bins: number of histogram bins
    :param output_file: output file
    :return: nada
    """
    from wgd.modeling import filter_group_data, reflected_kde
    df = pd.read_csv(ks_distribution, index_col=0, sep='\t')
    df = filter_group_data(df, filters[0], filters[1], filters[2],
                           ks_range[0], ks_range[1])
    reflected_kde(df, ks_range[0], ks_range[1], bandwidth, bins, output_file)
    pass


# MIXTURE MODELING -------------------------------------------------------------
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument(
        'ks_distribution', type=click.Path(exists=True), default=None
)
@click.option(
        '--filters', '-f', nargs=3, type=float, default=(0., 300, 0.),
        help="Data frame filters (alignment identity, length and coverage)",
        show_default=True
)
@click.option(
        '--ks_range', '-r', nargs=2, default=(0.005, 3), show_default=True,
        type=float,
        help='Ks range to use for modeling'
)
@click.option(
        '--bins', '-b', default=50, show_default=True, type=int,
        help="Number of histogram bins."
)
@click.option(
        '--output_dir', '-o', default="wgd_mix", show_default=True,
        help='output directory'
)
@click.option(
        '--method', type=click.Choice(['gmm', 'bgmm']), default='gmm',
        show_default=True, help="mixture modeling method"
)
@click.option(
        '--components', '-n', nargs=2, default=(1, 4), show_default=True,
        help='range of number of components to fit'
)
@click.option(
        '--gamma', '-g', default=1e-3, show_default=True,
        help='gamma parameter for bgmm models'
)
@click.option(
        '--n_init', '-ni', default=1, show_default=True,
        help='number of k-means initializations'
)
@click.option(
        '--max_iter', '-mi', default=1000, show_default=True,
        help='maximum number of iterations'
)
def mix(
        ks_distribution, filters, ks_range, bins, output_dir, method,
        components, gamma, n_init, max_iter
):
    """
    Mixture modeling of Ks distributions.

    This provides means to fit Gaussian mixture models (GMMs) using an
    expectation-maximization algorithm (EM, method=gmm) and variational Bayes
    (VB, method=bgmm) algorithm. In the later case, regularization is performed
    with strength inversely proportional with the gamma parameter. A low gamma
    value will lead to less active components in the mixture, which can be
    observed by noting the weights for those components shrinking towards 0.
    Note that histogram weighting is done after applying specified filters.

    Please interpret mixture model results with caution, for more
    info, refer to https://wgd.readthedocs.io/en/latest/mix.html#a-note-on-mixtu
    re-models-for-ks-distributions

    Not supported for one-vs-one ortholog Ks distributions.

    wgd  Copyright (C) 2018 Arthur Zwaenepoel
    This program comes with ABSOLUTELY NO WARRANTY;
    """
    mix_(
            ks_distribution, filters, ks_range, method, components, bins,
            output_dir, gamma, n_init, max_iter
    )


def mix_(
        ks_distribution, filters, ks_range, method, components, bins,
        output_dir, gamma, n_init, max_iter
):
    """
    Mixture modeling tools.

    Note that histogram weighting is done after applying specified filters. Also
    note that mixture models are fitted to node-averaged (not weighted)
    histograms. Please interpret mixture model results with caution, for more
    info, refer to :ref:`note_on_gmms`.

    :param ks_distribution: Ks distribution data frame
    :param filters: alignment stats filters
    :param ks_range: Ks range used for models
    :param method: mixture modeling method, Bayesian/ordinary Gaussian mixtures
    :param components: number of components to use (tuple: (min, max))
    :param bins: number histogram bins for visualization
    :param output_dir: output directory
    :param gamma: gamma parameter for BGMM
    :param n_init: number of k-means initializations (best is kept)
    :param max_iter: number of iterations
    :return: nada
    """
    from wgd.modeling import filter_group_data, get_array_for_mixture, fit_gmm
    from wgd.modeling import inspect_aic, inspect_bic, plot_aic_bic
    from wgd.modeling import plot_all_models_gmm, get_component_probabilities
    from wgd.modeling import fit_bgmm, plot_all_models_bgmm

    # make output dir if needed
    if not os.path.exists(output_dir):
        logging.info("Making directory {}".format(output_dir))
        os.mkdir(output_dir)

    # prepare data frame
    logging.info("Preparing data frame")
    df = pd.read_csv(ks_distribution, index_col=0, sep='\t')
    df = filter_group_data(df, filters[0], filters[1], filters[2],
                           ks_range[0], ks_range[1])
    X = get_array_for_mixture(df)

    logging.info(" .. max_iter = {}".format(max_iter))
    logging.info(" .. n_init   = {}".format(n_init))

    # GMM method
    if method == "gmm":
        logging.info("Method is GMM, interpret best model with caution!")
        models, bic, aic, best = fit_gmm(
                X, components[0], components[1], max_iter=max_iter,
                n_init=n_init
        )
        inspect_aic(aic)
        inspect_bic(bic)
        logging.info("Plotting AIC & BIC")
        plot_aic_bic(aic, bic, components[0], components[1],
                     os.path.join(output_dir, "aic_bic.svg"))
        logging.info("Plotting mixtures")
        plot_all_models_gmm(models, X, ks_range[0], ks_range[1], bins=bins,
                            out_file=os.path.join(output_dir, "gmms.svg"))

    # BGMM method
    else:
        logging.info("Method is BGMM, weights are informative for best model")
        logging.info(" .. gamma    = {}".format(gamma))
        models = fit_bgmm(
                X, components[0], components[1], gamma=gamma,
                max_iter=max_iter, n_init=n_init
        )
        logging.info("Plotting mixtures")
        plot_all_models_bgmm(models, X, ks_range[0], ks_range[1], bins=bins,
                             out_file=os.path.join(output_dir, "bgmms.svg"))
        logging.warning("Method is BGMM, unable to choose best model!")
        logging.info("Taking model with most components for the component-wise"
                     "probability output file.")
        logging.info("To get the output file for a particular number of "
                     "components, run wgd mix again ")
        logging.info("with the desired component number as maximum.")
        best = models[-1]

    # save component probabilities
    logging.info("Writing component-wise probabilities to file")
    new_df = get_component_probabilities(df, best)
    new_df.round(5).to_csv(os.path.join(
            output_dir, "ks_{}.tsv".format(method)), sep="\t")


@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option(
        '--ks_distributions', '-ks', default=None,
        help="comma-separated Ks distribution csv files, or a directory containing"
             " these (see `wgd ks`)"
)
@click.option(
        '--alpha_values', '-a', default=None,
        help="comma-separated alpha values (optional)"
)
@click.option(
        '--colors', '-c', default=None,
        help="comma-separated colors (optional)"
)
@click.option(
        '--labels', '-l', default=None,
        help="comma-separated labels (for legend, optional)")
@click.option(
        '--hist_type', '-ht', default='barstacked', show_default=True,
        type=click.Choice(['barstacked', 'step', 'stepfilled']),
        help="histogram type"
)
@click.option(
        '--title', '-t', default='', show_default=True,
        help="plot title"
)
@click.option(
        '--output_file', '-o', default='wgd_hist.svg', show_default=True,
        help="output file"
)
@click.option(
        '--interactive', '-i', is_flag=True,
        help="interactive visualization with bokeh"
)
@click.option(
        '--filters', '-f', nargs=3, type=float, default=(0., 300, 0.),
        help="Data frame filters (alignment identity, length and coverage)",
        show_default=True
)
@click.option(
        '--ks_range', '-r', nargs=2, default=(0.005, 3), show_default=True,
        type=float,
        help='Ks range to use for modeling'
)
@click.option(
        '--bins', '-b', default=50, show_default=True, type=int,
        help="Number of histogram bins."
)
@click.option(
        '--weighted', is_flag=True,
        help="Plot node-weighted histograms instead of node-averaged [NA if "
             "using --interactive]."
)
def viz(
        ks_distributions, alpha_values, colors, labels, hist_type, title,
        output_file, interactive, filters, ks_range, bins, weighted
):
    """
    Plot histograms/densities (interactively).

    Requires a running bokeh server for interactive visualization.
    Run a bokeh serve instance (in the background) with `bokeh serve &`.

    more information about bokeh: https://bokeh.pydata.org/en/latest/index.html

    wgd  Copyright (C) 2018 Arthur Zwaenepoel
    This program comes with ABSOLUTELY NO WARRANTY;
    """
    viz_(
            ks_distributions, alpha_values, colors, labels, hist_type, title,
            output_file, filters, ks_range, bins, interactive, weighted
    )


def viz_(
        ks_distributions, alpha_values, colors, labels, hist_type, title,
        output_file, filters, ks_range, bins, interactive=False, weighted=False
):
    """
    Plot (stacked) histograms (interactively). Add option to plot node-weighted
    histograms in the same fashion.

    :param ks_distributions: a directory with ks distributions (other files are
        ignored) or a comma-separated string of file names
    :param alpha_values: alpha values for the different distributions (in the
        same order). Only relevant for non-interactive visualization.
    :param colors: as in ``alpha_values`` but for colors
    :param labels: as in ``alpha_values`` but for legend labels (by default the
        file names are used), this is also relevant for the interactive bokeh
        visualization (as opposed to ``alpha_values`` and ``colors``.
    :param hist_type: histogram type (matplotlib), either 'barstacked', 'step'
        or 'stepfilled'.
    :param title: plot title
    :param output_file: output file name
    :param interactive: render an interactive bokeh plot. This makes some of the
        above arguments redundant
    :return: nada
    """
    from wgd.viz import plot_selection

    # basic input check
    if not ks_distributions:
        logging.error('No Ks distributions provided, run `wgd viz -h` for '
                      'usage instructions')
        logging.error('You have to provide one or more computed Ks '
                      'distributions! See for example `wgd ks -h`')
        return 1

    # directory provided
    if os.path.isdir(ks_distributions):
        ks_distributions = [os.path.join(ks_distributions, x) for x in
                            os.listdir(ks_distributions)]

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
            _ = d[['Ks']]
            d = d.fillna(0)
            dists_files.append(i)
            dists.append(d)
        except:
            logging.info('Not a Ks distribution: {}'.format(i))

    # interactive bokeh visualization
    if interactive:
        from wgd.viz import histogram_bokeh
        histogram_bokeh(dists_files, labels)
        return

    # normal matplotlib plots
    # features
    if alpha_values:
        alpha_values = [float(x) for x in alpha_values.split(',')]
        if len(alpha_values) != len(dists):
            logging.error('Please provide as much alpha values as there are '
                          'distributions')
    if colors:
        colors = [x for x in colors.split(',')]
        if len(colors) != len(dists):
            logging.error(
                    'Please provide as much colors as there are distributions')
    if not labels:
        labels = dists_files
        if len(labels) != len(dists):
            logging.error(
                    'Please provide as much labels as there are distributions')
    else:
        labels = labels.split(',')

    # make plot
    logging.info('Plotting Ks distributions overlay')
    plot_selection(dists, alphas=alpha_values, colors=colors, labels=labels,
                   output_file=output_file, title=title, histtype=hist_type,
                   filters=filters, ks_range=ks_range, bins=bins,
                   weighted=weighted)
    return


@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('sequences', default=None)
@click.argument('output_dir', default=None)
@click.option(
        '--gff_file', '-gff', default=None,
        help="GFF file for co-linearity analysis"
)
@click.option(
        '--n_threads', '-n', default=4, show_default=True,
        help="number of threads to use"
)
def wf1(sequences, output_dir, gff_file, n_threads):
    """
    Standard workflow whole paranome Ks.

    Parameters used are the default parameters in `wgd mcl`, `wgd ksd` and if
    relevant `wgd syn`.

    Example:

        wgd pipeline_1 -gff snail.gff -n 16 snail.fasta snail_wgd_out

    wgd  Copyright (C) 2018 Arthur Zwaenepoel
    This program comes with ABSOLUTELY NO WARRANTY;
    """
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # wgd blast
    blast_dir = os.path.join(output_dir, 'wgd_mcl')
    mcl_out = blast_mcl(cds=True, mcl=True, sequences=sequences,
                        n_threads=n_threads, output_dir=blast_dir)

    # wgd ks
    ks_dir = os.path.join(output_dir, 'wgd_ksd')
    ks_results = ksd_(gene_families=mcl_out, sequences=[sequences],
                      n_threads=n_threads, output_directory=ks_dir)

    # wgd syn
    if gff_file:
        syn_dir = os.path.join(output_dir, 'wgd_syn')
        syn_(gff_file=gff_file, families=mcl_out, output_dir=syn_dir,
             ks_distribution=ks_results)

    return


@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('sequences', default=None, nargs=2)
@click.argument('output_dir', default=None)
@click.option(
        '--n_threads', '-n', default=4, show_default=True,
        help="number of threads to use"
)
def wf2(sequences, output_dir, n_threads):
    """
    Standard workflow one-vs-one ortholog Ks.

    Parameters used are the default parameters in `wgd mcl`, and `wgd ksd`.

    Example:

        wgd pipeline_2 -n 8 snail.fasta whale.fasta snail_vs_whale_ks_out

    wgd  Copyright (C) 2018 Arthur Zwaenepoel
    This program comes with ABSOLUTELY NO WARRANTY;
    """
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # wgd blast
    blast_dir = os.path.join(output_dir, 'wgd_mcl')
    ovo_out = blast_mcl(cds=True, one_v_one=True, mcl=False,
                        sequences=sequences,
                        n_threads=n_threads,
                        output_dir=blast_dir)

    # wgd ks
    ks_dir = os.path.join(output_dir, 'wgd_ksd')
    ksd_(gene_families=ovo_out, sequences=sequences, n_threads=n_threads,
         output_directory=ks_dir, one_v_one=True)

    return


if __name__ == '__main__':
    cli()
