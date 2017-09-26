.. _cli:

Command line interface
**********************

The command line interface (CLI) and all command line utilities are written in python3.5 using the click library.
Upon installation, all available command line utilities bundled in the CLI can be listed by running::

    $ wgd

Which should, if installation went correctly, render the following response::

    Usage: wgd [OPTIONS] COMMAND [ARGS]...

      Welcome to the wgd command line interface!

                             _______
                             \  ___ `'.
             _     _ .--./)   ' |--.\  \
       /\    \\   ///.''\\    | |    \  '
       `\\  //\\ //| |  | |   | |     |  '
         \`//  \'/  \`-' /    | |     |  |
          \|   |/   /("'`     | |     ' .'
           '        \ '---.   | |___.' /'
                     /'""'.\ /_______.'/
                    ||     ||\_______|/
                    \'. __//
                     `'---'


    Options:
      --verbose [silent|info|debug]  Verbosity level, default = info.
      --help                         Show this message and exit.

    Commands:
      ks   Construct a Ks distribution.
      ls   List available commands with info.
      pos  Exploratory positive selection analysis.


To get a description for a command of interest (`e.g.` ls) run::

    $ wgd ls --help

Note that the verbose flag can be set before the subcommands::

    $ wgd --verbose [silent|info|debug] [COMMAND]

Which sets the verbosity for the logging.

Ks distribution analysis
========================

Command line utility for Ks, Ka and Ka/Ks analysis for whole paranomes.
For usage instructions run the following::

    $ wgd ks --help

Which will show::

    Usage: wgd ks [OPTIONS] CDS_FASTA_FILE OUTPUT_DIRECTORY

      Construct a Ks distribution. Enables running full analysis pipeline for a
      species of interest.

      ARGUMENTS (MINIMAL INPUT):
          - CDS_FASTA_FILE:
              File with CDS sequences in fasta format.
          - OUTPUT_DIRECTORY:
              Output directory name, species name is recommended (for plot titles etc.)

      OUTPUT:
          - csv files and histograms with Ks, Ka, w analysis output
          - modeling plots (Bayesian Gaussian mixture/kde)
          - html report

    Options:
      -ava, --all_vs_all TEXT         File with all-versus-all Blast results. Will
                                      be used for MCL clustering. Please format
                                      your file as `gene 1 [TAB] gene 2 [TAB]
                                      e-value`. Note that you can also provide
                                      precomputed gene families using the
                                      `--gene_families` flag.
      -gf, --gene_families TEXT       File with precomputed gene families, e.g.
                                      from MCL clustering of all-vs-all Blast
                                      results,OrthoMCL, INPARANOID, OrthoFinder
                                      etc. Please format your file with one tab
                                      separated gene family per line. Note that
                                      you can also provide raw all-vs-all Blast
                                      results using the `--all_vs_all` flag.
      -p, --prefix TEXT               Prefix for parsing out the genes of
                                      interest, for example 'orysa' for genes
                                      named like `orysa|OS1G11100`. Note that you
                                      can also provide a regex pattern
                                      (`--regex`). Note that if neither a prefix
                                      nor a regex pattern is provided, all genes
                                      are assumed to be relevant. (optional)
      -r, --regex TEXT                Regular expression pattern (Python/Perl
                                      type) for parsing out genes of interest.
                                      Especially useful when multi-species gene
                                      families are provided. Note that if neither
                                      a prefix nor a regex pattern is provided,
                                      all genes are assumed to be relevant.
                                      (optional)
      -I, --inflation_factor FLOAT    Inflation factor for MCL clustering, when
                                      blast results provided. (Default = 2)
      -ps, --protein_sequences TEXT   Protein sequences fasta file. Optional since
                                      by default the CDS file will be translated.
      -tmp, --tmp_dir TEXT            Path to store temporary files. (Default =
                                      ./)
      -m, --muscle TEXT               Path to muscle executable, not necessary if
                                      in PATH environment variable.
      -c, --codeml TEXT               Path to codeml executable, not necessary if
                                      in PATH environment variable.
      -t, --times INTEGER             Number of times to perform ML estimation
                                      (for more stable estimates). (Default = 1)
      --preserve / --no_preserve      Do you want to keep the multiple sequence
                                      alignment and codeml output? (Default =
                                      False)
      --prompt / --no_prompt          Prompt for directory clearing? (Default =
                                      True)
      -mix, --mixture [no|bayesian|gaussian]
                                      Mixture modeling method (requires
                                      sklearn.mixture). Set to no if not desired.
      --kde / --no_kde                Perform mixture modeling (requires
                                      sklearn.neighbors.KernelDensity)? (Default =
                                      True)
      -h, --help                      Show this message and exit.

Ks distribution analysis requires the following input data:

1. CDS sequences in fasta format
2. Output directory

If gene families or all-`versus`-all Blast results are provided separately, please take note of the following comments:

* All-`vs.`-all Blast output should be formatted as ``gene_1 [TAB] gene_2 [TAB]  e-value``
* The gene families file should contain one family per line with tab separated member genes.
* Optionally either a prefix or regular expression pattern can be provided using the ``--prefix`` or ``--regex`` flag respectively. These parameters are used to parse out the correct genes from the input all-`vs.`-all Blast results or gene families.
* Use ``--prefix`` when gene IDs from the nucleotide fasta are prefixed with a species identifier in the all-`vs.`-all/gene families file, `e.g.` use ``arath`` when the gene looks like ``arath|AT1G50700`` in the all-`vs.`-all output, and like ``AT1G50700`` in the CDS fasta file input.
* Use ``--regex`` when you want to use only genes that match a particular pattern (use Python/Perl regex patterns). This is particularly useful when using multi-species gene families as input. For example use ``AT[1-9].+`` for genes like ``AT1G50700``, ``AT4G09800`` `etc.`.

Note that if only genes of interest are included in the all-`vs.`-all/gene
families file and they are identical with those in the fasta file, these
options are not needed.

Example
-------
Using only a CDS file as input, ``wgd ks`` will perform all-`vs.`-all Blastp, MCL clustering and Ks analysis, `e.g.`::

    $ wgd ks aptenodytes_patagonicus.cds.fasta ./penguin_ks_out

Using *multi species gene families*, both a gene families file (`e.g.` OrthoMCL output) and a gene pattern (species specific regex
pattern) are required. For example with the data included in the `example` directory::

    $ wgd ks --regex 'AT.+' -gf ./example/input/gene_families.txt ./example/input/cds.fasta ./example/output

When using *all-vs-all Blast results* as input data, first MCL clustering
using ``mcl`` is performed (should be installed on your system, it is
available from most Linux distro's). A default inflation factor of 2 is
used (can be modified with the ``-I`` flag). An example command::

    $ wgd --verbose=debug ks -ava all_v_all_blast_results.tsv  --prefix 'arat' cds.fasta output_dir

Note the usage of the verbose flag, which works for all commands in ``wgd``.

Positive selection screen
=========================

Command line utility for running a positive selection screen for a species of interest.
For usage instructions run the following::

    $ wgd pos --help

Ks distribution analysis requires the following input data:

* Gene families file
* Nucleotide fasta (CDS) (one file with all sequences for all species, just ``cat >>`` your separate files (UNIX))
* Unique gene pattern (e.g. 'AT.+')

Example of a run on `Populus trichocarpa`::

    $ wgd pos -gd populus_descriptions.tsv 'PT.+' orthoMCL_out.txt nucleotide.fasta ./output_dir
