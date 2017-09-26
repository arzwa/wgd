# Whole genome duplication analysis

Python (3.5) package and command line interface (CLI) for the analysis
of whole genome duplications (WGDs).

To install: clone the repo, navigate to it and type

    $ pip install .

For the command line interface, upon installation run

    $ wgd

to get a list of the available commands. To get usage instructions for
a command (e.g. `ks`) run

    $ wgd ks --help

To use as a Python package as well as to find additional documentation
for the CLI, please consult the **extensive documentation** in the `doc`
directory or at https://arzwa.github.io/wgd/index.html.


## 1. _K<sub>S</sub>_ distributions

### Input data

_K<sub>S</sub>_ distribution analysis requires the following input data:

1. CDS sequences in fasta format
2. Output directory

If gene families/all-`versus`-all Blast results are provided separately, please take note of the following comments:

- All-_vs._-all Blast output should be formatted as
`gene_1 [TAB] gene_2 [TAB]  e-value`
- The gene families file should contain one family per line with tab
separated member genes.
- Optionally either a prefix or regular expression pattern can be provided
using the `--prefix` or `--regex` flag respectively. These parameters are
used to parse out the correct genes from the input all-_vs._-all Blast
results or gene families.
- Use `--prefix` when gene IDs from the nucleotide
fasta are prefixed with a species identifier in the all-_vs._-all/gene
families file, _e.g._ use `arath` when the gene looks like `arath|AT1G50700`
in the all-_vs._-all output, and like `AT1G50700` in the CDS fasta file input.
- Use `--regex` when you want to use only genes that match a particular
pattern (use Python/Perl regex patterns). This is particularly useful
when using multi-species gene families as input. For example use
`AT[1-9].+` for genes like `AT1G50700`, `AT4G09800` _etc._.

Note that if only genes of interest are included in the all-_vs._-all/gene
families file and they are identical with those in the fasta file, these
options are not needed.


### Usage

#### Command line utility

_K<sub>S</sub>_ distribution analysis can be performed using the command
line interface provided for the `wgd` package.
To use as a command line utility run `wgd ks --help` for further usage
instructions.

Basically usage of the command line utility is as follows:

    $ wgd ks [OPTIONS] <cds.fasta> <output directory>

The full help message reads as follows:

    $ wgd ks --help
    Usage: wgd ks_ [OPTIONS] CDS_FASTA_FILE OUTPUT_DIRECTORY

      Construct a Ks distribution. Enables running full analysis pipeline for a
      species of interest.

      ARGUMENTS (MINIMAL INPUT):
          - cds_fasta_file:
              File with CDS sequences in fasta format.
          - output_directory:
              Output directory name, species name is recommended (for plot titles etc.)

      OUTPUT:
          - csv files and histograms with Ks, Ka, w analysis output
          - modeling plots (Bayesian Gaussian mixture/kde)
          - html report

    Options:
      -ava, --all_vs_all TEXT         File with all-versus-all Blast results. Will
                                      e used for MCL clustering. Please format
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

###### Example: only CDS input
Using only a CDS file as input, ``wgd ks`` will perform all-`vs.`-all Blastp, MCL clustering and Ks analysis, `e.g.`:

    $ wgd ks aptenodytes_patagonicus.cds.fasta ./penguin_ks_out

This will take a while, the output will be in the directory `./penguin_ks_out`, in this case,
temporary file will be in `./penguin_ks_out.tmp`.

###### Example: multi species gene families
Using multi species gene families, both a gene families file
(_e.g._ OrthoMCL output) and a gene pattern (species specific regex
pattern) are required. For example with the data included in the
`example` directory:

    $ wgd ks --regex 'AT.+' -gf ./example/input/gene_families.txt ./example/input/cds.fasta ./example/output

This gives as output the results in a `csv` file, `png` images of the
histograms, plots of the mixture modeling and kernel density estimation
analysis and a summary report file in `html` format. The main file of
interest is this `html` file which gathers the results with some extra
information.

###### Example: from all-vs-all Blast results
When using all-vs-all Blast results as input data, first MCL clustering
using `mcl` is performed (should be installed on your system, it is
available from most Linux distro's). A default inflation factor of 2 is
used (can be modified with the `-I` flag). An example command:

    $ wgd --verbose=debug ks -ava all_v_all_blast_results.tsv  --prefix 'arat' cds.fasta output_dir

Note the usage of the verbose flag, which works for all commands in
`wgd`.

#### Python package

To use as a python package, just import the `wgd` module as usual. 
To perform the pipeline as implemented in `wgd ks`
manually in python do something like this:

	$ >>> ks_dist = KsDistribution(species=<species>, gene_families=<gene_families>,
	                               nucleotide=<nucleotide_sequences>)
	$ >>> ks_dist.ks_analysis_parallel()
	$ >>> ks_dist.mixture_modeling()
	$ >>> ks_dist.kde()
	$ >>> ks_dist.write_report()

Extensive documentation can be found in the `docs` directory.

### Third party software

The _K<sub>S</sub>_ distribution analysis uses Muscle (tested with
v3.8.31) for multiple sequence alignment and codeml (PAML package) for
_K<sub>S</sub>_ value estimation. It is recommended to add these
programs to the system PATHs.

