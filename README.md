[![Documentation Status](https://readthedocs.org/projects/wgd/badge/?version=latest)](http://wgd.readthedocs.io/en/latest/?badge=latest)
[![Hosted](https://img.shields.io/badge/hosted-singularity--hub-blue.svg)](https://singularity-hub.org/collections/2097)


Copyright (C) 2018 Arthur Zwaenepoel

VIB/UGent center for plant systems biology -
Bioinformatics & evolutionary genomics group https://www.vandepeerlab.org/

# wgd - simple command line tools for the analysis of ancient whole-genome duplications

## Installation

Python package and command line interface (CLI) for the analysis of
whole-genome duplications (WGDs). Tested with Python3 on Linux. If you don't have
python or pip installed a simple `sudo apt-get install python3-pip` should do.

To install, simply run 

```
git clone https://github.com/arzwa/wgd.git
cd wgd
pip install --user .
```

Note that depending on your python installation and whether you're in a
virtualenv, ``pip`` may default either to ``pip2`` or ``pip3``. If the
above installation step fails, please try to use ``pip3`` instead of
``pip``.

For the command line interface, upon installation run

    $ wgd

to get a list of the available commands. To get usage instructions for
a command (e.g. `ksd`) run

    $ wgd ksd --help

For **external software** requirements: please consult the relevant section
in the [docs](https://wgd.readthedocs.io/en/latest/index.html#external-software)

## Quick start

The main aim of `wgd` is computing whole-paranome and one-vs.-one ortholog Ks
distributions. For a whole-paranome distribution of a CDS sequence fasta file,
the minimal commands are:

    $ wgd dmd ath.cds.fasta
    $ wgd ksd wgd_dmd/ath.cds.fasta.mcl ath.cds.fasta

For one-vs.one orthologs the minimal commands are

    $ wgd dmd ath.cds.fasta vvi.cds.fasta
    $ wgd ksd wgd_dmd/ath1000.fasta_vvi1000.fasta.rbh ath.cds.fasta vvi.cds.fasta

For more information and these methods and other tools implemented in `wgd`,
please consult the [docs](https://wgd.readthedocs.io/en/latest/).

## NEW: diamond support

This has not been thoroughly tested yet, but may be of interest.
[Diamond](https://github.com/bbuchfink/diamond) is a sequence similarity search
tool that is much faster than Blast and performance similar on other measures.
If you have diamond installed, you can use `wgd dmd` for whole-paranome or
one-vs.-one ortholog delineation, e.g.

    $ wgd dmd ath.cds.fasta

for a whole-paranome, and

    $ wgd dmd ath.cds.fasta vvi.cds.fasta

for one-vs.-one orthologs.

This is part of a new `wgd` version that is forthcoming. In the future this will
replace `wgd mcl`. Please let me know of any bugs!

## Singularity container

A singularity container is available for ``wgd``, allowing all to use
all tools in ``wgd`` except ``wgd syn``, without having to install all
required software on your system. To install Singularity follow
the instructions [here](https://www.sylabs.io/docs/)

If you have singulaity installed (and you're in the virtual machine when
running on Windows or Mac), you can run the following to get the container

    singularity pull --name wgd.simg shub://arzwa/wgd

Then you can use ``wgd`` as follows

    singularity exec wgd.simg wgd <command>

## Notes

**Bug tracking:** If the program crashes, exits unexpectedly or some
unexpected results are obtained, please run it again with the
``--verbosity debug`` flag *before* the subcommand of interest (*e.g.*
``wgd --verbosity debug ksd gf.mcl cds.fasta``). If the anomaly persists,
please open an issue on this GitHub site.

**Note on input data:** while the input data is rather straightforward
(a CDS fasta file will do for most analyses) it may be of interest that
the wgd suite was extensively tested with data from the PLAZA platform,
so for examples of the right input data formats (in particular CDS fasta
files for sequence data and GFF files for structural annotation), please
have a look [there](https://bioinformatics.psb.ugent.be/plaza/versions/plaza_v4_dicots/download/).
It is generally advised not to include pipe characters (`|`) in your gene
IDs, since these can have special meanings in certain parts of `wgd`.

**Note on virtualenv:** you can install wgd in a _virtual environment_
(using [`virtualenv`](https://virtualenv.pypa.io/en/stable/)). If you
would however encounter problems with running the executable directly
(e.g. `wgd --help` doesn't work) you can circumvent this by directly
calling the CLI, using `python3 ./wgd_cli.py --help` (assuming you are
currently in the directory where you cloned wgd).

## Citation

Please cite us at https://doi.org/10.1093/bioinformatics/bty915

```
Zwaenepoel, A., and Van de Peer, Y. wgd - simple command line tools for the analysis of ancient whole genome duplications. Bioinformatics., bty915, https://doi.org/10.1093/bioinformatics/bty915
```

For citation of the tools used in wgd, please consult the documentation at
https://wgd.readthedocs.io/en/latest/index.html#citation.
