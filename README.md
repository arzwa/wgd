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

## Third party software

`wgd` requires the following third party software (preferably these
should also be in the `PATH` environment variable):

- `blast`, from which it uses the `blastp` and `makeblastdb` commands,
`sudo apt-get install ncbi-blast+` will often suffice for installation
- `muscle` [Edgar2004] for multiple sequence alignment. `muscle` can be \
retrieved from http://www.drive5.com/muscle/.
- `codeml` from the PAML software package [Yang]. PAML can be downloaded
from the following link: http://abacus.gene.ucl.ac.uk/software/paml.html
- `I-ADHoRe 3.0` [Proost2012] for co-linearity analysis. `I-ADHoRe 3.0`
can be obtained from http://bioinformatics.psb.ugent.be/beg/tools/i-adhore30

## Example: K<sub>S</sub> distribution from CDS input

First we infer paralogous families using all-vs-all `blastp` and `mcl`:

    $ wgd blast --cds --mcl -s aptenodytes_patagonicus.cds.fasta -o ./penguin.wgd_out

Then we can construct a K<sub>S</sub> distribution using `codeml`:

    $ wgd ks -gf penguin.wgd_out/out.mcl -s aptenodytes_patagonicus.cds.fasta -o ./penguin.wgd_out

If we have structural annotation data available, we can look for acnhor
pairs and make a K<sub>S</sub> distribution for the obtained co-linear
gene pairs:

    $ wgd coll -ks penguin.wgd_out/all.csv -gff annotation.gff -gf penguin.wgd_out/out.mcl -o ./penguin.wgd_out

That's it! Quite easy no?





