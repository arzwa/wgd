.. _cli:

Command line interface
**********************

The command line interface (CLI) and all command line utilities are written in python3.5 using the click library.
Upon installation, all available command line utilities bundled in the CLI can be listed by running::

    $ wgd --help

Note that the verbose flag can be set before the subcommands::

    $ wgd --verbose [silent|info|debug] [COMMAND]

Which sets the verbosity for the logging.

Ks distribution analysis from a CDS fasta file - example
========================================================

First we infer paralogous families using all-vs-all `blastp` and `mcl`::

    $ wgd blast --cds --mcl -s aptenodytes_patagonicus.cds.fasta -o ./penguin.wgd_out

Then we can construct a K<sub>S</sub> distribution using `codeml`::

    $ wgd ks -gf penguin.wgd_out/out.mcl -s aptenodytes_patagonicus.cds.fasta -o ./penguin.wgd_out

If we have structural annotation data available, we can look for acnhor
pairs and make a K<sub>S</sub> distribution for the obtained co-linear
gene pairs::

    $ wgd coll -ks penguin.wgd_out/all.csv -gff annotation.gff -gf penguin.wgd_out/out.mcl -o ./penguin.wgd_out

That's it! Quite easy no?
