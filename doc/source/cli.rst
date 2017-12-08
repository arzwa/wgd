.. _cli:

Command line interface
**********************

The command line interface (CLI) and all command line utilities are written in python3.5 using the click library.
Upon installation, all available command line utilities bundled in the CLI can be listed by running::

    $ wgd --help

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

      Arthur Zwaenepoel - 2017

    Options:
      --verbose [silent|info|debug]  Verbosity level, default = info.
      -h, --help                     Show this message and exit.

    Commands:
      blast       All-vs.-all blastp (+ MCL) analysis.
      hist        Plot (stacked) histograms.
      ks          Ks distribution construction.
      mix         Mixture modeling of Ks distributions.
      pipeline_1  Standard workflow whole paranome Ks.
      pipeline_2  Standard workflow one-vs-one ortholog Ks.
      syn         Co-linearity analyses.

Note that the verbose flag can be set before the subcommands::

    $ wgd --verbose [silent|info|debug] [COMMAND]

Which sets the verbosity for the logging.

Examples
********

Whole paranome |Ks| distribution construction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here a recipe is shown for |Ks| distribution construction using the provided commands in the ``wgd``
command line interface. For this analysis we could also use the provided pipeline ``wgd pipeline_1``
command, which provides less control, but is quite convenient in usage for this standard workflow.

**1) Full control**

First we infer paralogous families using all-vs-all `blastp` and `mcl`::

    $ wgd blast --cds --mcl -s penguin.cds.fasta -o ./

Then we can construct a |Ks| distribution for the paranome we just generated
in the file ``penguin.cds.fasta.mcl`` using `wgd ks`::

    $ wgd ks -gf penguin.cds.fasta.mcl -s penguin.cds.fasta -o ./

If we have adequate structural annotation data available, we can look for anchor
pairs and make a |Ks| distribution for the obtained co-linear
gene pairs::

    $ wgd coll -ks penguin.cds.fasta.ks.tsv -gff annotation.gff -gf penguin.cds.fasta.mcl -o ./

That's it! Quite easy no?

**2) Pipeline**

The same can be accomplished using the ``pipeline_1`` command::

    $ wgd pipeline_1 -gff annotation.gff penguin.cds.fasta penguin.wgd_out

However this dos not provide any further control, as opposed to the subcommands listed above.
Note that this ability to for more detailed control was not illustrated above.

One-`vs.`-one ortholog |Ks| distribution construction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Again, the longer way, with more customization options::

   $ wgd blast --cds --one_v_one -s ath.cds.fasta,aly.cds.fasta
   $ wgd ks --one_v_one -s ath.cds.fasta,aly.cds.fasta -gf aly.cds.fasta_ath.cds.fasta.ovo.tsv

And the shorter option::

   $ wgd pipeline_2 -n 8 ath.cds.fasta,aly.cds.fasta ./aly_ath_one-vs_one_out

Intragenomic co-linearity analysis and anchor-base |Ks| distributions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Example command::

   $ wgd syn -gff ath.gff -gf ath.paranome.mcl -o ./wgd_syn

Mixture modeling of |Ks| distributions and extraction of WGD-specific paralogs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``wgd mix`` command makes it easy to perform mixture model analysis of |Ks|
distributions. Also, it provides the possibility to extract paralogs for a given
mixture component of interest, which can than be used for functional analysis
(_e.g._ GO enrichment, see ``example`` directory in the ``wgd`` repository) and
relative or absolute phylogenomic data.

Example command::

   $ wgd mix --method bgmm -ks ./ath.ks.tsv -n 1,4 -s ath.cds.fasta --pairs

Will fit 1 to 4 components using a Bayesian Gaussian mixture model to the |Ks|
distribution in the file ``ath.ks.tsv`` and will subsequently write fasta files
for the pairs for each component of the best fitting mixture (which could then be
used for dating studies).

Interactive visualization of |Ks| distributions and kernel density estimates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Make sure you have a ``bokeh`` server running (https://bokeh.pydata.org/en/latest/) ::

   $ bokeh serve &

If you have some |Ks| distributions in a directory, say ``cool_plants_ks`` you can
plot them interactively using ``wgd viz``::

   $ wgd viz -ks ./cool_plants_ks -i

Where the ``-i`` flag triggers interactive mode. In this mode you can dynamically alter
the range, scale, colors, _etc._ as well hide and show different distributions and plot
kernel density estimates.

Command line interface
~~~~~~~~~~~~~~~~~~~~~~

More detailed documentation of the main ``wgd`` commands can be obtained here:

.. toctree::
   :maxdepth: 1

   wgd_cli <wgd_cli>


.. |Ks| replace:: K\ :sub:`S`
.. |Ka| replace:: K\ :sub:`A`