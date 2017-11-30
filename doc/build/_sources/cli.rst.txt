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
---------------------------------------------

Here  recipe is shown for |Ks| distribution construction using the provided commands in th ``wgd``
command line interface. For this analysis we could also use the provided pipeline ``wgd pipeline_1``
command, which provides less control, but is quite convenient in usage for this standard workflow.

1) Full control
~~~~~~~~~~~~~~~
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

2) Pipeline
~~~~~~~~~~~

The same cn be accomplished using the ``pipeline_1`` command::

    $ wgd pipeline_1 -gff annotation.gff penguin.cds.fasta penguin.wgd_out

However this dos not provide any further control, as opposed to the subcommands listed above.
Note that this ability to for more detailed control was not illustrated above.

More detailed documentation of the main ``wgd`` commands can be obtained here:

.. toctree::
   :maxdepth: 1

   wgd_cli <wgd_cli>


.. |Ks| replace:: K\ :sub:`S`
.. |Ka| replace:: K\ :sub:`A`