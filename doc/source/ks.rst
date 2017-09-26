Ks distribution analysis
************************

Ks distribution construction
============================

.. automodule:: wgd.ks_distribution
    :members:
    :private-members:
    :special-members: __init__

Modeling of Ks distributions
============================

.. automodule:: wgd.modeling
    :members:

Example: performing Ks analysis
===============================

To perform a full Ks analysis on a gene family file, one can do something like this::

    from wgd import ks_distribution

    ks_dist = KsDistribution(species='Arabidopsis thaliana', gene_pattern='AT.+',
                             gene_families='gene_families.txt', nucleotide='nucleotide.fasta')
    ks_dist.ks_analysis_parallel()
    ks_dist.mixture_modeling()
    ks_dist.kde()
    ks_dist.write_report()

In fact, this is the pipeline as implemented in the command line utility ``ks.py`` (see :ref:`cli`)