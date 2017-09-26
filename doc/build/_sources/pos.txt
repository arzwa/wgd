Positive selection screening
****************************

.. automodule:: wgd.positive_selection
    :members:
    :private-members:
    :special-members: __init__

Example: performing screen
==========================

To perform a full positive selection screen for a gene family file and a species of interest,
one can do something like this::

    from wgd import positive_selection

    pos = PositiveSelection(species='Arabidopsis thaliana', gene_pattern='AT.+', gene_families='OrthoMCL_out.txt',
                            nucleotide='nucleotide.fasta', output_dir='./test_w', tmp_dir='./tmp')
    pos.positive_selection()

In fact, this is the pipeline as implemented in the command line utility ``wgd pos`` (see :ref:`cli`)