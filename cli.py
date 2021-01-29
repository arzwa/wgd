#!/usr/bin/python3
import click
import logging
import sys
import os
import warnings
import pandas as pd
import subprocess as sp
import pkg_resources  # part of setuptools
from rich.logging import RichHandler
__version__ = pkg_resources.require("wgd")[0].version
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")


# CLI entry point
@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--verbosity', '-v', type=click.Choice(['info', 'debug']),
    default='info', help="Verbosity level, default = info.")
def cli(verbosity):
    """
    wgd - Copyright (C) 2018-2020 Arthur Zwaenepoel\n
    Contact: arzwa@psb.vib-ugent.be
    """
    logging.basicConfig(
        format='%(message)s',
        handlers=[RichHandler()],
        datefmt='%H:%M:%S',
        level=verbosity.upper())
    logging.info("This is wgd v{}".format(__version__))
    pass


# Diamond and gene families
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('sequences', nargs=-1, type=click.Path(exists=True))
@click.option('--outdir', '-o', default='wgd_dmd', show_default=True,
    help='output directory')
@click.option('--tmpdir', '-t', default=None, show_default=True,
    help='tmp directory')
@click.option('--inflation', '-I', default=2.0,
    help="inflation factor for MCL")
@click.option('--eval', '-e', default=1e-10,
    help="e-value cut-off for similarity")
@click.option('--to_stop', is_flag=True, 
    help="don't translate through STOP codons")
@click.option('--cds', is_flag=True,
    help="enforce proper CDS sequences")
def dmd(**kwargs):
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
    """
    _dmd(**kwargs)

def _dmd(sequences, outdir, tmpdir, inflation, eval, to_stop, cds):
    from wgd.core import SequenceData
    s = [SequenceData(s, out_path=outdir, tmp_path=tmpdir,
        to_stop=to_stop, cds=cds) for s in sequences]
    if len(s) == 0:
        logging.error("No sequences provided!")
        return
    logging.info("tmp_dir = {}".format(s[0].tmp_path))
    if len(s) == 1:
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
    return s


# Ks distribution construction
@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('families', type=click.Path(exists=True))
@click.argument('sequences', nargs=-1, type=click.Path(exists=True))
@click.option('--outdir', '-o', default="wgd_ksd", show_default=True,
    help='tmp directory')
@click.option('--tmpdir', '-t', default=None, show_default=True,
    help='tmp directory')
@click.option('--to_stop', is_flag=True, 
    help="don't translate through STOP codons")
@click.option('--cds', is_flag=True,
    help="enforce proper CDS sequences")
@click.option('--pairwise', is_flag=True,
    help="run codeml on all gene pairs separately")
@click.option('--strip_gaps', is_flag=True,
    help="remove all gap-containing columns in the alignment")
def ksd(**kwargs):
    """
    Paranome and one-to-one ortholog Ks distribution inference pipeline.

    Example 1 - whole-paranome:

        wgd ksd families.mcl cds.fasta

    Example 2 - one-to-one orthologs (RBH):

        wgd ksd orthologs.rbh cds1.fasta cds2.fasta
    """
    _ksd(**kwargs)

def _ksd(families, sequences, outdir, tmpdir, to_stop, cds, pairwise, strip_gaps):
    from wgd.core import get_gene_families, SequenceData, KsDistributionBuilder
    from wgd.core import read_gene_families
    from wgd.viz import default_plot, apply_filters
    s = [SequenceData(s, tmp_path=tmpdir, out_path=outdir, to_stop=to_stop, cds=cds)
         for s in sequences]
    logging.info("tmpdir = {}".format(s[0].tmp_path))
    fams = read_gene_families(families)
    fams = get_gene_families(s, fams, 
            pairwise=pairwise, 
            strip_gaps=strip_gaps)
    ksdb = KsDistributionBuilder(fams)
    ksdb.get_distribution()
    prefix = os.path.basename(families)
    outfile = os.path.join(outdir, "{}.ks.tsv".format(prefix))
    logging.info("Saving to {}".format(outfile))
    ksdb.df.to_csv(outfile)
    logging.info("Making plots")
    df = apply_filters(ksdb.df, [("dS", 1e-4, 5.), ("S", 10, 1e6)])
    fig = default_plot(df, title=prefix, rwidth=0.8, bins=50)
    fig.savefig(os.path.join(outdir, "{}.ksd.svg".format(prefix)))
    fig.savefig(os.path.join(outdir, "{}.ksd.pdf".format(prefix)))


if __name__ == '__main__':
    cli()
