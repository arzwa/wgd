#!/usr/bin/python3.5
import click
import logging
import sys
import os
import warnings
import pandas as pd
import subprocess as sp
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

__version__ = "2.0.0"

# CLI entry point
@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--verbosity', '-v', type=click.Choice(['info', 'debug']),
    default='info', help="Verbosity level, default = info.")
@click.option('--logfile', '-l', default=None,
    help="File to write logs to (optional)")
@click.option('--version', is_flag=True, help="Print version number")
def cli(verbosity, logfile, version):
    """
    wgd - Copyright (C) 2019 Arthur Zwaenepoel\n
    Contact: arzwa@psb.vib-ugent.be
    """
    logging.basicConfig(filename=logfile,
        filemode='a',
        format='%(asctime)s: %(levelname)s %(message)s',
        datefmt='%H:%M:%S',
        level=verbosity.upper())
    if version:
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
@click.option('--ignorestop', is_flag=True,
    help="translate through STOP codons")
@click.option('--nostrictcds', is_flag=True,
    help="don't enforce proper CDS sequences")
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

def _dmd(sequences, outdir, tmpdir, inflation, eval, ignorestop, nostrictcds):
    from wgd.core import SequenceData
    s = [SequenceData(s, out_path=outdir, tmp_path=tmpdir,
        to_stop=not ignorestop, cds=nostrictcds) for s in sequences]
    if len(s) == 0:
        logging.error("No sequences provided!")
        return
    elif len(s) == 1:
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
@click.option('--tmpdir', '-t', default=None, show_default=True,
    help='tmp directory')
def ksd(families, sequences, tmpdir):
    from wgd.core import get_gene_families, SequenceData, KsDistributionBuilder
    s = [SequenceData(s, tmp_path=tmpdir, to_stop=False, cds=False)
         for s in sequences]
    with open(families, "r") as f:
        fams = [x.strip().split("\t") for x in f.readlines()]
    fams = get_gene_families(s, fams)
    ksdb = KsDistributionBuilder(fams)
    ksdb.get_distribution()


if __name__ == '__main__':
    cli()
