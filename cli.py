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
@click.option('--nthreads', '-n', default=4, show_default=True,
    help="number of threads to use")
@click.option('--to_stop', is_flag=True, 
    help="don't translate through STOP codons")
@click.option('--cds', is_flag=True,
    help="enforce proper CDS sequences")
@click.option('--pairwise', is_flag=True,
    help="run codeml on all gene pairs separately")
@click.option('--strip_gaps', is_flag=True,
    help="remove all gap-containing columns in the alignment")
@click.option('--tree_method', '-tree', 
    type=click.Choice(['cluster', 'fasttree', 'iqtree']), 
    default='cluster', show_default=True,
    help="Tree inference method for node weighting")
def ksd(**kwargs):
    """
    Paranome and one-to-one ortholog Ks distribution inference pipeline.

    Example 1 - whole-paranome:

        wgd ksd families.mcl cds.fasta

    Example 2 - one-to-one orthologs (RBH):

        wgd ksd orthologs.rbh cds1.fasta cds2.fasta
    """
    _ksd(**kwargs)

def _ksd(families, sequences, outdir, tmpdir, nthreads, to_stop, cds, pairwise,
        strip_gaps, tree_method):
    from wgd.core import get_gene_families, SequenceData, KsDistributionBuilder
    from wgd.core import read_gene_families, merge_seqs
    from wgd.viz import default_plot, apply_filters
    if len(sequences) == 2:
        tree_method = "cluster"  # for RBH others don't make sense (and crash)
    s = merge_seqs([SequenceData(s, tmp_path=tmpdir, out_path=outdir,
            to_stop=to_stop, cds=cds) for s in sequences])
    logging.info("tmpdir = {}".format(s.tmp_path))
    fams = read_gene_families(families)
    fams = get_gene_families(s, fams, 
            pairwise=pairwise, 
            strip_gaps=strip_gaps,
            tree_method=tree_method)
    ksdb = KsDistributionBuilder(fams, s, n_threads=nthreads)
    ksdb.get_distribution()
    prefix = os.path.basename(families)
    outfile = os.path.join(outdir, "{}.ks.tsv".format(prefix))
    logging.info("Saving to {}".format(outfile))
    ksdb.df.to_csv(outfile)
    logging.info("Making plots")
    df = apply_filters(ksdb.df, [("dS", 1e-4, 5.), ("S", 10, 1e6)])
    ylabel = "Duplications"
    if len(sequences) == 2:
        ylabel = "RBH orthologs"
    fig = default_plot(df, title=prefix, rwidth=0.8, bins=50, ylabel=ylabel)
    fig.savefig(os.path.join(outdir, "{}.ksd.svg".format(prefix)))
    fig.savefig(os.path.join(outdir, "{}.ksd.pdf".format(prefix)))
    

@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('gff_files', nargs=-1, type=click.Path(exists=True))
@click.argument('families', type=click.Path(exists=True))
@click.option('--ks_distribution', '-ks', default=None,
    help="ks distribution tsv file (optional, see `wgd ksd`)")
@click.option('--outdir', '-o', default='./wgd_syn', show_default=True, 
    help='output directory')
@click.option('--feature', '-f', default='gene', show_default=True,
    help="keyword for parsing the genes from the GFF file (column 3)")
@click.option('--attribute', '-a', default='ID', show_default=True,
    help="keyword for parsing the gene IDs from the GFF file (column 9)")
@click.option('--min_length', '-l', default=250, show_default=True,
    help="minimum length of a genomic element to be included in dotplot.")
@click.option('--ks_range', '-r', nargs=2, default=(0.05, 5), show_default=True,
    type=float, help='Ks range to use for colored dotplot')
@click.option('--iadhore_options', default="",
    help="other options for I-ADHoRe, as a comma separated string, "
         "e.g. gap_size=30,q_value=0.75,prob_cutoff=0.05")
def syn(**kwargs):
    _syn(**kwargs)

def _syn(gff_files, families, ks_distribution, outdir, feature, attribute,
        min_length, ks_range, iadhore_options):
    """
    Co-linearity and anchor inference using I-ADHoRe.
    """
    from wgd.syn import make_gene_table, configure_adhore, run_adhore
    from wgd.syn import get_anchors, get_anchor_ksd, get_segments_profile
    from wgd.viz import default_plot, apply_filters, syntenic_depth_plot
    # non-default options for I-ADHoRe
    iadhore_opts = {x.split("=")[0].strip(): x.split("=")[1].strip()
               for x in iadhore_options.split(",") if x != ""}
    if len(iadhore_opts) > 0:
        logging.info("I-ADHoRe 3.0 options: {}".format(iadhore_opts))
    # read families and make table
    prefix = os.path.basename(families)
    families = pd.read_csv(families, index_col=0, sep="\t")
    table = make_gene_table(gff_files, families, feature, attribute)
    logging.info("Configuring I-ADHoRe co-linearity search")
    conf, out_path = configure_adhore(table, outdir, **iadhore_opts)
    table.to_csv(os.path.join(outdir, "gene-table.csv"))
    logging.info("Running I-ADHoRe")
    run_adhore(conf)

    # general post-processing
    logging.info("Processing I-ADHoRe output")
    anchors = get_anchors(out_path)
    anchors.to_csv(os.path.join(outdir, "anchors.csv"))
    segprofile = get_segments_profile(out_path)
    segprofile.to_csv(os.path.join(outdir, "segprofile.csv"))
    fig = syntenic_depth_plot(segprofile)
    fig.savefig(os.path.join(outdir, "{}.syndepth.svg".format(prefix)))
    fig.savefig(os.path.join(outdir, "{}.syndepth.pdf".format(prefix)))

    # anchor Ks distributions
    if ks_distribution:
        ylabel = "Duplications"
        if len(gff_files) == 2:
            ylabel = "RBH orthologs"
        ksd = pd.read_csv(ks_distribution, index_col=0)
        anchor_ks = get_anchor_ksd(ksd, anchors)
        anchor_ks.to_csv(os.path.join(outdir, "{}.anchors.ks.csv".format(prefix)))
        a = apply_filters(ksd,       [("dS", 1e-4, 5.), ("S", 10, 1e6)])
        b = apply_filters(anchor_ks, [("dS", 1e-4, 5.), ("S", 10, 1e6)])
        fig = default_plot(a, b, title=prefix, rwidth=0.8, bins=50, ylabel=ylabel)
        fig.savefig(os.path.join(outdir, "{}.ksd.svg".format(prefix)))
        fig.savefig(os.path.join(outdir, "{}.ksd.pdf".format(prefix)))


if __name__ == '__main__':
    cli()
