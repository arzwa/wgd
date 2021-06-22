import os
import logging
import tempfile
import pandas as pd
import subprocess as sp
from collections import defaultdict
from operator import itemgetter

gff_header = ["gene", "scaffold", "start", "orientation"] 

# we construct first a table with each row a gene containing it's family,
# scaffold, location and orientation
def make_gene_table(gffs, families, feature, attribute):
    """
    Construct a table from a bunch of gff files and gene families such that all
    information for synteny and co-linearity related analyses is in one place.
    
    :param gffs: a list of GFF *file names*
    :param families: a pandas data frame with the gene families
    :param feature: feature tag for the GFF files
    :param attribute: attribute tag for parsing out gene IDs from last column
    
    Note: currently we assume the same attribute/feature is required for all
    GFF files, and the burden is on the user to make sure this is the case. We
    may wish, for flexibility, to allow a list of feature/attribute arguments
    in the future, one for each GFF file.
    """
    gfftables = [gff2table(gff, feature, attribute) for gff in gffs]
    familytable = gene2family(families)
    df = pd.concat(gfftables)
    df = familytable.join(df)
    return df

def getattr(s, attribute):
    for x in s.split(";"):
        y = x.split("=")
        if y[0] == attribute:
            return y[1].strip()
    return ""

def gff2table(gff, feature, attribute):
    """
    Read a GFF file to a pandas data frame, from a filename.
    """
    rows = []
    with open(gff, "r") as f:
        for l in f.readlines():
            if l.startswith("#"):
                continue
            x = l.split("\t")
            if x[2] == feature:
                a = getattr(x[-1], attribute)
                rows.append({"gene": a, "scaffold": x[0], "start": int(x[3]), "or": x[6]})
    df = pd.DataFrame.from_dict(rows).set_index("gene")
    return df

def gene2family(families):
    """
    Get a gene to family ID mapping from a gene families data frame in
    OrthoFinder format.
    """
    rows = []
    for fam in families.index:
        for sp in families.columns:
            x = families.loc[fam, sp]
            if type(x) != str:
                continue
            for gene in x.split(", "):
                rows.append({"gene": gene.strip(), 
                    "species": sp, "family": fam})
    df = pd.DataFrame.from_dict(rows).set_index("gene")
    return df

def configure_adhore(table, outdir, **kwargs):
    """
    Write out all required files for running I-ADHoRe.
    """
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    famfile = os.path.join(outdir, "families.tsv")
    logging.info("Writing families file")
    write_families(famfile, table)
    logging.info("Writing gene lists")
    genelists = write_genelists(outdir, table)
    logging.info("Writing config file")
    config_file = os.path.join(outdir, "iadhore.conf")
    out_path = os.path.abspath(os.path.join(outdir, "iadhore-out"))
    write_config_adhore(config_file, out_path, genelists, os.path.abspath(famfile), **kwargs)
    return config_file, out_path

def write_families(fname, table):
    table["family"].to_csv(fname, sep="\t", header=None)

def write_genelists(outdir, table):
    genomes = {}
    for sp, df in table.groupby("species"):
        gdir = os.path.join(outdir, sp + ".lists")
        if not os.path.exists(gdir):
            os.mkdir(gdir)
        lists = {}
        for scaffold, sdf in df.groupby("scaffold"):
            if len(sdf.index) <= 2: continue
            fname = os.path.join(gdir, scaffold)
            with open(fname, "w") as o:
                for g in sdf.sort_values(by=["start"]).index:
                    o.write(g + sdf.loc[g,"or"] + "\n")
            lists[scaffold] = os.path.abspath(fname)
        genomes[sp] = lists
    return genomes

def write_config_adhore(
        config_file, output_path, genelists, families, gap_size=30,
        cluster_gap=35, q_value=0.75, prob_cutoff=0.01, anchor_points=3,
        alignment_method='gg2', level_2_only='false', table_type='family',
        multiple_hypothesis_correction='FDR', visualize_ghm='false',
        visualize_alignment='true'):
    """
    Write out the config file for I-ADHoRe. See I-ADHoRe manual for information
    on parameter settings.

    :param gene_lists: directory with gene lists per chromosome
    :param families: file with gene to family mapping
    :param config_file_name: name for the config file
    :param output_path: output path name
    :param gap_size: see I-ADHoRe 3.0 documentation
    :param cluster_gap: see I-ADHoRe 3.0 documentation
    :param q_value: see I-ADHoRe 3.0 documentation
    :param prob_cutoff: see I-ADHoRe 3.0 documentation
    :param anchor_points: see I-ADHoRe 3.0 documentation
    :param alignment_method: see I-ADHoRe 3.0 documentation
    :param level_2_only: see I-ADHoRe 3.0 documentation
    :param table_type: see I-ADHoRe 3.0 documentation
    :param multiple_hypothesis_correction: see I-ADHoRe 3.0 documentation
    :param visualize_ghm: see I-ADHoRe 3.0 documentation
    :param visualize_alignment: see I-ADHoRe 3.0 documentation
    :return: configuration file see I-ADHoRe 3.0 documentation
    """
    with open(config_file, 'w') as o:
        for k,v in genelists.items():
            o.write('genome= {}\n'.format(k))
            for scaffold, fname in v.items():
                o.write("{} {}\n".format(scaffold, fname))
            o.write("\n")

        o.write('blast_table= {}\n'.format(families))
        o.write('output_path= {}\n'.format(output_path))
        o.write('gap_size= {}\n'.format(gap_size))
        o.write('q_value= {}\n'.format(q_value))
        o.write('cluster_gap= {}\n'.format(cluster_gap))
        o.write('prob_cutoff= {}\n'.format(prob_cutoff))
        o.write('anchor_points= {}\n'.format(anchor_points))
        o.write('alignment_method= {}\n'.format(alignment_method))
        o.write('level_2_only= {}\n'.format(level_2_only))
        o.write('table_type= {}\n'.format(table_type))
        o.write('multiple_hypothesis_correction= {}\n'.format(
                multiple_hypothesis_correction))
        o.write('visualizeGHM= {}\n'.format(visualize_ghm))
        o.write('visualizeAlignment= {}\n'.format(visualize_alignment))
    return os.path.abspath(output_path)

def run_adhore(config_file):
    """
    Run I-ADHoRe for a given config file

    :param config_file: path to I-ADHoRe configuration file
    """
    completed = sp.run(['i-adhore', config_file], stderr=sp.PIPE, stdout=sp.PIPE)
    logging.warning(completed.stderr.decode('utf-8'))
    logging.info(completed.stdout.decode('utf-8'))
    return

def get_anchors(out_path):
    anchors = pd.read_csv(os.path.join(out_path, "anchorpoints.txt"), sep="\t", index_col=0)
    anchors["pair"] = anchors[["gene_x", "gene_y"]].apply(lambda x: "__".join(sorted([x[0], x[1]])), axis=1)
    # there are duplicates, due to anchors being in multiple multiplicons
    return anchors[["pair", "multiplicon"]].drop_duplicates("pair").set_index("pair")

def get_anchor_ksd(ks_distribution, anchors):
    return ks_distribution.join(anchors).dropna()

def get_segments_profile(out_path):
    segs = pd.read_csv(os.path.join(out_path, "segments.txt"), sep="\t", index_col=0)
    le = pd.read_csv(os.path.join(out_path, "list_elements.txt"), sep="\t", index_col=0)
    segs = segs.join(le.set_index("segment"), how="inner")
    segs["segment"] = segs.index
    counted = segs.groupby(["multiplicon", "genome"])["segment"].aggregate(lambda x: len(set(x)))
    profile = counted.unstack(level=-1).fillna(0)
    return profile



