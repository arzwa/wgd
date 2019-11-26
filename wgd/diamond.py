import uuid
import os
import logging
import pandas as pd
import subprocess as sp
from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError


# helper functions
def _write_fasta(fname, seq_dict):
    with open(fname, "w") as f:
        for k, v in seq_dict.items():
            f.write(">{}\n{}\n".format(k, v.seq))
    return fname

def _mkdir(dirname):
    if os.path.isdir(dirname):
        logging.warning("dir {} exists!".format(dirname))
    else:
        os.mkdir(dirname)
    return dirname


# keep in dict with keys safe ids, with as values the full record, allowing at
# all time full recovery of gene names etc.?
class SequenceData:
    """
    Sequence data container for Ks distribution computation pipeline. A helper
    class that bundles sequence manipulation methods.
    """
    def __init__(self, cds_fasta, tmp_path=None, out_path="wgd_dmd",
            to_stop=True, cds=True):
        if tmp_path == None:
            tmp_path = str(uuid.uuid4())
        self.tmp_path  = _mkdir(tmp_path)
        self.out_path  = _mkdir(out_path)
        self.cds_fasta = cds_fasta
        self.prefix = os.path.basename(self.cds_fasta)
        self.pro_fasta = os.path.join(tmp_path, self.prefix + ".tfa")
        self.pro_db = os.path.join(tmp_path, self.prefix + ".db")
        self.cds_seqs = {}
        self.pro_seqs = {}
        self.dmd_hits = {}
        self.rbh = {}
        self.mcl = {}
        self.read_cds(to_stop=to_stop, cds=cds)
        _write_fasta(self.pro_fasta, self.pro_seqs)

    def read_cds(self, to_stop=True, cds=True):
        for i, seq in enumerate(SeqIO.parse(self.cds_fasta, 'fasta')):
            gid = "{0}_{1:0>5}".format(self.prefix, i)
            try:
                aa_seq = seq.translate(to_stop=to_stop, cds=cds, id=seq.id)
            except TranslationError as e:
                logging.error("Translation error ({}) in seq {}".format(
                    e, seq.id))
                continue
            self.cds_seqs[gid] = seq
            self.pro_seqs[gid] = aa_seq
        return

    def make_diamond_db(self):
        cmd = ["diamond", "makedb", "--in", self.pro_fasta, "-d", self.pro_db]
        out = sp.run(cmd, capture_output=True)
        logging.debug(out.stderr.decode())
        if out.returncode == 1:
            logging.error(out.stderr.decode())

    def run_diamond(self, seqs, eval=1e-10):
        self.make_diamond_db()
        run = "_".join([self.prefix, seqs.prefix + ".tsv"])
        outfile = os.path.join(self.tmp_path, run)
        cmd = ["diamond", "blastp", "-d", self.pro_db, "-q",
            seqs.pro_fasta, "-o", outfile]
        out = sp.run(cmd, capture_output=True)
        logging.debug(out.stderr.decode())
        df = pd.read_csv(outfile, sep="\t", header=None)
        df = df.loc[df[0] != df[1]]
        self.dmd_hits[seqs.prefix] = df = df.loc[df[10] <= eval]
        return df

    def get_rbh_orthologs(self, seqs, eval=1e-10):
        if self == seqs:
            raise ValueError("RBH orthologs only defined for distinct species")
        df = self.run_diamond(seqs, eval=eval)
        df1 = df.sort_values(10).drop_duplicates([0])
        df2 = df.sort_values(10).drop_duplicates([1])
        self.rbh[seqs.prefix] = df1.merge(df2)
        # self.rbh[seqs.prefix] = seqs.rbh[self.prefix] = df1.merge(df2)
        # write to file using original ids for next steps

    def get_paranome(self, inflation=1.5, eval=1e-10):
        df = self.run_diamond(self, eval=eval)
        gf = self.get_mcl_graph(self.prefix)
        mcl_out = gf.run_mcl(inflation=inflation)
        with open(mcl_out, "r") as f:
            for i, line in enumerate(f.readlines()):
                self.mcl[i] = line.strip().split()

    def get_mcl_graph(self, *args):
        # args are keys in `self.dmd_hits` to use for building MCL graph
        gf = os.path.join(self.tmp_path, "_".join([self.prefix] + list(args)))
        df = pd.concat([self.dmd_hits[x] for x in args])
        df.to_csv(gf, sep="\t", header=False, index=False, columns=[0,1,10])
        return SequenceSimilarityGraph(gf)

    def write_paranome(self):
        fname = os.path.join(self.out_path, "{}.mcl".format(self.prefix))
        with open(fname, "w") as f:
            for k, v in sorted(self.mcl.items()):
                f.write("\t".join([self.cds_seqs[x].id for x in v]))
                f.write("\n")
        return fname

    def write_rbh_orthologs(self, seqs):
        prefix = seqs.prefix
        fname = "{}_{}.rbh".format(self.prefix, prefix)
        fname = os.path.join(self.out_path, fname)
        df = self.rbh[prefix]
        df["x"] = df[0].apply(lambda x: seqs.cds_seqs[x].id)
        df["y"] = df[1].apply(lambda x: self.cds_seqs[x].id)
        df.to_csv(fname, columns=["x", "y"], header=None, index=False, sep="\t")
        # header=[prefix, self.prefix]

    def remove_tmp(self, prompt=True):
        if prompt:
            ok = input("Removing {}, sure? [y|n]".format(self.tmp_path))
            if ok != "y":
                return
        out = sp.run(["rm", "-r", self.tmp_path], capture_output=True)
        logging.debug(out.stderr.decode())


class SequenceSimilarityGraph:
    def __init__(self, graph_file):
        self.graph_file = graph_file

    def run_mcl(self, inflation=1.5):
        g1 = self.graph_file
        g2 = g1 + ".tab"
        g3 = g1 + ".mci"
        g4 = g2 + ".I{}".format(inflation*10)
        outfile = g1 + ".mcl"
        command = ['mcxload', '-abc', g1, '--stream-mirror',
            '--stream-neg-log10', '-o', g3, '-write-tab', g2]
        logging.debug(" ".join(command))
        out = sp.run(command, capture_output=True)
        logging.debug(out.stderr.decode())
        command = ['mcl', g3, '-I', str(inflation), '-o', g4]
        logging.debug(" ".join(command))
        out = sp.run(command, capture_output=True)
        logging.debug(out.stderr.decode())
        command = ['mcxdump', '-icl', g4, '-tabr', g2, '-o', outfile]
        logging.debug(" ".join(command))
        out = sp.run(command, capture_output=True)
        logging.debug(out.stderr.decode())
        return outfile





class KsDistribution:
    def __init__(self):
        # self.tmp_path = ...
        # self.out_path = ...
        # self.gene_families = ...
        pass


class Codeml:
    pass


class CodonAlignment:
    def __init__(self):
        # self.pro_sequences = ...
        # self.cds_sequences = ...
        # self.pro_alignment = ...
        # self.cds_alignment = ...
        # self.aligner = ...
        # self.prequal = ...
        pass
