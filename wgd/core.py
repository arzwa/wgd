# Arthur Zwaenepoel (2020)
import uuid
import os
import logging
import numpy as np
import pandas as pd
import subprocess as sp
import fastcluster
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Align import MultipleSeqAlignment
from Bio.Data.CodonTable import TranslationError
from Bio import Phylo
from joblib import Parallel, delayed
from wgd.codeml import Codeml
from wgd.cluster import cluster_ks


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

def _strip_gaps(aln):
    new_aln = aln[:,0:0]
    for j in range(aln.get_alignment_length()):
        if any([x == "-" for x in aln[:,j]]):
            continue
        else:
            new_aln += aln[:,j:j+1]
    return new_aln

def _pal2nal(pro_aln, cds_seqs):
    aln = {}
    for i, s in enumerate(pro_aln):
        cds_aln = ""
        cds_seq = cds_seqs[s.id].seq
        k = 0
        for j in range(pro_aln.get_alignment_length()):
            if pro_aln[i, j] == "-":
                cds_aln += "---"
            elif pro_aln[i, j] == "X":
                cds_aln += "???"  # not sure what best choice for codeml is
                k += 3
            else:
                cds_aln += cds_seq[k:k+3]
                k += 3
        aln[s.id] = cds_aln
    return MultipleSeqAlignment([SeqRecord(v, id=k) for k, v in aln.items()])

def _log_process(o, program=""):
    logging.debug("{} stderr: {}".format(program.upper(), o.stderr.decode()))
    logging.debug("{} stdout: {}".format(program.upper(), o.stdout.decode()))

def _label_internals(tree):
    for i, c in enumerate(tree.get_nonterminals()):
        c.name = str(i)

def _label_families(df):
    df.index = ["GF{:0>5}".format(i+1) for i in range(len(df.index))]


class SequenceData:
    """
    Sequence data container for Ks distribution computation pipeline. A helper
    class that bundles sequence manipulation methods.
    """
    def __init__(self, cds_fasta,
            tmp_path=None, out_path="wgd_dmd",
            to_stop=True, cds=True):
        if tmp_path == None:
            tmp_path = str(uuid.uuid4())
        self.tmp_path  = _mkdir(tmp_path)
        self.out_path  = _mkdir(out_path)
        self.cds_fasta = cds_fasta
        self.prefix    = os.path.basename(self.cds_fasta)
        self.pro_fasta = os.path.join(tmp_path, self.prefix + ".tfa")
        self.pro_db    = os.path.join(tmp_path, self.prefix + ".db")
        self.cds_seqs  = {}
        self.pro_seqs  = {}
        self.dmd_hits  = {}
        self.rbh       = {}
        self.mcl       = {}
        self.idmap     = {}  # map from the new safe id to the input seq id
        self.read_cds(to_stop=to_stop, cds=cds)
        _write_fasta(self.pro_fasta, self.pro_seqs)

    def read_cds(self, to_stop=True, cds=True):
        """
        Read a CDS fasta file. We give each input record a unique safe ID, and
        keep the full records in a dict with these IDs. We use the newly assigned
        IDs in further analyses, but can reconvert at any time.
        """
        for i, seq in enumerate(SeqIO.parse(self.cds_fasta, 'fasta')):
            gid = "{0}_{1:0>5}".format(self.prefix, i)
            try:
                aa_seq = seq.translate(to_stop=to_stop, cds=cds, id=seq.id,
                                       stop_symbol="")
            except TranslationError as e:
                logging.warning("Translation error ({}) in seq {}".format(
                    e, seq.id))
                continue
            self.cds_seqs[gid] = seq
            self.pro_seqs[gid] = aa_seq
            self.idmap[seq.id] = gid
        return

    def merge(self, other):
        """
        Merge other into self, keeping the paths etc. of self.
        """
        self.cds_seqs.update(other.cds_seqs)
        self.pro_seqs.update(other.pro_seqs)
        self.idmap.update(other.idmap)

    def make_diamond_db(self):
        cmd = ["diamond", "makedb", "--in",
               self.pro_fasta, "-d", self.pro_db]
        out = sp.run(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        logging.debug(out.stderr.decode())
        if out.returncode == 1:
            logging.error(out.stderr.decode())

    def run_diamond(self, seqs, eval=1e-10):
        self.make_diamond_db()
        run = "_".join([self.prefix, seqs.prefix + ".tsv"])
        outfile = os.path.join(self.tmp_path, run)
        cmd = ["diamond", "blastp", "-d", self.pro_db, "-q",
            seqs.pro_fasta, "-o", outfile]
        out = sp.run(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
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

    def write_paranome(self, fname=None):
        if not fname:
            fname = os.path.join(self.out_path, "{}.tsv".format(self.prefix))
        with open(fname, "w") as f:
            f.write("\t" + self.prefix + "\n")
            for i, (k, v) in enumerate(sorted(self.mcl.items())):
                # We report original gene IDs
                f.write("GF{:0>5}\t".format(i+1))
                f.write(", ".join([self.cds_seqs[x].id for x in v]))
                f.write("\n")
        return fname

    def write_rbh_orthologs(self, seqs):
        prefix = seqs.prefix
        fname = "{}_{}.rbh".format(self.prefix, prefix)
        fname = os.path.join(self.out_path, fname)
        df = self.rbh[prefix]
        df[prefix] = df[0].apply(lambda x: seqs.cds_seqs[x].id)
        df[self.prefix] = df[1].apply(lambda x: self.cds_seqs[x].id)
        _label_families(df)
        df.to_csv(fname, columns=[prefix, self.prefix], sep="\t")

    def remove_tmp(self, prompt=True):
        if prompt:
            ok = input("Removing {}, sure? [y|n]".format(self.tmp_path))
            if ok != "y":
                return
        out = sp.run(["rm", "-r", self.tmp_path], stdout=sp.PIPE, stderr=sp.PIPE)
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
        out = sp.run(command, stdout=sp.PIPE, stderr=sp.PIPE)
        _log_process(out)
        command = ['mcl', g3, '-I', str(inflation), '-o', g4]
        logging.debug(" ".join(command))
        out = sp.run(command, stdout=sp.PIPE, stderr=sp.PIPE)
        _log_process(out)
        command = ['mcxdump', '-icl', g4, '-tabr', g2, '-o', outfile]
        _log_process(out)
        out = sp.run(command, stdout=sp.PIPE, stderr=sp.PIPE)
        _log_process(out)
        return outfile


# Gene family i/o
def _rename(family, ids):
    return [ids[x] for x in family]

def read_gene_families(fname):
    return pd.read_csv(fname, sep="\t", index_col=0).apply(
            lambda y: [x.split(", ") for x in y if x != ""])

def merge_seqs(seqs):
    if type(seqs) == list:
        if len(seqs) > 2:
            raise ValueError("More than two sequence data objects?")
        if len(seqs) == 2:
            seqs[0].merge(seqs[1])
        seqs = seqs[0]
    return seqs

def get_gene_families(seqs, families, rename=True, **kwargs):
    """
    Get the `GeneFamily` objects from a list of families (list with lists of
    gene IDs) and sequence data. When `rename` is set to True, it is assumed
    the gene IDs in the families are the original IDs (not those assigned 
    when reading the CDS from file).

    Note: currently this is only defined for one or two genomes (paranome 
    and one-to-one orthologs), but it should easily generalize to arbitrary
    gene families.
    """
    seqs = merge_seqs(seqs)
    gene_families = []
    for fid in families.index:
        family = []
        for col in families.columns:
            ids = families.loc[fid][col]
            if rename:
                family += _rename(ids, seqs.idmap)
            else:
                family += ids
        if len(family) > 1:
            cds = {x: seqs.cds_seqs[x] for x in family}
            pro = {x: seqs.pro_seqs[x] for x in family}
            tmp = os.path.join(seqs.tmp_path, fid)
            gene_families.append(GeneFamily(fid, cds, pro, tmp, **kwargs))
        else:
            logging.debug("Skipping singleton family {}{}".format(fid,family))
    return gene_families


# NOTE: It would be nice to implement an option to do a complete approach
# where we use the tree in codeml to estimate Ks-scale branch lengths?
class GeneFamily:
    def __init__(self, gfid, cds, pro, tmp_path,
            aligner="mafft", tree_method="cluster", ks_method="GY94",
            eq_freq="F3X4", kappa=None, prequal=False, strip_gaps=False,
            min_length=3, codeml_iter=1, aln_options="--auto", 
            tree_options="-m LG", pairwise=False):
        self.id = gfid
        self.cds_seqs = cds
        self.pro_seqs = pro
        self.tmp_path = _mkdir(tmp_path)
        self.cds_fasta = os.path.join(self.tmp_path, "cds.fasta")
        self.pro_fasta = os.path.join(self.tmp_path, "pro.fasta")
        self.cds_alnf = os.path.join(self.tmp_path, "cds.aln")
        self.pro_alnf = os.path.join(self.tmp_path, "pro.aln")
        self.cds_aln = None
        self.pro_aln = None
        self.codeml_results = None
        self.no_codeml_results = None
        self.tree = None
        self.out = os.path.join(self.tmp_path, "{}.csv".format(gfid))

        # config
        self.aligner = aligner  # mafft | prank | muscle
        self.tree_method = tree_method  # iqtree | fasttree | alc
        self.ks_method = ks_method  # GY | NG
        self.kappa = kappa
        self.eq_freq = eq_freq
        self.prequal = prequal
        self.strip_gaps = strip_gaps  # strip gaps based on overall alignment
        self.codeml_iter = codeml_iter
        self.min_length = min_length  # minimum length of codon alignment
        self.aln_options = aln_options
        self.tree_options = tree_options
        self.pairwise = pairwise

    def get_ks(self):
        logging.info("Analysing family {}".format(self.id))
        self.align()
        self.run_codeml()
        if self.codeml_results is not None:
            self.get_tree()
            self.compile_dataframe()
        self.combine_results()

    def combine_results(self):
        if self.no_codeml_results is None:
            return
        self.codeml_results = pd.concat([self.codeml_results, self.no_codeml_results])
    
    def nan_result(self, pairs):
        """
        For a bunch of pairs obtain a data frame with missing data.
        """
        if len(pairs) == 0: return None
        data = []
        for pair in pairs:
            pairid = "__".join(sorted(pair))
            data.append({
                "pair": pairid, 
                "gene1": pair[0], 
                "gene2": pair[1], 
                "family": self.id})
        return pd.DataFrame.from_records(data).set_index("pair")

    # NOT TESTED
    def run_prequal(self):
        cmd = ["prequal", self.pro_fasta]
        out = sp.run(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        _log_process(out, program="prequal")
        self.pro_fasta = "{}.filtered".format(self.pro_fasta)

    def align(self):
        _write_fasta(self.pro_fasta, self.pro_seqs)
        if self.prequal:
            self.run_prequal()
        if self.aligner == "mafft":
            self.run_mafft(options=self.aln_options)
        else:
            logging.error("Unsupported aligner {}".format(self.aligner))
        self.get_codon_alignment()

    def run_mafft(self, options="--auto"):
        cmd = ["mafft"] + options.split() + ["--amino", self.pro_fasta]
        out = sp.run(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        with open(self.pro_alnf, 'w') as f: f.write(out.stdout.decode('utf-8'))
        _log_process(out, program="mafft")
        self.pro_aln = AlignIO.read(self.pro_alnf, "fasta")

    def get_codon_alignment(self):
        self.cds_aln = _pal2nal(self.pro_aln, self.cds_seqs)
        if self.strip_gaps:
            self.cds_aln = _strip_gaps(self.cds_aln)

    def run_codeml(self):
        codeml = Codeml(self.cds_aln, exe="codeml", tmp=self.tmp_path, prefix=self.id)
        # TODO, do something with `no_result`
        if self.pairwise:
            result, no_result = codeml.run_codeml_pairwise(
                    preserve=True, times=self.codeml_iter)
        else:
            result, no_result = codeml.run_codeml(
                    preserve=True, times=self.codeml_iter)
        self.codeml_results = result
        self.no_codeml_results = self.nan_result(no_result)

    def get_tree(self):
        # dispatch method
        # This likely will have to catch families of only two or three members.
        if self.tree_method == "cluster":
            tree = self.cluster()
        elif self.tree_method == "iqtree":
            tree = self.run_iqtree(options=self.tree_options)
        elif self.tree_method == "fasttree":
            tree = self.run_fasttree()
        self.tree = tree

    def run_iqtree(self, options="-m LG"):
        cmd = ["iqtree", "-s", self.pro_alnf] + options.split()
        out = sp.run(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        _log_process(out, program="iqtree")
        tree = Phylo.read(self.pro_alnf + ".treefile", format="newick")
        _label_internals(tree)
        return tree

    def run_fasttree(self):
        tree_pth = self.pro_alnf + ".nw"
        cmd = ["fasttree", '-out', tree_pth, self.pro_alnf]
        out = sp.run(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        _log_process(out, program="fasttree")
        tree = Phylo.read(self.pro_alnf + ".nw", format="newick")
        _label_internals(tree)
        return tree 

    def cluster(self):
        return cluster_ks(self.codeml_results)

    def compile_dataframe(self):
        n = len(self.cds_seqs)
        d = {}
        l = self.tree.get_terminals()
        for i in range(len(l)):
            gi = l[i].name
            for j in range(i+1, len(l)):
                gj = l[j].name
                pair = "__".join(sorted([gi, gj]))
                node = self.tree.common_ancestor(l[i], l[j])
                length = self.cds_aln.get_alignment_length()
                d[pair] = {"node": node.name, "family": self.id, "alnlen": length}
        df = pd.DataFrame.from_dict(d, orient="index")
        self.codeml_results = self.codeml_results.join(df)

def _get_ks(family):
    family.get_ks()
    family.codeml_results.to_csv(family.out)


class KsDistributionBuilder:
    def __init__(self, gene_families, n_threads=4):
        self.families = gene_families
        self.df = None
        self.n_threads = n_threads

    def get_distribution(self):
        Parallel(n_jobs=self.n_threads)(
            delayed(_get_ks)(family) for family in self.families)
        self.df = pd.concat([pd.read_csv(x.out, index_col=0) 
            for x in self.families], sort=True)
