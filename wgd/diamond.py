import uuid
from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError


def _write_fasta(fname, seq_dict):
    with open(fname, "w") as f:
        for k, v in seq_dict.items():
            f.write(">{}\n{}\n".format(k, v))
    return fname

def _mk_tmpdir(dirname):
    if os.path.isdir(dirname):
        logging.warning("tmp dir {} exists!".format(dirname))
    else:
        os.mkdir(dirname)
    return dirname


class SequenceData:
    """
    Sequence data container for Ks distribution computation pipeline. A helper
    class that bundles sequence manipulation methods.
    """
    def __init__(self, cds_fasta, tmp_path="wgd_tmp"):
        self.tmp_path  = _mk_tmpdir(tmp_path)
        self.cds_fasta = cds_fasta
        self.pro_fasta = os.path.join(tmp_path, os.path.basename(cds_fasta))
        self.gids = {}
        self.cds_seqs = {}
        self.pro_seqs = {}
        self.read_cds()
        _write_fasta(self.pro_fasta, self.pro_seqs)

    def read_cds(self):
        prefix = os.path.basename(self.cds_fasta)
        for i, seq in enumerate(SeqIO.parse(self.cds_fasta, 'fasta')):
            gid = "{0}_{1:0>5}".format(prefix, i)
            try:
                aa_seq = seq.seq.translate(to_stop=True, cds=True)
            except TranslationError as e:
                logging.error("Translation error ({}) in seq {}".format(
                    e, seq.id))
                continue
            self.gids[seq.id] = gid
            self.cds_seqs[gid] = seq.seq
            self.pro_seqs[gid] = aa_seq
        return

    def make_diamond_db(self):
        dbname = os.path.join(self.tmp_path, str(uuid.uuid1()))
        cmd = ["diamond", "makedb", "--in", self.pro_fasta, "-d", dbname]
        completed = subprocess.run(cmd, capture_output=True)
        return dbname


def run_diamond(seqs1, seqs2):
    db = seqs1.make_diamond_db()
    cmd = ["diamond", ... ]
    completed = subprocess.run(cmd, capture_output=True)
    return


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
