class SequenceData:
    """
    Sequence data container for Ks distribution computation pipeline. A helper
    class that bundles sequence manipulation methods.
    """
    def __init__(self, cds_fasta, tmp_path="wgd_tmp"):
        self.tmp_path = tmp_path
        self.cds_fasta = cds_fasta
        self.pro_fasta = os.path.join(tmp_path, os.path.basename(cds_fasta))
        self.cds_seqs = self.read_cds()
        self.pro_seqs = self.translate()

    def read_cds(self):
        pass

    def translate(self):
        pass


class MCL:
    def __init__(self, cds_fasta, out_path="wgd_mcl"):
        self.out_path = out_path

    def make_diamond_db(self):

        cmd = ["diamond", "makedb", "--in", self.pro_sequences, "-d",

    def run_diamond(self):
        pass

    def get_graph(self):
        pass

    def run_mcl(self):
        pass


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
