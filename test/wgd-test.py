import unittest
import os
import logging
import shutil
from wgd.core import SequenceData, get_gene_families, KsDistributionBuilder

# some config: set logging level, and get directory
logging.basicConfig(level=logging.ERROR)
thisdir = os.path.dirname(os.path.abspath(__file__))


# Not sure how to structure exactly...
class TestSeqData(unittest.TestCase):

    def setUp(self):
        self.datadir = os.path.join(thisdir, "data")
        self.tmpdir  = os.path.join(thisdir, "_tmpdir")
        try:
            os.mkdir(self.tmpdir)
        except FileExistsError:
            raise FileExistsError("Will not proceed with"
                " tests as tmpdir {} already exists"
                "!".format(self.tmpdir))
        self.fasta1 = os.path.join(self.datadir, "ugi1000.fasta")
        self.fasta2 = os.path.join(self.datadir, "egu1000.fasta")

    def test_seqio(self):
        logging.info("Testing `SequenceData` (IO)")
        s = SequenceData(self.fasta1,
                         out_path=self.tmpdir,
                         tmp_path=self.tmpdir,
                         cds=True, to_stop=True)
        self.assertEqual(len(s.cds_seqs), 224,
                         "Should only read proper CDS")
        s1 = SequenceData(self.fasta1,
                         out_path=self.tmpdir,
                         tmp_path=self.tmpdir,
                         cds=False, to_stop=False)
        self.assertEqual(len(s1.cds_seqs), 500,
                         "Should read full file")
        s2 = SequenceData(self.fasta2,
                         out_path=self.tmpdir,
                         tmp_path=self.tmpdir,
                         cds=False, to_stop=False)
        self.assertEqual(len(s2.pro_seqs), 500,
                         "Proper translation")

    def test_paranome(self):
        logging.info("Testing paranome (requires diamond+mcl)")
        s = SequenceData(self.fasta1,
                         out_path=self.tmpdir,
                         tmp_path=self.tmpdir,
                         cds=False, to_stop=False)
        s.get_paranome()
        self.assertEqual(len(s.mcl), 35, "Got the paranome")

        # round tripping for families?
        # write to file/read from file
        families = s.write_paranome()
        with open(families, "r") as f:
            fams = [x.strip().split("\t") for x in f.readlines()]
        fams = get_gene_families(s, fams)
        self.assertEqual(len(s.mcl), len(fams), "Round tripped")

        # construct gene families directly from `self.mcl`
        fams2 = get_gene_families(s, s.mcl.values(), rename=False)
        self.assertEqual(len(fams2), len(fams))

    def test_rbh(self):
        logging.info("Testing RBH orthologs (requires diamond)")
        s1 = SequenceData(self.fasta1,
                          out_path=self.tmpdir,
                          tmp_path=self.tmpdir,
                          cds=False, to_stop=False)
        s2 = SequenceData(self.fasta2,
                          out_path=self.tmpdir,
                          tmp_path=self.tmpdir,
                          cds=False, to_stop=False)
        s1.get_rbh_orthologs(s2, eval=1e-3)
        df = s1.rbh[s2.prefix]
        self.assertEqual(len(df.index), 19, "Got all RBHs")
        self.assertEqual(len(df.columns), 12, "Tight No. columns")

    def test_ksd(self):
        s = SequenceData(self.fasta1,
                         out_path=self.tmpdir,
                         tmp_path=self.tmpdir,
                         cds=True, to_stop=True)
        s.get_paranome()

        # alignment without gap-stripping
        fams = get_gene_families(s, s.mcl.values(), rename=False,
                                 prequal=False, strip_gaps=False)
        f = fams[0]
        f.align()
        cds_len = f.cds_aln.get_alignment_length()
        pro_len = f.pro_aln.get_alignment_length()
        self.assertEqual(cds_len, 558, "Alignment working?")
        self.assertEqual(cds_len, 3*pro_len, "AA -> codon fine?")

        # alignment with gap-stripping
        fams = get_gene_families(s, s.mcl.values(), rename=False,
                                 prequal=False, strip_gaps=True)
        f = fams[0]
        f.align()
        cds_len = f.cds_aln.get_alignment_length()
        pro_len = f.pro_aln.get_alignment_length()
        gaps = 0
        for i in range(pro_len):
            if "-" in f.pro_aln[:,i]:
                gaps += 1
        self.assertEqual(cds_len, 3*(pro_len - gaps))

        # Codeml

        # Tree inference

        # run Ks distribution construction
        #ksdb = KsDistributionBuilder(fams)
        #ksdb.get_distribution()

    def tearDown(self):
        # Cleaning up (remove tmpdir)
        shutil.rmtree(self.tmpdir)

if __name__ == '__main__':
    unittest.main()


