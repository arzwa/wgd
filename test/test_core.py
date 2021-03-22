import pytest
import os
import logging
import shutil
from wgd.core import SequenceData, get_gene_families, read_gene_families 
from wgd.core import KsDistributionBuilder

# some config: set logging level, and get directory
logging.basicConfig(level=logging.ERROR)
thisdir = os.path.dirname(os.path.abspath(__file__))


# Test the internals
class TestCore:
    @pytest.fixture
    def data(self):
        datadir = os.path.join(thisdir, "data")
        s1 = os.path.join(datadir, "ugi1000.fasta")
        s2 = os.path.join(datadir, "egu1000.fasta")
        return (s1, s2)

    def test_seqio(self, tmpdir, data):
        s1, s2 = data
        logging.info("Testing `SequenceData` (IO)")
        kwargs = {"out_path": tmpdir, "tmp_path": tmpdir}
        d0 = SequenceData(s1, cds=True, to_stop=True, **kwargs)
        d1 = SequenceData(s1, cds=False, to_stop=False, **kwargs)
        d2 = SequenceData(s2, cds=False, to_stop=False, **kwargs)
        assert len(d0.cds_seqs) == 224  # Should only read proper CDS
        assert len(d1.cds_seqs) == 500  # Should read full file
        assert len(d2.pro_seqs) == 500  # Proper translation

    def test_paranome(self, data, tmpdir):
        logging.info("Testing paranome (requires diamond+mcl)")
        s1, s2 = data
        d = SequenceData(s1, out_path=tmpdir, tmp_path=tmpdir, cds=False, to_stop=True)
        d.get_paranome()
        # MCL is not deterministic, but the number of families is somewhere around 35
        assert 30 < len(d.mcl) < 40  

        # round tripping for families?
        # write to file/read from file
        families = d.write_paranome()
        fams = read_gene_families(families)
        fams = get_gene_families(d, fams)
        assert len(d.mcl) == len(fams)  # Round tripped

        # construct gene families directly from `self.mcl`
        #fams2 = get_gene_families(d, d.mcl.values(), rename=False)
        #assert len(fams2) == len(fams)

    def test_rbh(self, data, tmpdir):
        logging.info("Testing RBH orthologs (requires diamond)")
        s1, s2 = data
        kwargs = {"out_path": tmpdir, "tmp_path": tmpdir, "cds":False, "to_stop":False}
        d1 = SequenceData(s1, **kwargs)
        d2 = SequenceData(s2, **kwargs)
        d1.get_rbh_orthologs(d2, eval=1e-3)
        df = d1.rbh[d2.prefix]
        # test below fails with PAML 4.9i (works with 4.9j)
        # test below fails with diamond 0.9.18 (works with 2.0.5)
        assert len(df.index) == 19  # Got all RBHs 
        assert len(df.columns) == 12  # Right No. columns

    def test_ksd(self, data, tmpdir):
        s1, s2 = data
        d = SequenceData(s1, out_path=tmpdir, tmp_path=tmpdir, cds=True, to_stop=True)
        d.get_paranome()

        # TODO: split in subtests (functions)
        # alignment without gap-stripping
        families = d.write_paranome()
        fams = read_gene_families(families)
        gfams = get_gene_families(d, fams, prequal=False, strip_gaps=False)
        f = gfams[0]
        f.align()
        cds_len = f.cds_aln.get_alignment_length()
        pro_len = f.pro_aln.get_alignment_length()
        # test below fails with PAML 4.9i (works with 4.9j)
        # test below fails with diamond 0.9.18 (works with 2.0.5)
        assert cds_len == 558  # Alignment working?
        assert cds_len == 3*pro_len  # "AA -> codon fine?

        # alignment with gap-stripping
        gfams = get_gene_families(d, fams, prequal=False, strip_gaps=True)
        f = gfams[0]
        f.align()
        cds_len = f.cds_aln.get_alignment_length()
        pro_len = f.pro_aln.get_alignment_length()
        gaps = 0
        for i in range(pro_len):
            if "-" in f.pro_aln[:,i]:
                gaps += 1
        assert cds_len == 3*(pro_len - gaps)

        # run Ks distribution construction
        ksdb = KsDistributionBuilder(gfams)
        ksdb.get_distribution()
        assert len(ksdb.df.index.unique()) == len(ksdb.df.index)
        assert pytest.approx(25., 1.) == ksdb.df["dS"].mean()
        assert pytest.approx(90., 1.) == ksdb.df["S"].mean()
        
        # pairwise mode
        gfams = get_gene_families(d, fams, pairwise=True,
                prequal=False, strip_gaps=False)
        ksdb = KsDistributionBuilder(gfams)
        ksdb.get_distribution()
        assert len(ksdb.df.index.unique()) == len(ksdb.df.index)
        assert pytest.approx(25., 1.) == ksdb.df["dS"].mean()
        assert pytest.approx(90., 1.) == ksdb.df["S"].mean()


# TODO: test CLI
