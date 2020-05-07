import unittest

class TestDiamondRelated(unittest.TestCase):
    def __init__(self):
        self.data1 = ""
        self.data2 = ""

    def test_seqdata(self):
        sd = SequenceData(self.data1, out_path=outdir, tmp_path=tmpdir)
