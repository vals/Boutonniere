"""Simple tests creating and querying minimal bloom filters.
"""
import os
import unittest
import tempfile

import boutonniere as bt
#from pybloomfilter import BloomFilter

class BloomTest(unittest.TestCase):
    """ Bloom filter tests, basic creation and querying,
        no performance or stress tests.
    """

    @classmethod
    def setUpClass(self):
        self.oneread = os.path.join(os.path.dirname(__file__),
                                    "data", "1read.fastq")
        self.nomatch = os.path.join(os.path.dirname(__file__),
                                    "data", "nomatch.fastq")


        self.bloom_fh, self.tmpfile_path = tempfile.mkstemp(suffix=".bloom", dir=os.path.dirname(self.oneread))

    def test_match_oneread(self):
        bt.create_ref_bloom_filter(self.oneread, 0.0005, self.tmpfile_path)
        assert os.path.exists(self.tmpfile_path)
        assert os.path.getsize(self.tmpfile_path) > 0

    def test_query_nomatch(self):
        checked, matched  = bt.count_matches(self.nomatch, self.tmpfile_path,
                                   bt.sampling(bt.total_reads(self.nomatch), 1))
        assert checked == 1
        assert isinstance(matched, dict), "Not a dict !: {}".format(type(matched))
        assert matched[self.tmpfile_path] == 0, "Unexpected matches found {}".format(matched[self.tmpfile_path])

    def test_query_oneread(self):
        _, matched = bt.count_matches(self.oneread, self.tmpfile_path,
                                   bt.sampling(bt.total_reads(self.oneread), 1))
        assert matched[self.tmpfile_path] == 1

    @classmethod
    def tearDownClass(self):
        os.remove(self.tmpfile_path)
