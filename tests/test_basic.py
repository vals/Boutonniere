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

    def setUpClass(self):
        self.oneread = os.path.join(os.path.dirname(__file__),
                                    "data", "1read.fastq")
        self.nomatch = os.path.join(os.path.dirname(__file__),
                                    "data", "nomatch.fastq")


        self.bloom_fh, self.tmpfile_path = tempfile.mkstemp(suffix=".bloom", dir=os.path.dirname(self.oneread))
        print "Created bloom filter at %s ..." % self.tmpfile_path

    def test_match_oneread(self):
        bt.create_ref_bloom_filter(self.oneread, 0.0005, self.tmpfile_path)
        assert os.path.exists(self.tmpfile_path)
        assert os.path.getsize(self.tmpfile_path) > 0

    def test_query_oneread(self):
        matched = bt.count_matches(self.oneread, self.tmpfile_path,
                                   bt.ratio(bt.total_reads(self.oneread), 10))
        assert matched == 1

    def test_query_nomatch(self):
        print self.tmpfile_path
        matched = bt.count_matches(self.nomatch, self.tmpfile_path,
                                   bt.ratio(bt.total_reads(self.nomatch), 10))
        assert matched == 0

    def tearDown(self):
        pass
        #os.unlink(self.tmpfile_path)
