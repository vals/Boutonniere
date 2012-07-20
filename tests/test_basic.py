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

        self.test_files = {}

        self.test_files["oneread"] = {
                "ref_file": self.test_files[os.path.join(os.path.dirname(__file__), "data", "1read.fastq")],
                "bloomfilter": None
        }
        self.test_files["nomatch"] = {
                "ref_file": self.test_files[os.path.join(os.path.dirname(__file__), "data", "nomatch.fastq")],
                "bloomfilter": None
        }

        self.get_r = lambda name: self.test_files[name]["ref_file"]
        self.get_bf = lambda name: self.test_files[name]["bloomfilter"]

#        self.bloom_fh, self.tmpfile_path = tempfile.mkstemp(suffix=".bloom",
#                                       dir=os.path.dirname(self.test_files[]))

        for f in self.test_files.values():
            _, f["bloomfilter"] = tempfile.mkstemp(suffix=".bloom", dir=os.path.dirname(test_file))
            bt.create_ref_bloom_filter(f["ref_file"], 0.0005, f["bloomfilter"], format=test_file.split(".")[-1])

    def test_create_bloomfilter(self):
        assert os.path.exists(self.get_bf("oneread"))
        assert os.path.getsize(self.get_bf("oneread")) > 0

    def test_query_nomatch(self):
        checked, matched  = bt.count_matches(self.get_r("nomatch"), self.get_bf("nomatch"),
                                   bt.sampling(bt.total_reads(self.get_bf("nomatch")), 1))
        assert checked == 1
        assert isinstance(matched, dict), "Not a dict !: {}".format(type(matched))
        assert matched[self.get_bf("nomatch")] == 0, "Unexpected matches found {}".format(matched[self.get_bf("nomatch")])

    def test_query_oneread(self):
        _, matched = bt.count_matches(self.oneread, self.tmpfile_path,
                                   bt.sampling(bt.total_reads(self.test_files["ref_file"]+oneread), 1))
        assert matched[self.tmpfile_path] == 1

    @classmethod
    def tearDownClass(self):
        for tmp_file in self.test_files.values():
            os.remove(tmp_file)
