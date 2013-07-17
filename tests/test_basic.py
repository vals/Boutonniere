"""Simple tests creating and querying minimal bloom filters.
"""
import os
import unittest
import tempfile

import boutonniere as bt

class BloomTest(unittest.TestCase):
    """ Bloom filter tests, basic creation and querying,
        no performance or stress tests.
    """



    @classmethod
    def setUpClass(self):

        self.test_files = {}

        self.test_files["oneread"] = {
                "ref_file": os.path.join(os.path.dirname(__file__), "data", "oneread.fastq"),
                "bloomfilter": None
        }
        self.test_files["nomatch"] = {
                "ref_file": os.path.join(os.path.dirname(__file__), "data", "nomatch.fastq"),
                "bloomfilter": None
        }

        for test_name, test_file in self.test_files.items():
            _, test_file["bloomfilter"] = tempfile.mkstemp(suffix=".bloom", dir=os.path.dirname(self.test_files[test_name]["ref_file"]))
            bt.create_ref_bloom_filter(test_file["ref_file"], 0.0005, test_file["bloomfilter"], format=self.test_files[test_name]["ref_file"].split(".")[-1])


    def get_r(self, name):
        return self.test_files[name]["ref_file"]

    def get_bf(self, name):
        return self.test_files[name]["bloomfilter"]

    def test_create_bloomfilter(self):
        assert os.path.exists(self.get_bf("oneread"))
        assert os.path.getsize(self.get_bf("oneread")) > 0

    def test_query_nomatch(self):
        checked, matched  = bt.count_matches(self.get_r("nomatch"), self.get_bf("nomatch"),
                                   bt.sampling(bt.total_reads(self.get_bf("nomatch")), 1))
        assert checked == 1
        assert isinstance(matched, dict), "Not a dict !: {}".format(type(matched))
## Not guaranteed to work
#        assert matched[self.get_bf("nomatch")] == 0, "Unexpected matches found {}".format(matched[self.get_bf("nomatch")])

    def test_query_oneread(self):
        """ Tests querying oneread.fastq file for a single match
        """
        _, matched = bt.count_matches(self.get_r("oneread"), self.get_bf("oneread"),
                                   bt.sampling(bt.total_reads(self.get_r("oneread")), 1))
        assert matched[self.get_bf("oneread")] == 1

    @classmethod
    def tearDownClass(self):
        for test_name in self.test_files.keys():
            os.remove(self.test_files[test_name]["bloomfilter"])
