from pybloomfilter import BloomFilter
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def create_reference_bloom_filter(fastq_file, capacity, error_rate, bf_file):
    """From a given FASTQ reference sequence creates a bloom filter file
    from each read.
    """
    with open(fastq_file) as fastq_handle:
        fastq_it = FastqGeneralIterator(fastq_handle)
        read_it = (read for _, read, _ in fastq_it)
        bf = BloomFilter(capacity, error_rate, bf_file)
        bf.update(read_it)
        bf.close()
