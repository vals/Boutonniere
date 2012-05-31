import argparse

from pybloomfilter import BloomFilter
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def create_ref_bloom_filter(reference_fastq_file, error_rate, bf_file):
    """From a given FASTQ reference sequence creates a bloom filter file
    from each read.
    """
    capacity = total_reads(reference_fastq_file)
    with open(reference_fastq_file) as fastq_handle:
        fastq_it = FastqGeneralIterator(fastq_handle)
        read_it = (read for _, read, _ in fastq_it)
        bf = BloomFilter(capacity, error_rate, bf_file)
        bf.update(read_it)
        bf.close()


def total_reads(fastq_file):
    """Get the total number of reads in a fastq file.
    (Just reads all the fastq tuples of the file and counts.)
    """
    with open(fastq_file) as fastq_handle:
        fastq_it = FastqGeneralIterator(fastq_handle)
        for total_num_reads, _ in enumerate(fastq_it):
            pass

    return total_num_reads


def ratio(total_reads, subset_size):
    """Gives a ratio to be used when iterating over the fastq file to check
    for matches in the screen references. That is as in "Check every 1 in
    (ratio) reads".
    """
    if subset_size >= total_reads:
        return 1

    return int(total_reads / subset_size)


def count_matches(fastq_file, bf_files, ratio):
    """Goes through a fastq file and checks a sample of reads if they
    occur in the specified bloom filter.
    """
    if isinstance(bf_files, basestring):
        bf_files = [bf_files]

    bf = {}
    observed = {}
    for bf_file in bf_files:
        bf[bf_file] = BloomFilter.open(bf_file)
        observed[bf_file] = 0

    fastq_handle = open(fastq_file)
    fastq_it = FastqGeneralIterator(fastq_handle)
    checked = 0
    ratio = int(ratio)
    for i, (_, read, _) in enumerate(fastq_it):
        if i % ratio:
            continue

        checked += 1
        for bf_file in bf_files:
            if read in bf[bf_file]:
                observed[bf_file] += 1

    fastq_handle.close()
    bf.close()

    return checked, observed


def main():
    parser = argparse.ArgumentParser(description="Screens for contamination " \
        "by checking if reads from a sample of a given fastq file matches " \
        "reads in a reference bloom filter.")

    parser.add_argument("--build", dest='build', action='store_true', \
                        default=False, \
                        help="Create a bloom filter file from the given" \
                        "reference fastq file.")

    parser.add_argument('--seq', dest='sequences', action='store',
                        help="File to be screened, " \
                        "or a list of file to be indexed as " \
                        "bloom filters in the case of --build.")

    parser.add_argument("--errors", action='store', default=0.0005,
                        help="The rate of false positives in the created " \
                        "bloom filters, in the case of --build.")

    parser.add_argument("--out", action='store', default=None,
                        help="The name of the bloom filter file to be " \
                        "created in the case of --build.")

    parser.add_argument("--bloom", action='store', default=None, nargs='+',
                        help="The bloom filters to screen against.")

    args = parser.parse_args()


if __name__ == "__main__":
    main()
