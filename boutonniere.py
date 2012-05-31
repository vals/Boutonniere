from pybloomfilter import BloomFilter
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def create_ref_bloom_filter(reference_fastq_file, capacity, error_rate, bf_file):
    """From a given FASTQ reference sequence creates a bloom filter file
    from each read.
    """
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
    if subset_size > total_reads:
        return 1

    return int(total_reads / subset_size)


def count_matches(fastq_file, bf_file, ratio):
    """Goes through a fastq file and checks a sample of reads if they
    occur in the specified bloom filter.
    """
    bf = BloomFilter.open(bf_file)
    fastq_handle = open(fastq_file)

    fastq_it = FastqGeneralIterator(fastq_handle)
    checked = 0
    observed = 0
    ratio = int(ratio)
    for i, (_, read, _) in enumerate(fastq_it):
        if i % ratio:
            continue

        checked += 1
        if read in bf:
            observed += 1

    fastq_handle.close()
    bf.close()

    return checked, observed
