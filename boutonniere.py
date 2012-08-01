import argparse

from pybloomfilter import BloomFilter
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import FastaIterator


def create_ref_bloom_filter(reference_file, error_rate, bf_file, format="fasta"):
    """From a given FASTA reference sequence creates a bloom filter file
    from each read.
    """

    if format == "fasta":
    	file_it = FastaIterator
        record = lambda it: (seq.seq for seq in it)
    elif format == "fastq":
        file_it = FastqGeneralIterator
        record = lambda it: (seq for _, seq, _ in it)

    capacity = total_reads(reference_file)
    with open(reference_file) as handle:
        it = file_it(handle)
        read_it = record(it)
        bf = BloomFilter(capacity, error_rate, bf_file)
        bf.update(read_it)
        bf.close()


def total_reads(fastq_file):
    """Get the total number of reads in a fastq file.
    (Just reads all the fastq tuples of the file and counts.)
    """

    total_num_reads = 0
    with open(fastq_file) as fastq_handle:
        fastq_it = FastqGeneralIterator(fastq_handle)
        for total_num_reads, _ in enumerate(fastq_it):
            pass

    return total_num_reads+1


def sampling(total_reads, subset_size):
    """Gives a ratio to be used when iterating over the fastq file to check
    for matches in the screen references. That is as in "Check every 1 in
    (sampling) reads".
    """
    if subset_size >= total_reads:
        return 1

    return int(total_reads / subset_size)


def count_matches(fastq_file, bf_files, sampling):
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
    sampling = int(sampling)
    for i, (_, read, _) in enumerate(fastq_it):
        if i+1 % sampling:
            continue

        checked += 1
        for bf_file in bf_files:
            if read in bf[bf_file]:
                observed[bf_file] += 1

    fastq_handle.close()

    return checked, observed


def main():
    parser = argparse.ArgumentParser(description="Screens for contamination " \
        "by checking if reads from a sample of a given fastq file matches " \
        "reads in a reference bloom filter.")

    parser.add_argument("--build", dest='build', action='store_true', \
                        default=False, \
                        help="Create a bloom filter file from the given" \
                        "reference fastq file.")

    parser.add_argument("--seq", dest='seq', action='store',
                        help="File to be screened, " \
                        "or a list of file to be indexed as " \
                        "bloom filters in the case of --build.")

    parser.add_argument("--errors", action='store', default=0.0005,
                        help="The rate of false positives in the created " \
                        "bloom filters, in the case of --build.")

    parser.add_argument("--out", dest='out', action='store', default=None,
                        help="The name of the bloom filter file to be " \
                        "created in the case of --build.")

    parser.add_argument("--bloom", dest='bloom', action='store', default=None,
                        nargs='+', help="The bloom filters to screen against.")

    args = parser.parse_args()



    if args.build:
        create_ref_bloom_filter(args.seq, 0.0005, args.out)

    if args.bloom:
    	print count_matches(args.seq, args.bloom, sampling(total_reads(args.seq), 10))

if __name__ == "__main__":
    main()
