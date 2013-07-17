Boutonniere
===========

INSTALL

pip install pybloomfiltermmap

Create bloom filter:

$ python boutonniere.py --seq test/1read.fastq --build --out 1read.bloom


Query bloom filter:

$ python boutonniere.py --bloom 1read.bloom --seq test/1read.fastq
