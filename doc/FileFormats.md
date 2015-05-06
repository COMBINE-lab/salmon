Transcript Index Format
=======================

The sorted list (array) is the structure used by our program
to count the k-mers from the reads and it relies on a transcript
index.  The index is (currently) simply a sorted array containing
all of the kmers, in encoded (uint64_t) numeric order, that were
seen in the trancsripts.  The file format is as follows

````
kmer_len[uint32_t]
num_kmers[uint32_t]
k_1[uint_64t] . . . k_{num_kmers}[uint64_t]
````