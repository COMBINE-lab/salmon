.. _FileFormats:

Salmon Output File Formats
==========================

Quantification File
-------------------

Salmon's main output is its quantification file.  This file is a plain-text, tab-separated file
with a single header line (which names all of the columns).  This file is named ``quant.sf`` and
appears at the top-level of Salmon's output directory. The columns appear in the following order:

+------+--------+-----------------+----+----------+
| Name | Length | EffectiveLength |TPM | NumReads |
+------+--------+-----------------+----+----------+

Each subsequent row describes a single quantification record.  The columns have
the following interpretation.

* **Name** --- 
  This is the name of the target transcript provided in the input transcript database (FASTA file). 

* **Length** ---
  This is the length of the target transcript in nucleotides.

* **EffectiveLength** ---
  This is the computed *effective* length of the target transcript.  It takes into account 
  all factors being modeled that will effect the probability of sampling fragments from
  this transcript, including the fragment length distribution and sequence-specific and 
  gc-fragment bias (if they are being modeled).

* **TPM** ---
  This is salmon's estimate of the relative abundance of this transcript in units of Transcripts Per Million (TPM).
  TPM is the recommended relative abundance measure to use for downstream analysis. 

* **NumReads** --- 
  This is salmon's estimate of the number of reads mapping to each transcript that was quantified.  It is an "estimate" 
  insofar as it is the expected number of reads that have originated from each transcript given the structure of the uniquely 
  mapping and multi-mapping reads and the relative abundance estimates for each transcript.


Command Information File
------------------------

In the top-level quantification directory, there will be a file called ``cmd_info.json``.  This is a
JSON format file that records the main command line parameters with which Salmon was invoked for the 
run that produced the output in this directory.


Auxiliary Files
---------------

The top-level quantification directory will contain an auxiliary directory called ``aux_info`` (unless 
the auxiliary directory name was overridden via the command line).  This directory will have a number
of files (and subfolders) depending on how salmon was invoked.

""""""""""""""""
Meta information
""""""""""""""""

The auxiliary directory will contain a JSON format file called
``meta_info.json`` which contains meta information about the run,
including stats such as the number of observed and mapped fragments,
details of the bias modeling etc.  If Salmon was run with automatic
inference of the library type (i.e. ``--libType A``), then one
particularly important piece of information contained in this file is
the inferred library type.  Most of the information recorded in this
file should be self-descriptive.

"""""""""""""""""""""""""""""""
Unique and ambiguous count file
"""""""""""""""""""""""""""""""

The auxiliary directory also contains 2-column tab-separated file called
``ambig_info.tsv``. This file contains information about the number of
uniquely-mapping reads as well as the total number of ambiguously-mapping reads
for each transcript.  This file is provided mostly for exploratory analysis of
the results; it gives some idea of the fraction of each transcript's estimated
abundance that derives from ambiguously-mappable reads.

""""""""""""""""""""""""""""""
Observed library format counts
""""""""""""""""""""""""""""""

When run in *mapping-based* mode, the quantification directory will 
contain a file called ``lib_format_counts.json``.  This JSON file 
reports the number of fragments that had at least one mapping compatible 
with the designated library format, as well as the number that didn't.
It also records the strand-bias that provides some information about 
how strand-specific the computed mappings were.

Finally, this file contains a count of the number of *mappings* that
were computed that matched each possible library type.  These are
counts of *mappings*, and so a single fragment that maps to the
transcriptome in more than one way may contribute to multiple library
type counts. **Note**: This file is currently not generated when Salmon
is run in alignment-based mode.


""""""""""""""""""""""""""""
Fragment length distribution
""""""""""""""""""""""""""""

The auxiliary directory will contain a file called ``fld.gz``.  This
file contains an approximation of the observed fragment length
distribution.  It is a gzipped, binary file containing integer counts.
The number of (signed, 32-bit) integers (with machine-native
endianness) is equal to the number of bins in the fragment length
distribution (1,001 by default --- for fragments ranging in length
from 0 to 1,000 nucleotides).

""""""""""""""""""""""""""""
Sequence-specific bias files
""""""""""""""""""""""""""""

If sequence-specific bias modeling was enabled, there will be 4 files
in the auxiliary directory named ``obs5_seq.gz``, ``obs3_seq.gz``,
``exp5_seq.gz``, ``exp5_seq.gz``.  These encode the parameters of the
VLMM that were learned for the 5' and 3' fragment ends.  Each file
is a gzipped, binary file with the same format.

It begins with 3 32-bit signed integers which record the length of the
context (window around the read start / end) that is modeled, follwed
by the length of the context that is to the left of the read and the
length of the context that is to the right of the read.

Next, the file contains 3 arrays of 32-bit signed integers (each of which
have a length of equal to the context length recorded above).  The first
records the order of the VLMM used at each position, the second records
the *shifts* and the *widths* required to extract each sub-context --- these
are implementation details.

Next, the file contains a matrix that encodes all VLMM probabilities.
This starts with two signed integers of type ``std::ptrdiff_t``.  This
is a platform-specific type, but on most 64-bit systems should
correspond to a 64-bit signed integer.  These numbers denote the number of
rows (*nrow*) and columns (*ncol*) in the array to follow.

Next, the file contains an array of (*nrow* * *ncol*) doubles which
represent a dense matrix encoding the probabilities of the VLMM.  Each
row corresponds to a possible preceeding sub-context, and each column
corresponds to a position in the sequence context.  Unused values
(values where the length of the sub-context exceed the order of the
model at that position) contain a 0.  This array can be re-shaped
into a matrix of the appropriate size.

Finally, the file contains the marginalized 0:sup:`th`-order
probabilities (i.e. the probability of each nucleotide at each
position in the context).  This is stored as a 4-by-context length
matrix.  As before, this entry begins with two signed integers that
give the number of rows and columns, followed by an array of doubles
giving the marginal probabilities.  The rows are in lexicographic
order.

""""""""""""""""""""""
Fragment-GC bias files
""""""""""""""""""""""

If Salmon was run with fragment-GC bias correction enabled, the
auxiliary directory will contain two files named ``expected_gc.gz``
and ``observed_gc.gz``.  These are gzipped binary files containing,
respectively, the expected and observed fragment-GC content curves.
These files both have the same form.  They consist of a 32-bit signed
int, *dtype* which specifies if the values to follow are in
logarithmic space or not.  Then, the file contains two signed integers
of type ``std::ptrdiff`` which give the number of rows and columns of
the matrix to follow.  Finally, there is an array of *nrow* by *ncol*
doubles.  Each row corresponds to a conditional fragment GC
distribution, and the number of columns is the number of bins in the
learned (or expected) fragment-GC distribution.


.. _eq-class-file:

""""""""""""""""""""""
Equivalence class file
""""""""""""""""""""""

If Salmon was run with the ``--dumpEq`` option, then a file called ``eq_classes.txt``
will exist in the auxiliary directory.  The format of that file is as follows:


::
   
   N (num transcripts)
   M (num equiv classes)
   tn_1
   tn_2
   ...
   tn_N
   eq_1_size t_11 t_12 ... count
   eq_2_size t_21 t_22 ... count

   
That is, the file begins with a line that contains the number of
transcripts (say N) then a line that contains the number of
equivalence classes (say M). It is then followed by N lines that list
the transcript names --- the order here is important, because the
labels of the equivalence classes are given in terms of the ID's of
the transcripts. The rank of a transcript in this list is the ID with
which it will be labeled when it appears in the label of an
equivalence class. Finally, the file contains M lines, each of which
describes an equivalence class of fragments. The first entry in this
line is the number of transcripts in the label of this equivalence
class (the number of different transcripts to which fragments in this
class map --- call this k). The line then contains the k transcript
IDs. Finally, the line contains the count of fragments in this
equivalence class (how many fragments mapped to these
transcripts). The values in each such line are tab separated.


