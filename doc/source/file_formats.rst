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


Auxiliary File
--------------

The top-level quantification directory will contain an auxiliary directory called ``aux`` (unless 
the auxiliary directory name was overridden via the command line).  This directory will have a number
of files (and subfolders) depending on how salmon was invoked.

"""""""""""""""""""""
``fld.gz``
"""""""""""""""""""""

This file contains an approximation of the observed fragment length distribution.  It is a gzipped, binary file containing integer counts.  The number of (signed, 32-bit) integers (with machine-native endianness) is equal to the number of bins in the fragment length distribution (1,001 by default --- for fragments ranging in length from 0 to 1,000 nucleotides).


