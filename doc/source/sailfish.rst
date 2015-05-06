Sailfish
================

Sailfish is a tool for transcript quantification from RNA-seq data.  It
requires a set of target transcripts (either from a reference or *de-novo*
assembly) to quantify.  All you need to run sailfish is a fasta file containing
your reference transcripts and a (set of) fasta/fastq file(s) containing your
reads.  Sailfish runs in two phases; indexing and quantification.  The indexing
step is independent of the reads, and only needs to be run once for a particular
set of reference transcripts and choice of k (the k-mer size). The
quantification step, obviously, is specific to the set of RNA-seq reads and is
thus run more frequently. For a more complete description of all available
options in sailfish, see the manual.


Indexing
--------

To generate the sailfish index for your reference set of transcripts, you
should run the following command:

::

    > sailfish index -t <ref_transcripts> -o <out_dir> -k <kmer_len>


This will build a sailfish index for k-mers of length ``<kmer_len>`` for the
reference transcripts  provided in the file ``<ref_transcripts>`` and place the
index under the directory ``<out_dir>``.  There  are additional options that can
be passed to the sailfish indexer (e.g. the number of threads to use).  These
can be seen by executing the command ``sailfish index -h``.

Quantification
--------------

Now that you have generated the sailfish index (say that it's the directory
``<index_dir>`` --- this corresponds to the <out_dir> argument provided in the
previous step), you can quantify the transcript expression for a given set of
reads.  To perform the quantification, you run a command like the following:

::

    > sailfish quant -i <index_dir> -l "<libtype>" {-r <unmated> | -1 <mates1> -2 <mates2>} -o <quant_dir>

Where ``<index_dir>`` is, as described above, the location of the sailfish
index, ``<libtype>`` is a string describing the format of the fragment (read)
library (see :ref:`FragLibType`), ``<unmated>`` is a list of files
containing unmated reads, ``<mates{1,2}>`` are lists of files containg,
respectively, the first and second mates of paired-end reads. Finally,
``<quant_dir>`` is the directory where the output should be written. Just like the
indexing step, additional options are available, and can be viewed by running
``sailfish quant -h``.

When the quantification step is finished, the directory ``<quant_dir>`` will
contain a file named "quant.sf" (and, if bias correction is enabled, an
additional file names "quant_bias_corrected.sf").  This file contains the
result of the Sailfish quantification step.  This file contains a number of
columns (which are listed in the last of the header lines beginning with '#').
Specifically, the columns are (1) Transcript ID, (2) Transcript Length, (3)
Transcripts per Million (TPM), (4) Reads Per Kilobase per Million mapped reads
(RPKM), (5) K-mers Per Kilobase per Million mapped k-mers (KPKM), (6) Estimated
number of k-mers (an estimate of the number of k-mers drawn from this
transcript given the transcript's relative abundance and length) and (7)
Estimated number of reads (an estimate of the number of reads drawn from this
transcript given the transcript's relative abnundance and length).  The first
two columns are self-explanatory, the next four are measures of transcript
abundance and the final is a commonly used input for differential expression
tools.  The Transcripts per Million quantification number is computed as
described in [1]_, and is meant as an estimate of the number of transcripts, per
million observed transcripts, originating from each isoform.  Its benefit over
the K/RPKM measure is that it is independent of the mean expressed transcript
length (i.e. if the mean expressed transcript length varies between samples,
for example, this alone can affect differential analysis based on the K/RPKM.)
The RPKM is a classic measure of relative transcript abundance, and is an
estimate of the number of reads per kilobase of transcript (per million mapped
reads) originating from each transcript. The KPKM should closely track the
RPKM, but is defined for very short features which are larger than the chosen
k-mer length but may be shorter than the read length. Typically, you should
prefer the KPKM measure to the RPKM measure, since the k-mer is the most
natural unit of coverage for Sailfish.

References
----------

.. [1] Li, Bo, et al. "RNA-Seq gene expression estimation with read mapping uncertainty."
    Bioinformatics 26.4 (2010): 493-500.

.. _CMake : http://www.cmake.org
.. _Boost: http://www.boost.org
