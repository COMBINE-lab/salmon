Salmon
================

Salmon is a tool for **wicked-fast** transcript quantification from RNA-seq
data.  It requires a set of target transcripts (either from a reference or
*de-novo* assembly) to quantify.  All you need to run Salmon is a FASTA file
containing your reference transcripts and a (set of) FASTA/FASTQ file(s)
containing your reads.  Optionally, Salmon can make use of pre-computed
alignments (in the form of a SAM/BAM file) to the transcripts rather than the
raw reads.

The **lightweight-alignment**-based mode of Salmon runs in two phases; indexing and
quantification. The indexing step is independent of the reads, and only need to
be run one for a particular set of reference transcripts. The quantification
step, obviously, is specific to the set of RNA-seq reads and is thus run more
frequently. For a more complete description of all available options in Salmon,
see below.

The **alignment**-based mode of Salmon does not require indexing.  Rather, you can 
simply provide Salmon with a FASTA file of the transcripts and a SAM/BAM file
containing the alignments you wish to use for quantification.

Using Salmon
------------

As mentioned above, there are two "modes" of operation for Salmon.  The first,
requires you to build an index for the transcriptome, but then subsequently
processes reads directly.  The second mode simply requires you to provide a
FASTA file of the transcriptome and a ``.sam`` or ``.bam`` file containing a
set of alignments.

.. note:: Read / alignment order

    Salmon, like eXpress, uses a streaming inference method to perform 
    transcript-level quantification.  One of the fundamental assumptions 
    of such inference methods is that observations (i.e. reads or alignments)
    are made "at random".  This means, for example, that alignments should 
    **not** be sorted by target or position.  If your reads or alignments 
    do not appear in a random order with respect to the target transcripts,
    please randomize / shuffle them before performing quantification with 
    Salmon.

.. note:: Number of Threads

    The number of threads that Salmon can effectively make use of depends 
    upon the mode in which it is being run.  In alignment-based mode, the
    main bottleneck is in parsing and decompressing the input BAM file.
    We make use of the `Staden IO <http://sourceforge.net/projects/staden/files/io_lib/>`_ 
    library for SAM/BAM/CRAM I/O (CRAM is, in theory, supported, but has not been
    thorougly tested).  This means that multiple threads can be effectively used
    to aid in BAM decompression.  However, we find that throwing more than a 
    few threads at file decompression does not result in increased processing
    speed.  Thus, alignment-based Salmon will only ever allocate up to 4 threads
    to file decompression, with the rest being allocated to quantification.
    If these threads are starved, they will sleep (the quantification threads 
    do not busy wait), but there is a point beyond which allocating more threads
    will not speed up alignment-based quantification.  We find that allocating 
    8 --- 12 threads results in the maximum speed, threads allocated above this
    limit will likely spend most of their time idle / sleeping.

    For lightweight-alignment-based Salmon, the story is somewhat different.
    Generally, performance continues to improve as more threads are made
    available.  This is because the determiniation of the potential mapping
    locations of each read is, generally, the slowest step in
    lightweight-alignment-based quantification.  Since this process is
    trivially parallelizable (and well-parallelized within Salmon), more
    threads generally equates to faster quantification. However, there may
    still be a limit to the return on invested threads. Specifically, writing
    to the mapping cache (see `Misc`_ below) is done via a single thread.  With
    a huge number of quantification threads or in environments with a very slow
    disk, this may become the limiting step. If you're certain that you have
    more than the required number of observations, or if you have reason to
    suspect that your disk is particularly slow on writes, then you can disable
    the mapping cache (``--disableMappingCache``), and potentially increase the
    parallelizability of lightweight-alignment-based Salmon.

Lightweight-alignment-based mode (including quasimapping)
---------------------------------------------------------

One of the novel and innovative features of Salmon is its ability to accurately
quantify transcripts using *lightweight* alignments.  Lightweight alignments
are mappings of reads to transcript positions that are computed without
performing a base-to-base alignment of the read to the transcript.  Lightweight 
alignments are typically much faster to compute than traditional (or full)
alignments, and can sometimes provide superior accuracy by being more robust 
to errors in the read or genomic variation from the reference sequence.

Salmon currently supports two different methods for lightweight-alignment; 
SMEM-based mapping and quasi-mapping.  SMEM-based mapping is the original 
lightweight-alignment method used by Salmon, and quasi-mapping is a newer and 
considerably faster alternative.  Both methods are currently exposed via the 
same ``quant`` command, but the methods require different indices so that 
SMEM-based mapping cannot be used with a quasi-mapping index and vice-versa.

If you want to use Salmon in lightweight alignment-based mode, then you first
have to build an Salmon index for your transcriptome.  Assume that
``transcripts.fa`` contains the set of transcripts you wish to quantify. First,
you run the Salmon indexer:

::
    
    > ./bin/salmon index -q -k 31 -t transcripts.fa -i transcripts_index
    
This will build the quasi-mapping-based index, using an auxiliary k-mer hash
over k-mers of length 31.  While quasi-mapping will make used of arbitrarily 
long matches between the query and reference, the `k` size selected here will 
act as the *minimum* acceptable length for a valid match.  Thus, a smaller 
value of `k` may slightly improve sensitivty.  We find that a `k` of 31 seems
to work well for reads of 75bp or longer, but you might consider a smaller 
`k` if you plan to deal with shorter reads. Note that there is also a 
`k` parameter that can be passed to the ``quant`` command.  However, this has
no effect if one is using a quasi-mapping index, as the `k` value provided
during the index building phase overrides any `k` provided during
quantification in this case.

::
    
    > ./bin/salmon index -t transcripts.fa -i transcripts_index

This will build the SMEM-based mapping index.  Note that no value of `k` 
is given here.  However, the SMEM-based mapping index makes use of a parameter 
`k` that is passed in during the ``quant`` phase (the default value is `19`). 

Then, you can quantify any set of reads (say, paired-end reads in files
`reads1.fq` and `reads2.fq`) directly against this index using the Salmon
``quant`` command as follows:

::

    > ./bin/salmon quant -i transcripts_index -l <LIBTYPE> -1 reads1.fq -2 reads2.fq -o transcripts_quant

If you are using single-end reads, then you pass them to Salmon with 
the ``-r`` flag like:

::

    > ./bin/salmon quant -i transcripts_index -l <LIBTYPE> -r reads.fq -o transcripts_quant


This same ``quant`` command will work with either index (quasi-mapping or
SMEM-based), and Salmon will automatically determine the type of index being 
read and perform the appropriate lightweight mapping accordingly.

You can, of course, pass a number of options to control things such as the
number of threads used or the different cutoffs used for counting reads.
Just as with the alignment-based mode, after Salmon has finished running, there
will be a directory called ``salmon_quant``, that contains a file called
``quant.sf`` containing the quantification results.


Alignment-based mode
--------------------

Say that you've prepared your alignments using your favorite aligner and the
results are in the file ``aln.bam``, and assume that the sequence of the
transcriptome you want to quantify is in the file ``transcripts.fa``.  You
would run Salmon as follows:

::

    > ./bin/salmon quant -t transcripts.fa -l <LIBTYPE> -a aln.bam -o salmon_quant

The ``<LIBTYPE>`` parameter is described below and is shared between both modes
of Salmon.  After Salmon has finished running, there will be a directory called
``salmon_quant``, that contains a file called ``quant.sf``.  This contains the
quantification results for the run, and the columns it contains are similar to
those of Sailfish (and self-explanatory where they differ).

For the full set of options that can be passed to Salmon in its alignment-based
mode, and a description of each, run ``salmon quant --help-alignment``.

.. note:: Genomic vs. Transcriptomic alignments

    Salmon expects that the alignment files provided are with respect to the
    transcripts given in the corresponding fasta file.  That is, Salmon expects
    that the reads have been aligned directly to the transcriptome (like RSEM,
    eXpress, etc.) rather than to the genome (as does, e.g. Cufflinks).  If you
    have reads that have already been aligned to the genome, there are
    currently 3 options for converting them for use with Salmon.  First, you
    could convert the SAM/BAM file to a FAST{A/Q} file and then use the
    lightweight-alignment-based mode of Salmon described below.  Second, given the converted
    FASTA{A/Q} file, you could re-align these converted reads directly to the
    transcripts with your favorite aligner and run Salmon in alignment-based
    mode as described above.  Third, you could use a tool like `sam-xlate <https://github.com/mozack/ubu/wiki>`_
    to try and convert the genome-coordinate BAM files directly into transcript 
    coordinates.  This avoids the necessity of having to re-map the reads. However,
    we have very limited experience with this tool so far.

.. topic:: Multiple alignment files
    
    If your alignments for the sample you want to quantify appear in multiple 
    .bam/.sam files, then you can simply provide the Salmon ``-a`` parameter 
    with a (space-separated) list of these files.  Salmon will automatically 
    read through these one after the other quantifying transcripts using the 
    alignments contained therein.  However, it is currently the case that these
    separate files must (1) all be of the same library type and (2) all be
    aligned with respect to the same reference (i.e. the @SQ records in the 
    header sections must be identical).


What's this ``LIBTYPE``?
------------------------

Salmon, like sailfish, has the user provide a description of the type of
sequencing library from which the reads come, and this contains information
about e.g. the relative orientation of paired end reads.  However, we've
replaced the somewhat esoteric description of the library type with a simple
set of strings; each of which represents a different type of read library. This
new method of specifying the type of read library is being back-ported into
Sailfish and will be available in the next release.

The library type string consists of three parts: the relative orientation of
the reads, the strandedness of the library, and the directionality of the
reads.

The first part of the library string (relative orientation) is only provided if
the library is paired-end. The possible options are:

::

    I = inward
    O = outward
    M = matching

The second part of the read library string specifies whether the protocol is
stranded or unstranded; the options are:

::

    S = stranded
    U = unstranded

If the protocol is unstranded, then we're done.  The final part of the library
string specifies the strand from which the read originates in a strand-specific
protocol — it is only provided if the library is stranded (i.e. if the
library format string is of the form S).  The possible values are:

::

    F = read 1 (or single-end read) comes from the forward strand
    R = read 1 (or single-end read) comes from the reverse strand

An example of some library format strings and their interpretations are:

::

    IU (an unstranded paired-end library where the reads face each other)

::

    SF (a stranded single-end protocol where the reads come from the forward strand)

::

    OSR (a stranded paired-end protocol where the reads face away from each other,
         read1 comes from reverse strand and read2 comes from the forward strand)


.. note:: Strand Matching

    Above, when it is said that the read "comes from" a strand, we mean that
    the read should align with / map to that strand.  For example, for
    libraries having the ``OSR`` protocol as described above, we expect that
    read1 maps to the reverse strand, and read2 maps to the forward strand. 


For more details on the library type, see :ref:`FragLibType`. 

Output
------

Salmon writes its output in a simple tab-delimited file format.  Any line that begins 
with a ``#`` is a comment line, and can be safely ignored.  Salmon records the files
and options passed to it in comments at the top of its output file.  The last comment 
line gives the names of each of the data columns. The columns appear in the following order: 

+------+--------+-----+----------+
| Name | Length | TPM | NumReads |
+------+--------+-----+----------+

Each subsequent row described a single quantification record.  The columns have
the following interpretation.

* **Name** --- 
  This is the name of the target transcript provided in the input transcript database (FASTA file). 

* **Length** ---
  This is the length of the target transcript in nucleotides.

* **TPM** ---
  This is salmon's estimate of the relative abundance of this transcript in units of Transcripts Per Million (TPM).
  TPM is the recommended relative abundance measure to use for downstream analysis. 

* **NumReads** --- 
  This is salmon's estimate of the number of reads mapping to each transcript that was quantified.  It is an "estimate" 
  insofar as it is the expected number of reads that have originated from each transcript given the structure of the uniquely 
  mapping and multi-mapping reads and the relative abundance estimates for each transcript.  You can round these values 
  to the nearest integer and use them directly as input to count-based methods like 
  `Deseq2 <http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ and 
  `EdgeR <http://master.bioconductor.org/packages/release/bioc/html/edgeR.html>`_, among others.

Misc
----

Salmon deals with reading from compressed read files in the same way as
sailfish --- by using process substitution.  Say in the
lightweigh-alignment-based salmon example above, the reads were actually in the
files ``reads1.fa.gz`` and ``reads2.fa.gz``, then you'd run the following
command to decompress the reads "on-the-fly":

::

    > ./bin/salmon quant -i transcripts_index -l <LIBTYPE> -1 <(gzcat reads1.fa.gz) -2 <(gzcat reads2.fa.gz) -o transcripts_quant

and the gzipped files will be decompressed via separate processes and the raw
reads will be fed into salmon.

.. note:: The Mapping Cache 

    Salmon requires a specific number of observations (fragments) to
    be observed before it will report its quantification results.  If it 
    doesn't see enough fragments when reading through the read files the 
    first time, it will process the information again (don't worry; it's not 
    double counting. The results from the first pass essentially become 
    a "prior" for assigning the proper read counts in subsequent passes).

    The first time the file is processed, the set of potential mappings for
    each fragment is written to a temporary file in an efficient binary format
    --- this file is called the mapping cache.  As soon as the required number
    of obvservations have been seen, salmon stops writing to the mapping cache
    (ensuring that the file size will not grow too large).  However, for
    experiments with fewer than the required number of observations, the
    mapping cache is a significant optimization over reading through the raw
    set of reads multiple times.  First, the work of determining the potential
    mapping locations for a read is only performed once, during the inital pass
    through the file.  Second, since the mapping cache is implemented as a
    regular file on disk, the information contained within a file can be
    processed multiple times, even if the file itself is being produced via
    e.g. process substitution as in the example above.
    
    You can control the required number of observations and thus, indirectly,
    the maximum size of the mapping cache file, via the ``-n`` argument.
    Note that the cache itself is considered a "temporary" file, and it is
    removed from disk by salmon before the program terminates.  If you are
    certain that your read library is large enough that you will observe the
    required number of fragments in the first pass, or if you have some other 
    reason to avoid creating the temporary mapping cache, it can disabled with
    the ``--disableMappingCache`` flag.

**Finally**, the purpose of making this beta executable (as well as the Salmon
code) available is for people to use it and provide feedback.  A pre-print and
manuscript are in the works, but the earlier we get feedback, thoughts,
suggestions and ideas, the better!  So, if you have something useful to report
or just some interesting ideas or suggestions, please contact us
(`rob.patro@cs.stonybrook.edu` and/or `carlk@cs.cmu.edu`).  Also, please use
the same e-mail addresses to contact us with any *detailed* bug-reports (though
bug-support for these early beta versions may be slow).
