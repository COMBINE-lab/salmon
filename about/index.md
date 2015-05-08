---
layout: page
title: About Salmon
excerpt: "Salmon is an easy-to-use and ultrafast program for quantification from RNA-seq data"
modified: 2014-08-08T19:44:38.564948-04:00
image:
  feature: SalmonLogo.png
#  credit: WeGraphics
#  creditlink: http://wegraphics.net/downloads/free-ultimate-blurred-background-pack/
---

Salmon is a tool for**wicked-fast** transcript quantification from RNA-seq data.
It requires a set of target transcripts (either from a reference or de-novo
assembly) to quantify. All you need to run Salmon is a FASTA file containing
your reference transcripts and a (set of) FASTA/FASTQ file(s) containing your
reads. Optionally, Salmon can make use of pre-computed alignments (in the form
of a SAM/BAM file) to the transcripts rather than the raw reads.

The*lightweight-alignment-based* mode of Salmon runs in two phases; indexing and
quantification. The indexing step is independent of the reads, and only need to
be run one for a particular set of reference transcripts. The quantification
step, obviously, is specific to the set of RNA-seq reads and is thus run more
frequently. For a more complete description of all available options in Salmon,
see below.

The *alignment-based* mode of Salmon does not require indexing. Rather, you can
simply provide Salmon with a FASTA file of the transcripts and a SAM/BAM file
containing the alignments you wish to use for quantification.

