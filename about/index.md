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

Salmon is a **wicked-fast** program to produce a highly-accurate, transcript-level quantification estimates from RNA-seq data. Salmon achieves is accuracy and speed via a number of different innovations, including the use of quasi-mapping (an accurate but fast-to-compute proxy for traditional read alignments) and a two-phase inference procedure that makes use of massively-parallel stochastic collapsed variational inference. The result is a versatile tool that fits nicely into many different pipelines. For example, you can choose to make use of quasi-mapping by providing Salmon with raw sequencing reads, or, if it is more convenient, you can provide Salmon with regular alignments (e.g. computed with your favorite aligner), and it will use the same wicked-fast, state-of-the-art inference algorithm to estimate transcript-level abundances for your experiment.

All you need to run Salmon is a FASTA file containing your reference transcripts and a (set of) FASTA/FASTQ file(s) containing your reads. Optionally, Salmon can make use of pre-computed alignments (in the form of a SAM/BAM file) to the transcripts rather than the raw reads.

The *quasi-mapping-based* mode of Salmon runs in two phases; indexing and quantification. The indexing step is independent of the reads, and only need to be run one for a particular set of reference transcripts. The quantification step, obviously, is specific to the set of RNA-seq reads and is thus run more frequently. For a more complete description of all available options in Salmon, see below.

The *alignment-based* mode of Salmon does not require indexing. Rather, you can simply provide Salmon with a FASTA file of the transcripts and a SAM/BAM file containing the alignments you wish to use for quantification.

