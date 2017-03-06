---
layout: page
title: Frequently Asked Questions
excerpt: >-
  Salmon is an easy-to-use and ultrafast program for quantification from RNA-seq
  data
modified: {}
image: null
published: true
---
This page hosts a (potentially growing) list of frequently-asked questions about Salmon.  If you think you have a common question, let us know, and we'll consider adding it to this page.

### Q: What is Salmon?

*A:* Salmon is a program for quickly and accurately estimating transcript-level abundance from RNA-seq data.  There are more details about this in [about](https://combine-lab.github.io/salmon/about/) page.

### Q: Can I use Salmon to quantify my reads against a target genome instead of a transcriptome?

*A:* No.  Salmon quantifies a provided set of target transcripts given an RNA-seq sample (similar to tools like RSEM and eXpress).  If you have a well-annotated genome, consider using the reference set of transcripts.  It is also possible to use another program like StringTie, TransComb, etc. to identify novel transcripts, and then to add the corresponding assembled transcripts to the Salmon index for subsequent quantification.  If you are working with a *de novo* assembled transcriptome, you can use Salmon to quantify the assembled transcripts directly. 

### Q: I see that Salmon can either use quasi-mapping or alignments (computed with e.g., STAR or HISAT2). Which should I use?

*A:* It's really up to you!  If you haven't already aligned your reads, it's probably easiest to use Salmon with quasi-mapping; this provides a very fast and accurate way to bypass the traditional alignment step of quantification.  However, if you've already aligned your reads (to the transcriptome), or if you have a specific scenario where you've aligned your reads in a special way, then you can provide Salmon with your aligned reads.  If you do this, be sure that the reads are *not* sorted by position, and that all alignments for the same read appear consecutively in the alignment file.

### Q: Does Salmon support reading from gzipped files?

*A:* Yes.  As of version 0.7.0, Salmon can read directly from gzipped FASTQ files.  This is true of both the indexing and quantification sub-commands.

### Q: How can I control the speed of Salmon's bias correction?

*A:* Salmon's fragment-GC bias correction, by default, evauates the potential bias attributable to every position on every transcript from a fragment of every potential length (which occurs with a non-trivial probability).  For very large transcriptomes this can be a lot of evaluation.  The `--biasSpeedSamp` option will instead sample from the set of fragment sizes and interpolate intermediate results.  This can speed up fragment-GC bias correction by a multiplicative factor (at the cost of potentially reducing the fidelity of the bias model).  Generally, the effect of such sampling tends to be small, and sampling parameters of 10 or less should give very similar results to the full model.  Likewise, you can reduce the amount of memory Salmon uses to store pre-computed GC counts along transcripts with the `--gcSizeSamp` option.  This option has no effect on accuracy, but a larger sample size will save memory at the cost of increasing runtime.
