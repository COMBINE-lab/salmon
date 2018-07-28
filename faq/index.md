---
layout: page
title: Frequently Asked Questions
excerpt: >-
  salmon is an easy-to-use and ultrafast program for quantification from RNA-seq
  data
modified: {}
image: null
published: true
---
This page hosts a (potentially growing) list of frequently-asked questions about salmon.  If you think you have a common question, let us know, and we'll consider adding it to this page.

#### Q: What is salmon?

*A:* salmon is a program for quickly and accurately estimating transcript-level abundance from RNA-seq data.  There are more details about this in [about](https://combine-lab.github.io/salmon/about/) page.

#### Q: The results in salmon's quant.sf have fewer transcripts than the transcriptome I indexed.  Why is this?

*A* When salmon indexes transcripts it collapses _exact duplicates_ (i.e. transcripts with perfectly identical sequence).  For each group of sequence-identical transcripts, it will select only one representative to quantify.  Such sequence identical transcripts often arise from transcrips annotated on ALT contigs which nonetheless do not contain any variants with respect to the reference transcripts.  While salmon does collapse these exact duplicates, it records all of the decisions it makes.

  In the salmon index directory, there is a file called “duplicate_clusters.tsv”.  This is a 2-column tsv file where the first column lists a retained transcript and the right column lists a discarded duplicate.  A transcript can represent multiple discarded sequences, in which case it appears multiple times in the first column. This will allow you to see the sequence duplicate clusters. If quantified, all transcripts in these clusters would obtain identical quantification estimates.

  **If you really want to go through with quantification of sequence duplicates.** You can pass `--keepDuplicates` to the salmon indexing command. This will tell salmon not to discard these duplicates, and they will appear in the output quantifications.

#### Q: Can I use salmon to quantify my reads against a target genome instead of a transcriptome?

*A:* No.  salmon quantifies a provided set of target transcripts given an RNA-seq sample (similar to tools like RSEM and eXpress).  If you have a well-annotated genome, consider using the reference set of transcripts.  It is also possible to use another program like StringTie, TransComb, etc. to identify novel transcripts, and then to add the corresponding assembled transcripts to the salmon index for subsequent quantification.  If you are working with a *de novo* assembled transcriptome, you can use salmon to quantify the assembled transcripts directly. 

#### Q: What is mapping validation?  Should I use it?

*A:* Alignment validation is an algorithm to increase the sensitivity and specificity of the mapping performed in salmon (when not using alignment-based mode).  During the initial mapping process, the stringency is slightly decreased, leading to more potential mapping locations being reported.  Subsequently, the mappings are _scored/validated_ by using the vectorized, banded extension-alignment algorithm of [ksw2](https://github.com/lh3/ksw2).  When paired with salmon's [data-driven equivalence class factorization](https://academic.oup.com/bioinformatics/article/33/14/i142/3953977), this allows salmon to incorporate meaningful alignment scores to discard potentially spurious mappings, and to properly assign probabilites to high-quality mappings.  Salmon also implements a caching strategy to avoid the recomputation of redundant alignment scores, which turns out to be particularly important since most multimapping to the transcriptome results from exactly duplicated sequences.  Turning on mapping validation will lead to a moderate increase in the computational burden of the mapping step --- though it remains quite fast.  Unless you have a specific reason not to enable mapping validation, it is probably a good idea to turn it on.

#### Q: I see that salmon can either use quasi-mapping or alignments (computed with e.g., STAR or HISAT2). Which should I use?

*A:* It's really up to you!  If you haven't already aligned your reads, it's probably easiest to use salmon with quasi-mapping; this provides a very fast and accurate way to bypass the traditional alignment step of quantification.  However, if you've already aligned your reads (to the transcriptome), or if you have a specific scenario where you've aligned your reads in a special way, then you can provide salmon with your aligned reads.  If you do this, be sure that the reads are *not* sorted by position, and that all alignments for the same read appear consecutively in the alignment file.

#### Q: Does salmon support reading from gzipped files?

*A:* Yes.  As of version 0.7.0, salmon can read directly from gzipped FASTA/Q files.  This is true of both the indexing and quantification sub-commands.

#### Q: I'm currently using [Sailfish](https://github.com/kingsfordgroup/sailfish), should I switch to salmon?

*A:* Probably.  The majority of our development efforts have shifted to salmon, since it provides a richer framework and model for future extension.  Also, some aspects of the salmon model cannot be properly backported to Sailfish.  We recommend using salmon if there's nothing preventing you from doing so.

#### Q: How can I control the speed of salmon's bias correction?

*A:* salmon's fragment-GC bias correction, evauates the potential bias attributable to every position on every transcript from fragments of different potential lengths (which occurs with a non-trivial probability).  Sampling the fragment distribution at different frequencies allows to control the computational efficiency of this procedure at a potential cost to the bias modeling quality.  By default, salmon (>= v0.10) samples every 5th fragment size, and interpolates between these values. The `--biasSpeedSamp` option controls the sampling rate, sampling from the set of fragment sizes and interpolating intermediate results.  This can speed up fragment-GC bias correction by a multiplicative factor (at the cost of potentially reducing the fidelity of the bias model).  Generally, the effect of such sampling tends to be small, and sampling parameters of 10 or less should give very similar results to the full model.  Additionally, you can reduce the amount of memory salmon uses to store pre-computed GC counts along transcripts with the `--reduceGCMemory` option.  This option has no effect on accuracy, but can reduce memory usage (for a _very small_ potential increase in runtime).
