[![Build Status](https://travis-ci.org/COMBINE-lab/salmon.svg?branch=master)](https://travis-ci.org/COMBINE-lab/salmon)
[![Documentation Status](https://readthedocs.org/projects/salmon/badge/?version=latest)](http://salmon.readthedocs.org/en/latest)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/salmon/README.html)


**Try out alevin (salmon's single-cell processing module)!  Get started with the [tutorial](https://combine-lab.github.io/alevin-tutorial/#blog)**

**Help guide the development of Salmon, [take our survey](https://docs.google.com/forms/d/e/1FAIpQLSeWhBNE_fA_0uVHvbAlAulDmfmowv7rAYla879DZpqCARyRTQ/viewform)**

### Pre-computed decoy transcriptomes 

Below, you can find download some pre-processed transcriptomes for common species (if you have a particular species or a particular annotation request, let us know).  You can easily build your own decoy-enhanced transcriptome using the the `genereateDecoyTranscriptome.sh` script from the [SalmonTools](https://github.com/COMBINE-lab/SalmonTools) repository.  However, we are providing these for convenience.  Below, you can simply download the file in the "transcriptome with decoys" column, and unzip the archive with `tar xzvf`.  Each decompressed directory contains a `gentrome.fa` file and `decoys.txt` file that can be provided to `salmon` as:

```
$ salmon index -t gentrome.fa -d decoys.txt -i combined_index
```

| description |transcriptome with decoys  |   link to base genome | link to annotation |
| -------- | ------------- | ------------ | -------------|
| human (ensembl)      | [**human.tar.gz**](http://bit.ly/2HUU7S6) | [Homo_sapiens.GRCh38.dna.toplevel.fa.gz](http://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz) | [Homo_sapiens.GRCh38.91.gtf.gz](http://ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf.gz) |
| mouse (ensembl)     | [**mouse.tar.gz**](http://bit.ly/2Xoop4X) | [Mus_musculus.GRCm38.dna.toplevel.fa.gz](http://ftp.ensembl.org/pub/release-91/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz) | [Mus_musculus.GRCm38.91.gtf.gz](http://ftp.ensembl.org/pub/release-91/gtf/mus_musculus/Mus_musculus.GRCm38.91.gtf.gz) |
| drosophila (ensembl) | [**drosophila.tar.gz**](http://bit.ly/2KrlCnF) | [Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz](http://ftp.ensembl.org/pub/release-91/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz) | [Drosophila_melanogaster.BDGP6.91.gtf.gz](http://ftp.ensembl.org/pub/release-91/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.91.gtf.gz) |

What is Salmon?
===============

Salmon is a **wicked**-fast program to produce a highly-accurate, transcript-level quantification estimates from 
RNA-seq data.  Salmon achieves its accuracy and speed via a number of different innovations, including the 
use of *quasi-mapping* (accurate but fast-to-compute proxies for traditional read alignments), and 
massively-parallel stochastic collapsed variational inference.  The result is a versatile tool that fits nicely
into many different pipelines.  For example, you can choose to make use of our *quasi-mapping* algorithm by providing Salmon with raw sequencing reads, or, if it is more convenient, you can provide Salmon with regular alignments (e.g. an **unsorted** BAM file produced with your favorite aligner), and it will use the same **wicked**-fast, state-of-the-art inference algorithm 
to estimate transcript-level abundances for your experiment.

Give salmon a try!  You can find the latest binary releases [here](https://github.com/COMBINE-lab/salmon/releases).

The current version number of the master branch of Salmon can be found [**here**](http://combine-lab.github.io/salmon/version_info/latest)

**NOTE**: Salmon works by (quasi)-mapping sequencing reads directly to the *transcriptome*.  This means the Salmon index should be built on a set of target transcripts, **not** on the *genome* of the underlying organism.  If indexing appears to be taking a very long time, or using a tremendous amount of memory (which it should not), please ensure that you are not attempting to build an index on the genome of your organism!

Documentation
==============

The documentation for Salmon is available on [ReadTheDocs](http://readthedocs.org), check it out [here](http://salmon.readthedocs.org).

Salmon is, and will continue to be, [freely and actively supported on a best-effort basis](https://oceangenomics.com/about/#open).
If you need industrial-grade technical support, please consider the options at [oceangenomics.com/support](http://oceangenomics.com/support).

Chat live about Salmon
======================

You can chat with the Salmon developers and other users via Gitter!

[![Join the chat at https://gitter.im/COMBINE-lab/salmon](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/COMBINE-lab/salmon?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
