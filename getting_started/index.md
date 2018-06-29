---
layout: page
title: Getting Started
excerpt: >-
  Salmon is an easy-to-use and ultrafast program for quantification from RNA-seq
  data
modified: {}
image: null
published: true
---

## Getting started with Salmon

This brief tutorial will explain how you can get started using Salmon to quantify your RNA-seq data.  This tutorial will walk you through installing salmon, building an index on a transcriptome, and then quantifying some RNA-seq samples for downstream processing.

Index of this tutorial:

1. [Obtaining Salmon](#obtaining-salmon)
2. [Indexing the transcriptome](#indexing-txome)
3. [Obtaining the reads](#obtaining-reads)
4. [Quantifying the samples](#quantifying-samples)
5. [After quantification](#after-quantification)

### Obtaining Salmon {#obtaining-salmon}

Salmon is a free (both as in "free beer" and "free speech") software tool for estimating transcript-level abundance from RNA-seq read data.  It is developed openly on GitHub.  You can visit Salmon's GitHub page [here](https://github.com/COMBINE-lab/salmon), and check out the Salmon source code, feature requests, known issues etc.  


**The easiest way to install salmon** is likely via [bioconda](https://anaconda.org/bioconda/salmon).

```
$ conda config --add channels conda-forge
$ conda config --add channels bioconda
$ conda create -n salmon salmon
```
This will install the latest salmon in its own conda environment.  The environment can then be activated via:

```
$ conda activate salmon
```

Another way to get started with Salmon is to download the pre-compiled binaries for your platform from the [releases page](https://github.com/COMBINE-lab/salmon/releases).  
If you've downloaded a specific binary, you simply decompress it like so:

```
$ tar xzvf salmon-0.10.2_linux_x86_64.tar.gz
```

then, the binary will be located in the `bin` directory inside of the uncompressed folder; for example `salmon-0.10.2_linux_x86_64/bin/salmon` in the example above.  You can either run salmon directly using the full path, or place it into your PATH variable for easier execution.  The rest of the tutorial below will assume that you've placed the `salmon` executable in your path, so that simply running `salmon` will invoke the program.  You can test that salmon is running on your system and get a list of available commands using the `-h` command; you should see output like the following

```
$ salmon -h
salmon v0.10.2

Usage:  salmon -h|--help or
        salmon -v|--version or
        salmon -c|--cite or
        salmon [--no-version-check] <COMMAND> [-h | options]

Commands:
     index Create a salmon index
     quant Quantify a sample
     alevin single cell analysis
     swim  Perform super-secret operation
     quantmerge Merge multiple quantifications into a single file
```

**Note**: OSX is frustratingly particular about how it looks for dynamic symbols in programs. If you're on OSX and you're getting an unresolved symbol error, you should run Salmon with the library directory in you `DYLD_FALLBACK_LIBRARY_PATH`, like this:

```
$ DYLD_FALLBACK_LIBRARY_PATH=<path_to_salmon_folder>/lib ./salmon -h
```

now, Salmon should find the appropriate symbols.

#### Obtaining a docker image for Salmon

Salmon is also available via Docker hub.  You can obtain a docker image of salmon using the command:

```
$ docker pull combinelab/salmon:latest
```

Then, if you wish, you can follow the tuorital below using this contanerized version of Salmon.

### Analyzing your RNA-seq data with Salmon

#### Obtaining a transcriptome and building an index {#indexing-txome}

In order to quantify transcript-level abundances, Salmon requires a target *transcriptome*.  This transcriptome is given to Salmon in the form of a (possibly compressed) multi-FASTA file, with each entry providing the sequence of a transcript[^1].  For this example, we'll be analyzing some *Arabidopsis thaliana* data, so we'll download and index the *A. thaliana* transcriptome.  First, create a directory where we'll do our analysis, let's call it `salmon_tutorial`:

```
$ mkdir salmon_tutorial
$ cd salmon_tutorial
```

Now, download the transcriptome:

```
$ curl ftp://ftp.ensemblgenomes.org/pub/plants/release-28/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.28.cdna.all.fa.gz -o athal.fa.gz
```

Here, we've used a reference transcriptome for *Arabadopsis*.  However, one of the benefits of performing quantification directly on the transcriptome (rather than via the host genome), is that one can easily quantify assembled transcripts as well (obtained via software such as [StringTie](https://ccb.jhu.edu/software/stringtie/) for organisms with a reference or [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) for *de novo* RNA-seq experiments).

Next, we're going to build an *index* on our transcriptome.  The index is a structure that salmon uses to [quasi-map](http://bioinformatics.oxfordjournals.org/content/32/12/i192.abstract) RNA-seq reads during quantification.  The index need only be constructed once per transcriptome, and it can then be reused to quantify many experiments.  We use the *index* command of salmon to build our index:

```
$ salmon index -t athal.fa.gz -i athal_index
```

There are a number of different options you can pass to the indexer to change its behavior (read more about those [here](http://salmon.readthedocs.io/en/latest/)), but the default should work well for most data.

#### Obtaining sequencing data {#obtaining-reads}

In addition to the *index*, salmon obviously requires the RNA-seq reads from the experiment to perform quantification.  In this tutorial, we'll be analyzing data from [this 4-condition experiment](https://www.ebi.ac.uk/ena/data/view/DRP001761) [accession PRJDB2508].  You can use the following shell script to obtain the raw data and place the corresponding read files in the proper locations.  Here, we're simply placing all of the data in a directory called `data`, and the left and right reads for each sample in a sub-directory labeled with that sample's ID (i.e. `DRR016125_1.fastq.gz` and `DRR016125_2.fastq.gz` go in a folder called `data/DRR016125`).

```
#!/bin/bash
mkdir data
cd data
for i in `seq 25 40`; 
do 
  mkdir DRR0161${i}; 
  cd DRR0161${i}; 
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR016/DRR0161${i}/DRR0161${i}_1.fastq.gz; 
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR016/DRR0161${i}/DRR0161${i}_2.fastq.gz; 
  cd ..; 
done
cd .. 
```

We'll place these commands in a script called [`dl_tut_reads.sh`](https://raw.githubusercontent.com/COMBINE-lab/salmon/gh-pages/assets/dl_tut_reads.sh).  To download the data, just run the script and wait for it to complete:

```
$ bash dl_tut_reads.sh
```

*Now might be a good time to grab a cup of coffee (or tea)*.

#### Quantifying the samples {#quantifying-samples}

Now that we have our index built and all of our data downloaded, we're ready to quantify our samples.  Since we'll be running the same command on each sample, the simplest way to automate this process is, again, a simple shell script ([`quant_tut_samples.sh`](https://raw.githubusercontent.com/COMBINE-lab/salmon/gh-pages/assets/quant_tut_samples.sh)):

```
#!/bin/bash
for fn in data/DRR0161{25..40};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i athal_index -l A \
         -1 ${fn}/${samp}_1.fastq.gz \
         -2 ${fn}/${samp}_2.fastq.gz \
         -p 8 -o quants/${samp}_quant
done 
```

This script simply loops through each sample and invokes `salmon` using fairly barebone options.  The `-i` argument tells salmon where to find the index `-l A` tells salmon that it should automatically determine the library type of the sequencing reads (e.g. stranded vs. unstranded etc.).  The `-1` and `-2` arguments tell salmon where to find the left and right reads for this sample (notice, salmon will accept gzipped FASTQ files directly).  Finally, the `-p 8` argument tells salmon to make use of 8 threads and the `-o` argument specifies the directory where salmon's quantification results sould be written.  Salmon exposes *many* different options to the user that enable extra features or modify default behavior.  However, the purpose and behavior of all of those options is beyond the scope of this introductory tutorial.  You can read about salmon's many  options in the [documentation](http://salmon.readthedocs.io/en/latest/).

After the salmon commands finish running, you should have a directory named `quants`, which will have a sub-directory for each sample.  These sub-directories contain the quantification results of salmon, as well as a lot of other information salmon records about the sample and the run.  The main output file (called `quant.sf`) is rather self-explanatory.  For example, take a peek at the quantification file for sample `DRR016125` in `quants/DRR016125/quant.sf` and you'll see a simple TSV format file listing the name (`Name`) of each transcript, its length (`Length`), effective length (`EffectiveLength`) (more details on this in the documentation), and its abundance in terms of Transcripts Per Million (`TPM`) and estimated number of reads (`NumReads`) originating from this transcript.

#### After quantification {#after-quantification}

That's it!  Quantifying your RNA-seq data with salmon is that simple (and fast).
Once you have your quantification results you can use them for downstream analysis with differential expression tools like 
[DESeq2](https://bioconductor.org/packages/DESeq2),
[edgeR](https://bioconductor.org/packages/edgeR),
[limma](https://bioconductor.org/packages/limma), or
[sleuth](http://pachterlab.github.io/sleuth/).
Using the [tximport](http://bioconductor.org/packages/tximport) package, 
you can import salmon's transcript-level quantifications 
and optionally aggregate them to the gene level for gene-level differential expression analysis.
You can read more about how to import salmon's results into DESeq2 by reading the `tximport` section
of the excellent [DESeq2 vignette](https://bioconductor.org/packages/DESeq2).
For instructions on importing for use with edgeR or limma, see the 
[tximport vignette](http://bioconductor.org/packages/tximport).
For preparing salmon output for use with sleuth, 
see the [wasabi](https://github.com/COMBINE-lab/wasabi) package.

[^1]:
	When you are building a salmon index, **please do not build the index on the genome of the organism whose transcripts you want to quantify**, this is almost certainly not want you want to do and will not provide you with meaningful results.
