---
layout: page
title: Getting Started
excerpt: >-
  Salmon is an easy-to-use and ultrafast program for quantification from RNA-seq
  data
modified: {}
image:
  feature: SalmonLogo.png
published: true
---

## Getting started with Salmon

This brief tutorial will explain how you can get started using Salmon to quantify your RNA-seq data.  This tutorial will walk you through installing salmon, building an index on a transcriptome, and then quantifying some RNA-seq samples for downstream processing.

### Obtaining Salmon

Salmon is a free (both as in "free beer" and "free speech") software tool for estimating transcript-level abundance from RNA-seq read data.  It is developed openly on GitHub.  You can visit Salmon's GitHub page [here](https://github.com/COMBINE-lab/salmon), and check out the Salmon source code, feature requests, known issues etc.  However, the easiest way to get started with Salmon is to download the pre-compiled binaries for your platfor from the [releases page](https://github.com/COMBINE-lab/salmon/releases).  We provide binaries for both 64-bit Linux and MacOS.  

Once you've downloaded the appropriate binary (e.g. Salmon-0.7.2_linux_x86_64.tar.gz for a 64-bit Linux system), you simply decompress it like so:

```
$ tar xzvf Salmon-0.7.2_linux_x86_64.tar.gz
```

then, the binary will be located in the `bin` directory inside of the uncompressed folder; for example `Salmon-0.7.2_linux_x86_64/bin/salmon` in the example above.  You can either run salmon directly using the full path, or place it into your PATH variable for easier execution.  The rest of the tutorial below will assume that you've placed the `salmon` executable in your path, so that simply running `salmon` will invoke the program.  You can test that salmon is running on your system and get a list of available commands using the `-h` command; you should see output like the following

```
$ salmon -h
Salmon v0.7.2

Usage:  salmon -h|--help or
        salmon -v|--version or
        salmon -c|--cite or
        salmon [--no-version-check] <COMMAND> [-h | options]

Commands:
     cite  Show salmon citation information
     index Create a salmon index
     quant Quantify a sample
     swim  Perform super-secret operation
```

### Analyzing your RNA-seq data with Salmon

**Obtaining a transcriptome and building an index**

In order to quantify transcript-level abundances, Salmon requires a target *transcriptome*.  This transcriptome is given to Salmon in the form of a (possibly compressed) multi-FASTA file, with each entry providing the sequence of a transcript[^1]



[^1]:
	When you are building a salmon index, **please do not build the index on the genome of the organism whose transcripts you want to quantify**, this is almost certainly not want you want to do and will not provide you with meaningful results.
