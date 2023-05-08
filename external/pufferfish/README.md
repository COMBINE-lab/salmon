![C/C++ CI](https://github.com/COMBINE-lab/pufferfish/workflows/C/C++%20CI/badge.svg?branch=master)

# Table of Contents
 * [What is puffaligner?](#puffaligner)
 * [What is pufferfish?](#whatis)
 * [How to install?](#building)
 * [How to use?](#using)

## What is Puffaligner <a name="puffaligner"></a>

Puffaligner is a fast, sensitive and accurate aligner built on top of the [Pufferfish index](#whatis).
It tries to occupy a less-well-explored position in the space of read aligners, typically 
using more memory than BWT-based approaches (unless there are _highly_ repetitive references), but
considerably less than very fast but memory-hungry aligners like STAR.  Puffaligner is based on 
hashing relatively long seeds and then extending them to MEMs, and so it is very fast (typically 
much faster than approaches based on arbitrary pattern matching in the BWT).  It takes a 
seed -> chain -> align approach similar to [BWA-MEM](https://github.com/lh3/bwa) and [minimap2](https://github.com/lh3/minimap2).

It supports aligning to transcriptome as well as genome, but currently supports only contiguous 
alignments (i.e. spliced-alignment is not yet implemented).  However, the design means that 
one can use Puffalign to align reads to a collection of genomes (or/and transcriptomes), as well
as to a joint index that contains both the genome and the spliced transcripts.

Another feature of the aligner is introducing a list of decoys along with the main sequences to the index. A read is discarded if it aligns better to a decoy sequence rather than a sequence in the main list. One of the usecases of such feature is in improving the transcript alignment accuracy in case of retained introns, processed pseudogenes, etc..


The main steps in the aligner are:
1. find the first unmapped kmer from the read in the pufferfish index
2. extend the mapping to a [uni-MEM](https://github.com/HongzheGuo/deBGA) (Maximal Extended Match on a unitig) between read and index
3. repeat 1 and 2 for the next uni-MEM until reaching the end of the read
4. project uni-MEMs to reference-based MEMs and compact them.
5. find the best chain of the MEMs (adopted from [minimap2](https://github.com/lh3/minimap2) chaining)
6. align the gaps between the MEMs and at the edges of the read
7. find the best pair of the reads in case of paired-end
8. recover orphans

There are a series of heuristics and best-practices used to improve both the performance and accuracy of the results.
 
## What is Pufferfish <a name="whatis"></a>

**short answer** : Pufferfish is a new time and memory-efficient data structure for indexing a compacted, colored de Bruijn graph (ccdBG).  You can read more about pufferfish in the [paper](https://academic.oup.com/bioinformatics/article/34/13/i169/5045749), which appeared at ISMB 2018.

**long answer** : 
Though the de Bruijn Graph (dBG) has enjoyed tremendous popularity as an assembly and sequence comparison data structure, it has only relatively recently begun to see use as an index of the reference sequences (e.g. [deBGA](https://github.com/HongzheGuo/deBGA), [kallisto](https://github.com/pachterlab/kallisto)). Particularly, these tools index the compacted dBG (cdBG), in which all non-branching paths are collapsed into individual nodes and labeled with the string they spell out. This data structure is particularly well-suited for representing repetitive reference sequences, since a single contig in the cdBG represents all occurrences of the repeated sequence. The original positions in the reference can be recovered with the help of an auxiliary "contig table" that maps each contig to the reference sequence, position, and orientation where it appears as a substring. The [deBGA paper](https://academic.oup.com/bioinformatics/article-abstract/32/21/3224/2415064/deBGA-read-alignment-with-de-Bruijn-graph-based?redirectedFrom=fulltext) has a nice description how this kind of index looks (they call it a unipath index, because the contigs we index are unitigs in the cdBG), and how all the pieces fit together to be able to resolve the queries we care about.  Moreover, the cdBG can be built on multiple reference sequences (transcripts, chromosomes, genomes), where each reference is given a distinct color (or colour, if you're of the British persuasion). The resulting structure, which also encodes the relationships between the cdBGs of the underlying reference sequences, is called the compacted, colored de Bruijn graph (ccdBG).  This is not, of course, the only variant of the dBG that has proven useful from an indexing perspective. The (pruned) dBG has also proven useful as a graph upon which to build a path index of arbitrary variation / sequence graphs, which has enabled very interesting and clever indexing schemes like that adopted in [GCSA2](https://github.com/jltsiren/gcsa2) (which we won't discuss further here, but which I hope to cover in more detail in a future post).  Also, thinking about sequence search in terms of the dBG has led to interesting representations for variation-aware sequence search backed by indexes like the vBWT (implemented in the excellent [gramtools](https://github.com/iqbal-lab-org/gramtools) package).

While existing hash-based indices based on the cdBG (and ccdBG) are very efficient for search, they typically occupy a large amount of space in memory (both during construction and even when built). As a result, to make use of such data structures on large reference sequences (e.g., the human genome) or collections of reference sequences (e.g., in a metagenomic context), one typically requires a very large memory machine — if the structures can be built at all. Pufferfish implements a new and much more compact data structure for indexing the ccdBG. While maintaining very efficient queries, this allows Pufferfish to index reference sequences while reducing the memory requirements considerably (by an order-of-magnitude or more). This greatly reduces the memory burden for indexing reference sequences and makes it possible to build hash-based indexes of sequences of size that were not previously feasible.

**about pufferfish development:**
Currently, Pufferfish is the software implementing this efficient ccdBG index, and allowing point (i.e., k-mer) queries.  Pufferfish is under active development, but we want to be as open (and as useful to as many people) as possible early on. 


**branches:**
The **master** branch of pufferfish is _not_ necessarily stable, but it should, at any given time contain a working version of the index.  That is, breaking changes should not be pushed to master.  The **develop** branch of pufferfish is guaranteed to be neither stable nor working at any given point, but a best-faith effort will be made to not commit broken code to this branch.  For feature branches, all bets are off.

For more details about pufferfish, please check out our [paper](https://academic.oup.com/bioinformatics/article/34/13/i169/5045749), as well as the blog post [here](http://robpatro.com/blog/?p=494).

## How to Install <a name="building"></a>
To build the pufferfish/puffaligner do the following,

```
> git clone git@github.com:COMBINE-lab/pufferfish.git
> cd pufferfish
> mkdir build
> cd build
> cmake ../
> make
```

## How to Use <a name="using"></a>

**Programs used within pufferfish**

Building a pufferfish index requires first having a compacted de Bruijn graph, for which we use a modified version of [TwoPaCo](https://github.com/medvedevgroup/TwoPaCo). However, some modification of the TwoPaCo output is required for pufferfish to properly index the graph (e.g. a k-mer must appear at most once in the graph and palindromic contigs output by TwoPaCo must be removed). Thus we rely on a modified version of TwoPaCo which we bundle with pufferfish in the `external` directory.

To choose an appropriate filter size to pass to TwoPaCo to build the compacted dBG, we make use the the hyper-log-log implementation of [ntCard](https://github.com/bcgsc/ntCard). Because we use this as a library instead of an executable, and to avoid an external dependency to simply call one function, we bundle a modified version of that code with pufferfish and also include it in the `external` directory.

We are also dependent on [SeqLib](https://github.com/walaj/SeqLib) and hence all the libraries that it is dependent on such as `bz2`, `lzma`, and `z` for mapping part. So it is required to install these libraries on the system.
However, we also have the selected libraries from seqlib that we use bundled with pufferfish repo, so the installation should work without any difficulties.

#### Core Operations

**Building a pufferfish index**

To build a pufferfish index, you can use the `index` command.  It is used like so:

```
pufferfish index -r <fasta_file_to_index> -o <pufferfish index directory>
```

There are also optional parameters including `-k` (setting the kmer size -- default:31)
, `-s` (the ability to build a sparser and smaller index), `-p` (control the number of threads used during construction), and `-f` (to provide an explicit filter size for TwoPaCo dBG construction).

**Aligning via Puffaligner**

To align a set of paired-end reads to the reference one can use
the following command:

```
pufferfish align -i <pufferfish_index> -1 <readfile1> -2 <readfile2> -o <outputfile> 
```

The input read files can be compressed or uncompressed `fastq` files

Puffaligner can generate different types of output including [SAM format](https://samtools.github.io/hts-specs/SAMv1.pdf).
There is also an efficient binary format, which we call pam that can be generated using the option `-p`.

There are a variety of optional choices for changing the default thresholds for allowing more alignments, higher or lower scored alignments, only the best, or only one best alignment, orphans, discordants etc. 

---

***Pufferfish* is now the main (and only) index used in [Salmon](https://github.com/COMBINE-lab/salmon.git) when
it is run in mapping-based mode (i.e. with selective-alignment).**
