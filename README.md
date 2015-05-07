[![Build Status](https://travis-ci.org/COMBINE-lab/salmon.svg?branch=master)](https://travis-ci.org/COMBINE-lab/salmon)
[![Documentation Status](https://readthedocs.org/projects/salmon/badge/?version=latest)](https://readthedocs.org/projects/salmon/?badge=latest)

What is Salmon?
===============

Salmon is a **wicked**-fast program to produce a highly-accurate, transcript-level quantification estimates from 
RNA-seq data.  Salmon achieves is accuracy and speed via a number of different innovations, including the 
use of *lightweight* alignments (accurate but fast-to-compute proxies for traditional read alignments) and 
massively-parallel stochastic collapsed variational inference.  The result is a versatile tool that fits nicely
into many differnt pipelines.  For example, you can choose to make use of our *lightweight* alignments by providing Salmon with raw sequencing reads, or, if it is more convenient, you can provide Salmon with regular alignments (e.g. 
computed with your favorite aligner), and it will use the same **wicked**-fast, state-of-the-art inference algorithm 
to estimate transcript-level abundances for your experiment.

Documentation
==============

The documentation for Salmon is available on [ReadTheDocs](http://readthedocs.org), check it out [here](http://salmon.readthedocs.org).

Chat live about Salmon
======================

You can chat with the Salmon developers and other users via Gitter!

[![Join the chat at https://gitter.im/COMBINE-lab/salmon](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/COMBINE-lab/salmon?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
