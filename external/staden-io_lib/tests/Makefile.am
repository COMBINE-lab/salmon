# Copyright (c) 2008-2009, 2013 Genome Research Ltd.
# Author(s): James Bonfield
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# 
#    1. Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.
# 
#    2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
# 
#    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
#    Institute nor the names of its contributors may be used to endorse
#    or promote products derived from this software without specific
#    prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
# LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
## Makefile.am -- Process this file with automake to produce Makefile.in

EXTRA_DIST              = $(TESTS) data compare_sam.pl generate_data.pl cram_io_test.c
MAINTAINERCLEANFILES    = Makefile.in

noinst_PROGRAMS = cram_io_test

test_outdir              = test.out

TESTS_ENVIRONMENT       = \
			top_builddir=$(top_builddir) \
			outdir=$(test_outdir)

TESTS                   = \
                        srf2fasta.test \
                        srf2fastq.test \
			srf_info.test \
			srf_filter.test \
			srf_index.test \
			scram.test \
			scram_mt.test \
			scram_v3.test \
			cram_io.test \
			java.test

cram_io_test_SOURCES = cram_io_test.c
cram_io_test_LDADD = $(top_builddir)/io_lib/libstaden-read.la

INCLUDES= -I${top_srcdir}

dist-hook:
	rm -rf `find $(distdir)/data -name .svn`
	rm -f $(distdir)/data/ce#sorted.*am $(distdir)/data/ce#unsorted.*am
	rm -f $(distdir)/data/*~
	chmod u+w $(distdir)/data

distclean-local:
	-rm -rf $(test_outdir)

