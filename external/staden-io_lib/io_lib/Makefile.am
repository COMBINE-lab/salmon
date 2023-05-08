# Copyright (c) 2003, 2005-2007, 2009-2010, 2013 Genome Research Ltd.
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
lib_LTLIBRARIES = libstaden-read.la

libstaden_read_la_SOURCES = \
	Read.c \
	Read.h \
	scf_extras.c \
	scf_extras.h \
	translate.c \
	translate.h \
	compression.c \
	compression.h \
	ztr.c \
	ztr.h \
	ztr_translate.c \
	deflate_interlaced.c \
	deflate_interlaced.h \
	abi.h \
	fpoint.c \
	fpoint.h \
	seqIOABI.c \
	seqIOABI.h \
	seqIOALF.c \
	alf.h \
	expFileIO.c \
	expFileIO.h \
	plain.h \
	seqIOPlain.c \
	misc_scf.c \
	read_scf.c \
	scf.h \
	write_scf.c \
	sff.c \
	sff.h \
	srf.c \
	srf.h \
	array.c \
	array.h \
	compress.c \
	compress.h \
	error.c \
	error.h \
	files.c \
	find.c \
	mach-io.c \
	mach-io.h \
	misc.h \
	open_trace_file.c \
	open_trace_file.h \
	read_alloc.c \
	strings.c \
	tar_format.h \
	traceType.c \
	traceType.h \
	xalloc.c \
	xalloc.h \
	vlen.c \
	vlen.h \
	hash_table.c \
	hash_table.h \
	jenkins_lookup3.c \
	jenkins_lookup3.h \
	mFILE.c \
	mFILE.h \
	stdio_hack.h \
	pooled_alloc.c \
	pooled_alloc.h \
	bam.h \
	bam.c \
	sam_header.h \
	sam_header.c \
	cram.h \
	rANS_static.c \
	rANS_static.h \
	cram_codecs.c \
	cram_codecs.h \
	cram_encode.c \
	cram_encode.h \
	cram_decode.c \
	cram_decode.h \
	cram_stats.c \
	cram_stats.h \
	cram_io.c \
	cram_io.h \
	cram_index.c \
	cram_index.h \
	cram_structs.h \
	zfio.c \
	zfio.h \
	dstring.c \
	dstring.h \
	string_alloc.c \
	string_alloc.h \
	md5.c \
	md5.h \
	crc32.c \
	crc32.h \
	scram.c \
	scram.h \
	thread_pool.c \
	thread_pool.h \
	binning.h \
	binning.c \
	cram_bambam.c

libstaden_read_la_CPPFLAGS = @LIBCURL_CPPFLAGS@

libstaden_read_la_LDFLAGS = -version-info @VERS_CURRENT@:@VERS_REVISION@:@VERS_AGE@ 
libstaden_read_la_LIBADD = @LIBZ@ @LIBCURL@

INCLUDES = -I${top_srcdir}
