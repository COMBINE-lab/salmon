# Copyright (c) 2003-2010, 2013 Genome Research Ltd.
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
AUTOMAKE_OPTIONS = foreign no-dependencies

SUBDIRS = io_lib progs tests

man_MANS = \
	man/man1/srf_index_hash.1 \
	man/man1/srf_info.1 \
	man/man1/srf_list.1 \
	man/man1/srf2fasta.1 \
	man/man1/srf2fastq.1 \
	man/man1/scramble.1 \
	man/man3/read_scf_header.3 \
	man/man3/read_allocate.3 \
	man/man3/ExperimentFile.3 \
	man/man3/write_reading.3 \
	man/man3/write_scf.3 \
	man/man3/read_deallocate.3 \
	man/man3/fwrite_reading.3 \
	man/man3/read2scf.3 \
	man/man3/scf2read.3 \
	man/man3/fwrite_scf.3 \
	man/man3/write_scf_header.3 \
	man/man3/read_reading.3 \
	man/man3/read_scf.3 \
	man/man3/exp2read.3 \
	man/man3/fread_reading.3 \
	man/man3/read2exp.3 \
	man/man3/fread_scf.3 \
	man/man4/Read.4

pkginclude_HEADERS = \
	io_lib/Read.h \
	io_lib/scf_extras.h \
	io_lib/translate.h \
	io_lib/abi.h \
	io_lib/fpoint.h \
	io_lib/seqIOABI.h \
	io_lib/alf.h \
	io_lib/expFileIO.h \
	io_lib/plain.h \
	io_lib/scf.h \
	io_lib/sff.h \
	io_lib/array.h \
	io_lib/compress.h \
	io_lib/error.h \
	io_lib/mach-io.h \
	io_lib/misc.h \
	io_lib/open_trace_file.h \
	io_lib/tar_format.h \
	io_lib/traceType.h \
	io_lib/xalloc.h \
	io_lib/mFILE.h \
	io_lib/stdio_hack.h \
	io_lib/vlen.h \
	io_lib/hash_table.h \
	io_lib/jenkins_lookup3.h \
	io_lib/os.h \
	io_lib/compression.h \
	io_lib/ztr.h \
	io_lib/deflate_interlaced.h \
	io_lib/srf.h \
	io_lib/pooled_alloc.h \
	io_lib/cram.h \
	io_lib/cram_structs.h \
	io_lib/cram_io.h \
	io_lib/cram_encode.h \
	io_lib/cram_decode.h \
	io_lib/cram_codecs.h \
	io_lib/cram_index.h \
	io_lib/cram_stats.h \
	io_lib/zfio.h \
	io_lib/scram.h \
	io_lib/bam.h \
	io_lib/sam_header.h \
	io_lib/dstring.h \
	io_lib/string_alloc.h \
	io_lib_config.h \
	io_lib/md5.h \
	io_lib/thread_pool.h \
	io_lib/binning.h

bin_SCRIPTS = io_lib-config

EXTRA_DIST = README COPYRIGHT ChangeLog CHANGES man options.mk bootstrap \
	docs/ZTR_format docs/Hash_File_Format io_lib-config.in io_lib/os.h.in

dist-hook:
	rm -rf `find $(distdir)/man -name .svn`
