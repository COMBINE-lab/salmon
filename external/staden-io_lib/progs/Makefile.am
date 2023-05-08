# Copyright (c) 2003-2009, 2011, 2013 Genome Research Ltd.
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
bin_PROGRAMS = convert_trace makeSCF extract_seq extract_qual extract_fastq index_tar scf_dump scf_info scf_update get_comment hash_tar hash_extract hash_list trace_dump hash_sff append_sff ztr_dump srf_dump_all srf_index_hash srf_extract_linear srf_extract_hash srf2fastq srf2fasta srf_filter srf_info srf_list hash_exp cram_dump cram_index scramble scram_merge scram_pileup scram_flagstat scram_test cram_size

convert_trace_SOURCES = convert_trace.c
convert_trace_LDADD = $(top_builddir)/io_lib/libstaden-read.la

get_comment_SOURCES = get_comment.c
get_comment_LDADD = $(top_builddir)/io_lib/libstaden-read.la

makeSCF_SOURCES = makeSCF.c
makeSCF_LDADD = $(top_builddir)/io_lib/libstaden-read.la

extract_seq_SOURCES = extract_seq.c
extract_seq_LDADD = $(top_builddir)/io_lib/libstaden-read.la

extract_qual_SOURCES = extract_qual.c
extract_qual_LDADD = $(top_builddir)/io_lib/libstaden-read.la

extract_fastq_SOURCES = extract_fastq.c
extract_fastq_LDADD = $(top_builddir)/io_lib/libstaden-read.la

index_tar_SOURCES = index_tar.c

scf_dump_SOURCES = scf_dump.c
scf_dump_LDADD = $(top_builddir)/io_lib/libstaden-read.la

scf_info_SOURCES = scf_info.c
scf_info_LDADD = $(top_builddir)/io_lib/libstaden-read.la

scf_update_SOURCES = scf_update.c
scf_update_LDADD = $(top_builddir)/io_lib/libstaden-read.la

hash_exp_SOURCES = hash_exp.c
hash_exp_LDADD = $(top_builddir)/io_lib/libstaden-read.la

hash_tar_SOURCES = hash_tar.c
hash_tar_LDADD = $(top_builddir)/io_lib/libstaden-read.la

hash_list_SOURCES = hash_list.c
hash_list_LDADD = $(top_builddir)/io_lib/libstaden-read.la

hash_sff_SOURCES = hash_sff.c
hash_sff_LDADD = $(top_builddir)/io_lib/libstaden-read.la

append_sff_SOURCES = append_sff.c
append_sff_LDADD = $(top_builddir)/io_lib/libstaden-read.la

hash_extract_SOURCES = hash_extract.c
hash_extract_LDADD = $(top_builddir)/io_lib/libstaden-read.la

trace_dump_SOURCES = trace_dump.c
trace_dump_LDADD = $(top_builddir)/io_lib/libstaden-read.la

ztr_dump_SOURCES = ztr_dump.c
ztr_dump_LDADD = $(top_builddir)/io_lib/libstaden-read.la

srf2fastq_SOURCES = srf2fastq.c
srf2fastq_LDADD = $(top_builddir)/io_lib/libstaden-read.la

srf2fasta_SOURCES = srf2fasta.c
srf2fasta_LDADD = $(top_builddir)/io_lib/libstaden-read.la

srf_dump_all_SOURCES = srf_dump_all.c
srf_dump_all_LDADD = $(top_builddir)/io_lib/libstaden-read.la

srf_extract_linear_SOURCES = srf_extract_linear.c
srf_extract_linear_LDADD = $(top_builddir)/io_lib/libstaden-read.la

srf_extract_hash_SOURCES = srf_extract_hash.c
srf_extract_hash_LDADD = $(top_builddir)/io_lib/libstaden-read.la

srf_index_hash_SOURCES = srf_index_hash.c
srf_index_hash_LDADD = $(top_builddir)/io_lib/libstaden-read.la

srf_filter_SOURCES = srf_filter.c
srf_filter_LDADD = $(top_builddir)/io_lib/libstaden-read.la

srf_info_SOURCES = srf_info.c
srf_info_LDADD = $(top_builddir)/io_lib/libstaden-read.la

srf_list_SOURCES = srf_list.c
srf_list_LDADD = $(top_builddir)/io_lib/libstaden-read.la

cram_dump_SOURCES = cram_dump.c
cram_dump_LDADD = $(top_builddir)/io_lib/libstaden-read.la

cram_size_SOURCES = cram_size.c
cram_size_LDADD = $(top_builddir)/io_lib/libstaden-read.la

cram_index_SOURCES = cram_index.c
cram_index_LDADD = $(top_builddir)/io_lib/libstaden-read.la

#cram_to_sam_SOURCES = cram_to_sam.c
#cram_to_sam_LDADD = $(top_builddir)/io_lib/libstaden-read.la
#
#sam_to_cram_SOURCES = sam_to_cram.c
#sam_to_cram_LDADD = $(top_builddir)/io_lib/libstaden-read.la
#
#sam_convert_SOURCES = sam_convert.c
#sam_convert_LDADD = $(top_builddir)/io_lib/libstaden-read.la

scramble_SOURCES = scramble.c
scramble_LDADD = $(top_builddir)/io_lib/libstaden-read.la

scram_merge_SOURCES = scram_merge.c
scram_merge_LDADD = $(top_builddir)/io_lib/libstaden-read.la

scram_pileup_SOURCES = scram_pileup.c scram_pileup.h
scram_pileup_LDADD = $(top_builddir)/io_lib/libstaden-read.la

scram_flagstat_SOURCES = scram_flagstat.c
scram_flagstat_LDADD = $(top_builddir)/io_lib/libstaden-read.la

scram_test_SOURCES = scram_test.c
scram_test_LDADD = $(top_builddir)/io_lib/libstaden-read.la

INCLUDES= -I${top_srcdir}
