/*
 * Copyright (c) 2005-2007, 2013 Genome Research Ltd.
 * Author(s): James Bonfield
 * 
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are met:
 * 
 *    1. Redistributions of source code must retain the above copyright notice,
 *       this list of conditions and the following disclaimer.
 * 
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
 *    Institute nor the names of its contributors may be used to endorse
 *    or promote products derived from this software without specific
 *    prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
 * LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Author(s): James Bonfield
 * 
 * Copyright (c) 2001 MEDICAL RESEARCH COUNCIL
 * All rights reserved
 * 
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are met:
 * 
 *    1 Redistributions of source code must retain the above copyright notice, 
 *      this list of conditions and the following disclaimer.
 * 
 *    2 Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 * 
 *    3 Neither the name of the MEDICAL RESEARCH COUNCIL, THE LABORATORY OF
 *      MOLECULAR BIOLOGY nor the names of its contributors may be used
 *      to endorse or promote products derived from this software without
 *      specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef _OPEN_TRACE_FILE_H_
#define _OPEN_TRACE_FILE_H_

#include "io_lib/mFILE.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Tokenises the search path splitting on colons (unix) or semicolons
 * (windows).
 * We also  explicitly add a "./" to the end of the search path
 *
 * Returns: A new search path with items separated by nul chars. Two nul
 *          chars in a row represent the end of the tokenised path.
 * Returns NULL for a failure.
 *
 * The returned data has been malloced. It is up to the caller to free this
 * memory.
 */
char *tokenise_search_path(char *searchpath);

/*
 * Opens a trace file named 'file'. This is initially looked for as a
 * pathname relative to a file named "relative_to". This may (for
 * example) be the name of an experiment file referencing the trace
 * file. In this case by passing relative_to as the experiment file
 * filename the trace file will be picked up in the same directory as
 * the experiment file. Relative_to may be supplied as NULL.
 *
 * 'file' is looked for at relative_to, then the current directory, and then
 * all of the locations listed in 'path' (which is a colon separated list).
 * If 'path' is NULL it uses the RAWDATA environment variable instead.
 *
 * Returns a mFILE pointer when found.
 *           NULL otherwise.
 */
mFILE *open_path_mfile(char *file, char *path, char *relative_to);

/*
 * Returns a mFILE containing the entire contents of the url;
 *         NULL on failure.
 */
mFILE *find_file_url(char *file, char *url);


/*
 * Opens a trace file named 'file'. This is initially looked for as a
 * pathname relative to a file named "relative_to". This may (for
 * example) be the name of an experiment file referencing the trace
 * file. In this case by passing relative_to as the experiment file
 * filename the trace file will be picked up in the same directory as
 * the experiment file. Relative_to may be supplied as NULL.
 *
 * 'file' is looked for at relative_to, then the current directory, and then
 * all of the locations listed in RAWDATA (which is a colon separated list).
 *
 * Returns a mFILE pointer when found.
 *           NULL otherwise.
 */
mFILE *open_trace_mfile(char *file, char *relative_to);
FILE *open_trace_file(char *file, char *relative_to);
mFILE *open_exp_mfile(char *file, char *relative_to);
FILE *open_exp_file(char *file, char *relative_to);

void  iolib_set_trace_path(char *path);
char *iolib_get_trace_path(void);
void  iolib_set_exp_path  (char *path);
char *iolib_get_exp_path  (void);

#ifdef __cplusplus
}
#endif

#endif /* _OPEN_TRACE_FILE_H_ */
