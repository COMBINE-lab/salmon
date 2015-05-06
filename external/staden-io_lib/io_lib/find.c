/*
 * Copyright (c) 2007, 2010 Genome Research Ltd.
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
 * Author(s): James Bonfield, John Taylor
 * 
 * Copyright (c) 1996, 1999-2000 MEDICAL RESEARCH COUNCIL
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

/*
 * Copyright (c) Medical Research Council 1994. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies.
 *
 * This file was written by James Bonfield, Simon Dear, Rodger Staden,
 * as part of the Staden Package at the MRC Laboratory of Molecular
 * Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
 *
 * MRC disclaims all warranties with regard to this software.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include "io_lib/misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* 19/3/99 johnt - added Corba support */
#ifdef USE_CORBA
#include "stcorba.h"
#endif

#ifdef USE_BIOLIMS
#include "spBiolims.h"
#endif

char *myfind(char *file, char* searchpath, int (*found) (char *) )
{
    static char wholePath[1024];
    char *path;
    char *f;

    f = NULL;
    if (found(file)) {
	strcpy(wholePath,file);
	f = wholePath;
    } else if (searchpath != NULL) {
	char *paths;
	char *next;

	paths = (char *) malloc(strlen(searchpath)+1);
	strcpy(paths,searchpath);

	path = paths;
	next = strchr(path,':');
        while( next && (*(next+1) == ':' )){
    	   /* 26/03/99 johnt - allow : to be entered into path by using :: */
 	   memmove(next,next+1,strlen(next+1)+1); /* shuffle up data [including \0]*/
	   next = strchr(next+1,':');
	}
        if(next)
	    *next = '\0';

	while (path!= NULL) {

#ifdef USE_CORBA	
	  /* 19/03/99 johnt - if it is a corba path - look there */
          if( !strncmp( CORBATAG,path,strlen(CORBATAG))){
	    if(corba_found(wholePath,path+strlen(CORBATAG),file)){
	      f = wholePath;
	      break;
	    }
	  }
	  else
#endif
#ifdef USE_BIOLIMS
	  if( !strncmp( BIOLIMS_TAG,path,strlen(BIOLIMS_TAG))){
	    if(biolims_found(wholePath,path+strlen(BIOLIMS_TAG),file)){
	      f = wholePath;
	      break;
	    }
	  }
	  else
#endif    
	  {
	    (void) strcpy(wholePath,path);
	    (void) strcat(wholePath,"/");
	    (void) strcat(wholePath,file);
	    if (found(wholePath)) {
	      f = wholePath;
	      break;
	    }
	  }
	  path = next;
	  if( path ){
	      path++;
	      next = strchr(path,':');
	      while( next && (*(next+1) == ':' )){
    		 /* 26/03/99 johnt - allow : to be entered into path by using :: */
 		 memmove(next,next+1,strlen(next+1)+1); /* shuffle up data */
		 next = strchr(next+1,':');
	      }
	      if(next)
		*next='\0';
	  }
	}
	free(paths);
    }

    return f;
}
