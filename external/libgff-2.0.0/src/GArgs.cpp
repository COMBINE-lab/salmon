#include "GBase.h"
#include "GArgs.h"
#include <ctype.h>

GArgs::GArgs(int argc, char* argv[], const char* format, bool nodigitopts) {
   /* format can be:
      <string>{;|=} e.g. disable-test;PID=S= for --disable-test PID=50 (or --PID 50) S=3.5 etc.
      <letter>[:]  e.g. p:hT  for -p testing (or -ptesting) -h -T
   */
const char* fstr=format;
fmtcount=0;
count=0;
nonOptCount=0;
nonOptPos=0;
optPos=0;
errarg=0;
err_valmissing=false;
args=NULL;
fmt=NULL;
_argc=argc;
_argv=argv;
int fmtlen=strlen(format);
//---- first parse the format string
while (fstr-format < fmtlen ) {
  int l=strcspn(fstr, ";=:");
  if (fstr[l]==0) { //end of string reached
      //all previous chars are just switches:
       GREALLOC(fmt, (fmtcount+l)*sizeof(fmtdef));
       //store each switch
       for (int i=0; i<l;i++) { 
         fmt[fmtcount+i].longopt=NULL;
         fmt[fmtcount+i].opt=fstr[i];
         fmt[fmtcount+i].req_value = false;
         fmt[fmtcount+i].code=fmtcount+i+1;
         }
       fmtcount+=l;
       break;
     }
   else {
     if (fstr[l]==':') {
         //fstr[l-1] is an argument, but all the previous are just switches
         GREALLOC(fmt, (fmtcount+l)*sizeof(fmtdef));
         //store each switch AND the option
         for (int i=0; i<l;i++) { 
           fmt[fmtcount+i].longopt=NULL; //one char length
           fmt[fmtcount+i].opt=fstr[i];
           fmt[fmtcount+i].req_value = (i==l-1);
           fmt[fmtcount+i].code=fmtcount+i+1;
           }
         fmtcount+=l;
         }
      else { // fstr[l]=='=' or ';' 
         GREALLOC(fmt, (fmtcount+1)*sizeof(fmtdef));
         fmt[fmtcount].longopt=Gstrdup(fstr, fstr+l-1);
         fmt[fmtcount].opt=0;
         fmt[fmtcount].req_value=(fstr[l]=='=');
         fmt[fmtcount].code=fmtcount+1;
         fmtcount++;
         }
     fstr+=l+1;
     }
  }
 //-- now parse the arguments based on the given format specification
 parseArgs(nodigitopts);
 }

void dbg_dequote(char* &p) {
   int alen=strlen(p);
   if (alen>1 && p[0]=='\'' && p[alen-1]=='\'') {
	 p++;
	 p[alen-2 ]='\0';
   }
}

int GArgs::parseArgs(bool nodigitopts) {
  int p=1; //skip program name
  int f=0;
  while (p<_argc) {
   // silly patch for annnoying MacOS gdb/eclipse issues:
   #if defined(__APPLE__) && defined(DEBUG)
	   dbg_dequote(_argv[p]);
   #endif
   //--
   if (_argv[p][0]=='-' && (_argv[p][1]==0 || _argv[p][1]!='-')) { 
     //single-dash argument
     int cpos=1;
     char c=_argv[p][cpos];
     if (c==0 || (nodigitopts && isdigit(c)) || 
            (c=='.' && isdigit(_argv[p][cpos+1]))) { 
        //special case: plain argument '-' or just a negative number
        GREALLOC(args, (count+1)*sizeof(argdata));
        args[count].opt=NULL;
        args[count].fmti=-1;
        if (c==0) {
          GCALLOC(args[count].value, 2);
          args[count].value[0]='-';
          }
         else { //negative number given
          args[count].value=Gstrdup(_argv[p]);
          }
        count++;
        nonOptCount++;
        }
      else { //single-dash argument or switch
       COLLAPSED:
        if ((f=validShortOpt(c))>=0) {
          GREALLOC(args, (count+1)*sizeof(argdata));
          GCALLOC(args[count].opt, 2);
          args[count].opt[0]=c;
          args[count].fmti=f;
          if (!fmt[f].req_value) {//switch type
            GCALLOC(args[count].value,1);//so getOpt() functions would not return NULL
            count++;
            // only switches can be grouped with some other switches or options
            if (_argv[p][cpos+1]!='\0') {
               cpos++;
               c=_argv[p][cpos];
               goto COLLAPSED;
               }
            }
           else {
              //single-dash argument followed by a value
            if (_argv[p][cpos+1]=='\0') {
                if (p+1<_argc && _argv[p+1][0]!=0) { //value is the whole next argument
                   p++;
				#if defined(__APPLE__) && defined(DEBUG)
				   dbg_dequote(_argv[p]);
				#endif
                   args[count].value=Gstrdup(_argv[p]);
                }
                  else {
                   errarg=p;
                   err_valmissing=true;
                   return errarg;
                   }
                }
               else { //value immediately follows the dash-option
                args[count].value=Gstrdup(_argv[p]+cpos+1);
                }
            count++;
            }
          } //was validShortOpt
         else { //option not found in format definition!
           errarg=p;
           return errarg;
           }
        }
     } //-single-dash
   else {//not a single-dash argument
     char* ap=_argv[p];
     bool is_longopt=false;
     if (*ap=='-' && ap[1]=='-') { //double-dash option
        is_longopt=true;
        ap+=2;
        }
     char* e=strchr(ap+1,'=');
     while (e!=NULL && *(e-1)=='\\') e=strchr(e,'=');
     if (e==NULL && is_longopt) {
        e=ap;
        while (*e!=0 && *e!=' ') e++;
        //e will be on eos or next space
        }
     if (e!=NULL && e>ap) {
       //this must be a long option
       //e is on eos, space or '='
       if ((f=validLongOpt(ap,e-1))>=0) {
            GREALLOC(args, (count+1)*sizeof(argdata));
            args[count].opt=Gstrdup(ap,e-1);
            args[count].fmti=f;
            if (fmt[f].req_value) {
               if (*e==0) {
                   //value is the next argument
                   if (p+1<_argc && _argv[p+1][0]!=0) {
                      p++;
					#if defined(__APPLE__) && defined(DEBUG)
                      dbg_dequote(_argv[p]);
   	   	   	   	   	#endif
                      args[count].value=Gstrdup(_argv[p]);
                      }
                    else {
                      errarg=p;
                      err_valmissing=true;
                      return errarg;
                      }
                   }
                else { //value is in the same argument
                   //while (*e!=0 && (*e==' ' || *e=='=')) e++;
                   if (*e=='=') e++;
                   if (*e==0) {
                      errarg=p;
                      err_valmissing=true;
                      return errarg;
                      }
                   args[count].value=Gstrdup(e);
                   }
               } //value required
              else { //no value expected
               GCALLOC(args[count].value,1); //do not return NULL
               }
            count++;
            }
          else { //error - this long argument not recognized
           errarg=p;
           return errarg;
           }
        }
      else { //just a plain non-option argument
       if (e==ap) { //i.e. just "--"
          errarg=p;
          return errarg;
          }
       GREALLOC(args, (count+1)*sizeof(argdata));
       args[count].opt=NULL; //it's not an option
       args[count].value=Gstrdup(_argv[p]);
       args[count].fmti=-1;
       count++;
       nonOptCount++;
       }
     }
   p++;//check next arg string
   } //while arguments
 return errarg;
}

void GArgs::printError(FILE* fout, const char* usage, bool exitProgram) {
 if (errarg==0) return;
 if (usage) fprintf(fout, "%s\n", usage);
 if (err_valmissing) 
     fprintf(fout, "Error: value required for option '%s'\n", _argv[errarg]);
    else 
     fprintf(fout, "Error: invalid argument '%s'\n", _argv[errarg]);
 if (exitProgram)
     exit(1);
}

void GArgs::printError(const char* usage, bool exitProgram) {
 printError(stderr, usage, exitProgram);
}

void GArgs::printCmdLine(FILE* fout) {
 if (_argv==NULL) return;
 for (int i=0;i<_argc;i++) {
   fprintf(fout, "%s%c", _argv[i], (i==_argc-1)?'\n':' ');
   }
}

GArgs::GArgs(int argc, char* argv[], const GArgsDef fmtrecs[], bool nodigitopts) {
 fmtcount=0;
 count=0;
 nonOptCount=0;
 nonOptPos=0;
 optPos=0;
 errarg=0;
 err_valmissing=false;
 args=NULL;
 fmt=NULL;
 _argc=argc;
 _argv=argv;
 if (fmtrecs==NULL) return;
 
 const GArgsDef* frec=fmtrecs;
 while ((frec->longopt || frec->opt) && fmtcount<255) {
     fmtcount++;
     frec=&(fmtrecs[fmtcount]);
     }
 GCALLOC(fmt, fmtcount*sizeof(fmtdef));
 for (int i=0;i<fmtcount;i++) {
   fmt[i].longopt=Gstrdup(fmtrecs[i].longopt); //do we need to use Gstrdup here?
   fmt[i].opt=fmtrecs[i].opt;
   fmt[i].req_value=fmtrecs[i].req_value;
   fmt[i].code=fmtrecs[i].code;
   }
 parseArgs(nodigitopts);
}


GArgs::~GArgs() {
 int i;
 for (i=0; i<fmtcount; i++)
    GFREE(fmt[i].longopt);
 GFREE(fmt);
 for (i=0; i<count; i++) {
  GFREE(args[i].opt);
  GFREE(args[i].value);
  }
 GFREE(args);
}

int GArgs::validShortOpt(char o) {
 for (int i=0; i<fmtcount; i++) 
  if (fmt[i].opt==o) return i;
 return -1; 
}

int GArgs::validLongOpt(char* o, char* to) {
 char* pstr=Gstrdup(o,to);
 for (int i=0; i<fmtcount; i++) {
  if (fmt[i].longopt && strcmp(fmt[i].longopt, pstr)==0) {
       GFREE(pstr);
       return i;
       }
  }
 GFREE(pstr); 
 return -1;
}

int GArgs::validOpt(int code) {
 for (int i=0; i<fmtcount; i++) 
   if (fmt[i].code==code) return i;
 return -1;
}


int GArgs::isError() { // returns the offending argv position or 0 if no error
 return errarg;
 }

char* GArgs::getOpt(const char* o) { /* retrieve the value for option o
                   returns 
                       NULL    if option not given at all
                     !=NULL    if boolean option was given
                     opt.value if value option was given
                     */
 for (int i=0; i<count; i++) 
  if (args[i].opt!=NULL && strcmp(args[i].opt, o)==0) 
           return args[i].value;
 return NULL;
}

char* GArgs::getOpt(const char o) {
 for (int i=0; i<count; i++) 
  if (args[i].opt!=NULL && args[i].opt[0]==o && args[i].opt[1]=='\0') 
      return args[i].value;
 return NULL;
}

char* GArgs::getOpt(int c) {
 for (int i=0; i<count; i++) 
  if (args[i].fmti>=0 && fmt[args[i].fmti].code==c)
      return args[i].value;
 return NULL;
}

char* GArgs::getOptName(int c) {
 for (int i=0; i<count; i++) 
  if (args[i].fmti>=0 && fmt[args[i].fmti].code==c)
      return args[i].opt;
 return NULL;
}


int GArgs::startNonOpt(){ //reset iteration through non-option arguments
   //returns the number of non-option arguments
nonOptPos=0;
return nonOptCount;   
}
   
   
char* GArgs::nextNonOpt() { //get the next non-dashed argument
               //or NULL if no more 
for (int i=nonOptPos;i<count;i++)
 if (args[i].opt==NULL) {
      nonOptPos=i+1;
      return args[i].value;
      }
return NULL;
}

int GArgs::startOpt(){ //reset iteration through option arguments
   //returns the number of option arguments
optPos=0;
return count-nonOptCount;
}
   
   
char* GArgs::nextOpt() { //get the next non-dashed argument
               //or NULL if no more 
for (int i=optPos;i<count;i++)
 if (args[i].opt!=NULL) {
      optPos=i+1;
      return args[i].opt;
      }
return NULL;
}

int GArgs::nextCode() { //get the next non-dashed argument
               //or NULL if no more 
for (int i=optPos;i<count;i++)
 if (args[i].opt!=NULL && args[i].fmti>=0) {
      optPos=i+1;
      return fmt[args[i].fmti].code;
      }
return 0; //must make sure that codes are > 0 for this to work properly
}

