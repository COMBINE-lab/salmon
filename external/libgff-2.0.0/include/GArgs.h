/*
GArgs is a quick'n'dirty object oriented replacement for the standard 
   getopts library available on many unix platforms;
   it accepts the regular single letter, single-dash style options 
     -<letter>[ ][<value>] 
   but also attr=value style options:
     <optname>=<value>
     or
     --<optname>[=]<value>
*/

#ifndef G_ARGS_DEFINED
#define G_ARGS_DEFINED

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

struct GArgsDef {
  const char* longopt;
  char opt; //equivalent one-char option, if any
  bool req_value; //true if the string that follows must be a value
  int code; //an enum code to be associated with this option
};

class GArgs {
   //structure for parsing arguments format definition
   struct fmtdef {
     char* longopt;
     char opt; //equivalent one-char option, if any
     bool req_value; //true if the string that follows must be a value
     int code; //an enum code to be associated with this option
     };
   int fmtcount;
   fmtdef* fmt; //this will store format definition after parsing it
   struct argdata {
     char*  opt; // this is NULL for non-dashed arguments
                 // a single character for single dash style arguments
                //  a string for ARG=VALUE or --long_option style arguments
     char* value; // is NULL for switches (dashed flags)
     int fmti; //index in fmt table
     //int code; // if GArgsDef[] constructor was used, for getOpt
     };
   int _argc;            
   char** _argv; //the original main() values
   argdata* args; //arguments table after parsing it
   int count; //total count of elements in 'args' array
   int nonOptCount; //count of non-dashed, non= arguments
   int nonOptPos; //current position for nonOpt arguments iterator
   int optPos; //current position for options iterator
   int errarg; //argv error position after parsing
   bool err_valmissing; //if the error is strictly about missing value for errarg option
   int parseArgs(bool nodigitopts=false);
   //parsing helper functions
   int validOpt(int c);  
   int validShortOpt(char o);  
   int validLongOpt(char* o, char* to);
 public:
 
   GArgs(int argc, char* argv[], const char* format, bool nodigitopts=false);
   /* format can be:
       <string>{;|=} e.g. disable-test;PID=S= for --disable-test PID=50 (or --PID 50) S=3.5 etc.
       <letter>[:]  e.g. p:hT  for -p testing (or -ptesting) -h -T
   This means that the long options, if present, should be given at the beginning
   of the format string, before the single-dash, single-char options
   */
   GArgs(int argc, char* argv[], const GArgsDef fmtrecs[], bool nodigitopts=false);
   
   ~GArgs();
   int isError(); // returns the offending argv position or 0 if no error
   int getCount() { return count; } //total number of arguments given
   int getFmtCount() { return fmtcount; } //total number of option definitions
   int getNonOptCount() { return nonOptCount; } //total number of non-option arguments
   char* getOpt(const char* o); /* retrieve the value for option o
                   returns 
                       NULL    if option not given at all
                     !=NULL    if boolean option was given
                     opt's value if value option was given
                     */
   char* getOpt(const char o);
   char* getOpt(int c); //retrieve value by enum code
   char* getOptName(int c); //retrieve name of by enum code
   int startOpt(); //init iteration through option arguments
       // returns number of option args
       
   char* nextOpt(); //get next option argument's string
   int nextCode(); //get next option argument's code

   int startNonOpt(void); //init iteration through non-option arguments
             // returns the number of non-option arguments
   void printError(FILE* fout, const char* usage=NULL,
                      bool exitProgram=false);
   void printError(const char* usage=NULL, bool exitProgram=false);
   void printCmdLine(FILE* fout);
   char* nextNonOpt(); //get the next non-option argument
};

#endif
