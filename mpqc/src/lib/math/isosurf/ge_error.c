/****************************************************************************
 * ge_error.c
 * Author Joel Welling
 * Copyright 1988, Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Permission use, copy, and modify this software and its documentation
 * without fee for personal use or use within your organization is hereby
 * granted, provided that the above copyright notice is preserved in all
 * copies and that that copyright and this permission notice appear in
 * supporting documentation.  Permission to redistribute this software to
 * other organizations or individuals is not granted;  that must be
 * negotiated with the PSC.  Neither the PSC nor Carnegie Mellon
 * University make any representations about the suitability of this
 * software for any purpose.  It is provided "as is" without express or
 * implied warranty.
 *****************************************************************************/

/*
This module provides error handling.
*/

/* If USE_VARARGS below is set, varargs.h-type parameter list handling
 * will be done.  Otherwise, a method which 'usually' works is used.
 */
#define USE_VARARGS 1

#ifdef VMS
#include stdio
#else
#include <stdio.h>
#include <string.h>
#endif
#ifdef USE_VARARGS
#include <varargs.h>
#endif

/* boolean values */
#define FALSE 0
#define TRUE 1
typedef int boolean;

/* Storage for name of executable */
#define maxname 20
static char name[maxname+1];

#ifdef USE_VARARGS
static va_list ap;
static char *p, *sval;
static int ival;
static double dval;
#else
/* Error message buffer */
#define maxstring 128
static char messagebuf[maxstring+maxname+3];
#endif

/* The following macro simulates the action of fprintf in an USE_VARARGS
 * environment.  Unfortunately it currently doesn't implement the full
 * functionality of fprintf. */
#define parsing_macro {\
  va_start( ap ); \
  for (p = msg; *p; p++) { \
    if (*p != '%') { \
      putc( *p, stderr ); \
      continue; \
    } \
    switch (*++p) { \
    case 'd': \
      ival= va_arg(ap,int); \
      fprintf(stderr,"%d",ival); \
      break; \
    case 'f': \
      dval= va_arg(ap,double); \
      fprintf(stderr,"%f",dval); \
      break; \
    case 's': \
      for (sval= va_arg(ap,char*); *sval; sval++) putc( *sval, stderr ); \
      break; \
    case 'c': \
      putc( va_arg(ap,char), stderr ); \
      break; \
    default: \
      putc( *p, stderr ); \
      break; \
    } \
  } \
  va_end(ap); \
} /* end of parsing macro */

/* Debugging on if the following variable is 'TRUE'. */
static boolean debugmode= FALSE;

void ger_init(exename)
char *exename;
/*  
This routine initializes the error handling module.  The parameter
'exename' should be passed the name of the executable, for inclusion in
error messages.
*/
{
        (void) strncpy(name,exename,maxname);
        name[maxname]= '\0';
}

void ger_toggledebug()
/*  This routine toggles debugging on and off. */
{
        if (debugmode) debugmode= FALSE;
        else debugmode= TRUE;
}

#ifdef USE_VARARGS
void ger_error(msg,va_alist)
char *msg;
va_dcl
/* This routine prints an error message.  */
{
  (void) fputs(name, stderr);
  (void) fputs(": ", stderr);
  parsing_macro;
  putc( '\n', stderr );
}

void ger_debug(msg, va_alist)
char *msg;
va_dcl
/* This routine prints a debugging message, if debugging is on. */
{
  if (debugmode) {
    (void) fputs(name, stderr);
    (void) fputs(": ", stderr);
    parsing_macro;
    putc( '\n', stderr );
  }
}

void ger_fatal(msg,va_alist)
char *msg;
va_dcl
/* This routine prints an error message and exits. */
{
  (void) fputs(name, stderr);
  (void) fputs(": ", stderr);
  parsing_macro;
  putc( '\n', stderr );
  (void) exit(2);
}

#else /* USE_VARARGS not set */
void ger_debug(msg,p1,p2,p3,p4,p5)
char *msg;
int p1,p2,p3,p4,p5;
/* This routine prints a debugging message, if debugging is on. */
{
        if (debugmode)
                {
                (void) strcpy(messagebuf,name);
                (void) strncat(messagebuf,": ",2);
                (void) strncat(messagebuf,msg,maxstring);
                (void) strcat(messagebuf,"\n");

                (void) fprintf(stderr,messagebuf,p1,p2,p3,p4,p5);
                }
}

void ger_error(msg,p1,p2,p3,p4,p5)
char *msg;
int p1,p2,p3,p4,p5;
/* This routine prints an error message */
{
        (void) strcpy(messagebuf,name);
        (void) strncat(messagebuf,": ",2);
        (void) strncat(messagebuf,msg,maxstring);
        (void) strcat(messagebuf,"\n");

        (void) fprintf(stderr,messagebuf,p1,p2,p3,p4,p5);
}

void ger_fatal(msg,p1,p2,p3,p4,p5)
char *msg;
int p1,p2,p3,p4,p5;
/* This routine prints an error message and exits. */
{
        ger_error(msg,p1,p2,p3,p4,p5);
        (void) exit(2);
}
#endif /* section not using varargs */






