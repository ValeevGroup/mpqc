
/* $Log$
 * Revision 1.1  1993/12/29 12:53:59  etseidl
 * Initial revision
 *
 * Revision 1.3  1993/04/28  18:31:08  jannsen
 * generated files are placed in the current directory
 *
 * Revision 1.2  1992/06/17  23:11:00  jannsen
 * all tmpl generated includes now automatically include <tmpl.h> if
 * it hasn't been included already
 *
 * Revision 1.1.1.1  1992/03/17  18:09:52  seidl
 * template generator 2.0
 *
 * Revision 1.1  1992/03/17  18:09:51  seidl
 * Initial revision
 *
 * Revision 1.2  1992/01/15  12:30:31  seidl
 * add clj's patches
 *
 * Revision 1.4  1991/12/04  23:32:55  cljanss
 * moved includes from tmpl.h to main.c
 *
 * Revision 1.3  1991/12/02  17:41:03  cljanss
 * shortened name lengths
 *
 * Revision 1.2  91/08/29  16:06:36  cljanss
 * writes out the name of the source file
 * 
 * Revision 1.1  1991/06/17  01:29:36  janssen
 * Initial revision
 * */
static char rcsid[] = "$Id$";

#include <stdio.h>
#include <ctype.h>
#include <sys/types.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "tmpl.h"

extern FILE *yyin;
FILE *tmplglobal;
FILE *tmpllocal;

main(argc,argv)
int argc;
char *argv[];
{
  char *globaltmplname;
  char *localtmplname;
  char *sourcename;
  const char *globalsuffix="gbl";
  const char *localsuffix="lcl";
  const char *tmpglobaltmplname = "tmpl.tmp.gbl";
  const char *tmplocaltmplname = "tmpl.tmp.lcl";
  char *basesourcename;

  assert(argc == 2);

  sourcename = argv[1];

  /* Open the input file. */
  yyin = fopen(sourcename,"r");
  assert(yyin);

  /* The basesourcename is used to determine the names of the gbl and lcl
   * files. */
  if (!strchr(sourcename,'/')) {
    basesourcename = strdup(sourcename);
    }
  else {
    basesourcename = strrchr(sourcename,'/');
    basesourcename = strdup(&basesourcename[1]);
    }

  /* Obtain the name of the global template file. */
  globaltmplname = (char *) malloc(strlen(basesourcename) + strlen(globalsuffix));
  assert(globaltmplname);
  strcpy(globaltmplname,basesourcename);
  globaltmplname[strlen(globaltmplname)-1] = '\0';
  strcat(globaltmplname,globalsuffix);

  /* Obtain the name of the local template file. */
  localtmplname = (char *) malloc(strlen(basesourcename) + strlen(localsuffix));
  assert(localtmplname);
  strcpy(localtmplname,basesourcename);
  localtmplname[strlen(localtmplname)-1] = '\0';
  strcat(localtmplname,localsuffix);

  /* Open the template files. */
  tmplglobal = fopen(tmpglobaltmplname,"w");
  assert(tmplglobal);
  tmpllocal = fopen(tmplocaltmplname,"w");
  assert(tmpllocal);

  /* Initialize the global header file. */
  fprintf(tmplglobal,"#ifndef _TMPL_H\n");
  fprintf(tmplglobal,"#  include <tmpl.h>\n");
  fprintf(tmplglobal,"#endif\n");
  fprintf(tmplglobal,"#ifndef GLOBAL_FUNCTION\n");
  fprintf(tmplglobal,"#  define GLOBAL_FUNCTION\n");
  fprintf(tmplglobal,"#endif\n");
  fprintf(tmplglobal,"#ifndef GLOBAL_VA_FUNCTION\n");
  fprintf(tmplglobal,"#  define GLOBAL_VA_FUNCTION\n");
  fprintf(tmplglobal,"#endif\n");

  /* Initialize the local header file. */
  fprintf(tmpllocal,"\n");
  fprintf(tmpllocal,"#ifndef _TMPL_H\n");
  fprintf(tmpllocal,"#  include <tmpl.h>\n");
  fprintf(tmpllocal,"#endif\n");
  fprintf(tmpllocal,"#ifdef LOCAL_HEADER_USED\n");
  /* Error is not permitted in some cpp's (DEC for example). */
  fprintf(tmpllocal,"/* #  error included two local headers */");
  fprintf(tmpllocal,"#else\n");
  fprintf(tmpllocal,"#  define LOCAL_HEADER_USED\n");
  fprintf(tmpllocal,"#endif\n");
  fprintf(tmpllocal,"#ifndef LOCAL_STORAGE_CLASS\n");
  fprintf(tmpllocal,"#  define LOCAL_STORAGE_CLASS static\n");
  fprintf(tmpllocal,"#endif\n");
  fprintf(tmpllocal,"#ifndef LOCAL_FUNCTION\n");
  fprintf(tmpllocal,"#  define LOCAL_FUNCTION LOCAL_STORAGE_CLASS\n");
  fprintf(tmpllocal,"#endif\n");
  fprintf(tmpllocal,"#ifndef LOCAL_VA_FUNCTION\n");
  fprintf(tmpllocal,"#  define LOCAL_VA_FUNCTION LOCAL_STORAGE_CLASS\n");
  fprintf(tmpllocal,"#endif\n");

  /* Produce the temporary template files. */
  yylex();

  /* Close things up. */
  fclose(tmplglobal);
  fclose(tmpllocal);

  /* See if the templates have actually changed, if so update them. */
  printf("tmpl %s:",sourcename);
  if (updatetmpl(localtmplname,tmplocaltmplname)) {
    printf(" updated %s",localtmplname);
    }
  if (updatetmpl(globaltmplname,tmpglobaltmplname)) {
    printf(" updated %s",globaltmplname);
    }
  printf("\n");
  return(0);
  }

int
updatetmpl(old,new)
char *old, *new;
{
  if (filediff(old,new)) {
    unlink(old);
    assert(!rename(new,old));
    return(1);
    }
  else {
    assert(!unlink(new));
    return(0);
    }
  }

int
filediff(n1,n2)
char *n1,*n2;
{
  FILE *f1,*f2;
  int ch1, ch2;

  f1 = fopen(n1,"r");
  if (!f1) return(1);
  f2 = fopen(n2,"r");
  assert(f2);

  for (;;) {
    ch1 = getc(f1);
    ch2 = getc(f2);
    if (ch1 != ch2) {
      fclose(f1);
      fclose(f2);
      return(1);
      }
    else if (ch1 == EOF) {
      fclose(f1);
      fclose(f2);
      return(0);
      }
    }
  }

