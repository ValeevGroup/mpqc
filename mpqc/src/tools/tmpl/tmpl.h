
/* $Log$
 * Revision 1.1  1993/12/29 12:53:59  etseidl
 * Initial revision
 *
 * Revision 1.2  1992/03/30  23:12:35  seidl
 * merge in sandia changes
 *
 * Revision 1.4  1991/12/11  03:02:21  cljanss
 * Really move includes from here to main.
 *
 * Revision 1.3  1991/12/04  23:33:15  cljanss
 * moved includes from tmpl.h to main.c
 *
 * Revision 1.2  1991/12/02  17:41:27  cljanss
 * converted to the DEC
 *
 * Revision 1.1  91/06/17  01:29:36  janssen
 * Initial revision
 *  */

#ifndef _TMPL_H
#define _TMPL_H

#ifdef DEC
#define NO_CONST
#endif

#ifdef NCUBE_V2
#define NO_CONST
#endif

#ifdef NO_CONST
#define CONST
#else
#define CONST const
#endif

#define VOID void
#define VOID_PTR void *

#ifdef PROTO_LIBRARY

VOID_PTR malloc(size_t);
VOID_PTR realloc(VOID_PTR,size_t);
VOID free(VOID_PTR);

FILE *fopen(CONST char *, CONST char *);
int   fclose(FILE *);
int fprintf(FILE *, CONST char *, ... );
int printf(CONST char *, ... );
int sprintf(char *, CONST char *, ... );

VOID perror(CONST char *);

VOID exit(int);

int strcmp(CONST char *, CONST char *);
int strlen(CONST char *);
char *strcpy(char *, CONST char *);
char *strcat(char *, CONST char *);

/* These are needed for getc and putc (to avoid sun3 gcc warnings). */
int _filbuf(), _flsbuf();

int yylex();

int yyparse();

/* This is for gcc. */
void abort( void );

#endif /* PROTO_LIBRARY */

#endif /* _TMPL_H */

