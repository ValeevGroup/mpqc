%{

/* $Log$
 * Revision 1.2  1996/05/02 23:17:04  cljanss
 * SP2 update
 *
 * Revision 1.1.1.1  1993/12/29 12:54:00  etseidl
 * SC source tree 0.1
 *
 * Revision 1.2  1992/03/30  23:13:09  seidl
 * merge in sandia changes
 *
 * Revision 1.4  1992/02/04  01:55:59  cljanss
 * changed NO_TEMPLATES to _STDC__
 *
 * Revision 1.3  1991/12/04  23:33:34  cljanss
 * removed stdio.h include to eliminate a makedepend error message
 *
 * Revision 1.2  1991/12/04  17:35:21  cljanss
 * cpp derivatives around function arg declarations are allowed in a
 * primitive way (to allow varargs to work with either ansi or nonansi
 * compilers)
 *
 * Revision 1.1  91/06/17  01:29:36  janssen
 * Initial revision
 *  */
static char rcsid[] = "$Id$";

int ce_depth;
int global,vargs;
int might_need_comma = 0;
int prefer_newline = 1;
int foundname;
int noargs;
extern FILE *tmplglobal;
extern FILE *tmpllocal;

%}

ws         [ \t\n]*
identifier [a-zA-Z_][a-zA-Z0-9_]*
noarg      {ws}")"
arg        (({ws}{identifier}{ws}",")|({ws}{identifier}{ws}")")|({ws}")"))

%START S_FUNCTION_DEC S_FUNCTION_ARG S_FUNCTION_PAR CONSTANT_EXPRESSION

%%
^"GLOBAL_FUNCTION"			{ BEGIN S_FUNCTION_DEC;
							  yyout = tmplglobal;
							  fprintf(yyout,"extern\n#ifdef __cplusplus\n");
							  fprintf(yyout,"\"C\"\n#endif\n");
							  global = 1; vargs = 0; }
^"GLOBAL_VA_FUNCTION"		{ BEGIN S_FUNCTION_DEC;
							  yyout = tmplglobal;
							  fprintf(yyout,"extern\n#ifdef __cplusplus\n");
							  fprintf(yyout,"\"C\"\n#endif\n");
							  global = 1; vargs = 1; }
^"LOCAL_FUNCTION"			{ BEGIN S_FUNCTION_DEC;
							  yyout = tmpllocal;
							  fprintf(yyout,"LOCAL_STORAGE_CLASS");
							  global = 0; vargs = 0; }
^"LOCAL_VA_FUNCTION"		{ BEGIN S_FUNCTION_DEC;
							  yyout = tmpllocal;
							  fprintf(yyout,"LOCAL_STORAGE_CLASS");
							  global = 0; vargs = 1; }
<S_FUNCTION_DEC>"#"" "*"ifndef".*$		{}
<S_FUNCTION_DEC>"("/{noarg}	{ BEGIN S_FUNCTION_PAR; noargs = 1; ECHO; }
<S_FUNCTION_DEC>"("/{arg}	{ BEGIN S_FUNCTION_PAR; noargs = 0; ECHO; }
<S_FUNCTION_DEC>.			{ ECHO; }
<S_FUNCTION_PAR>","" "*"..."" "*")" { BEGIN S_FUNCTION_ARG;
							  if (!vargs) {
                                fprintf(stderr,"tmpl: inconsistent use");
							    fprintf(stderr," of \"VA_FUNCTION\"");
								fprintf(stderr," and \"...\"\n");
                                exit(1);
                                }
							  fprintf(yyout,"\n#if defined(__STDC__) || defined(__cplusplus)\n");
							  }
<S_FUNCTION_PAR>")"			{ BEGIN S_FUNCTION_ARG;
							  fprintf(yyout,"\n#if defined(__STDC__) || defined(__cplusplus)\n");
							  }
<S_FUNCTION_PAR>.			;
<S_FUNCTION_ARG>"["			{ BEGIN CONSTANT_EXPRESSION;
							  ce_depth=1;
							  ECHO;
							  }
<CONSTANT_EXPRESSION>"["	{ ce_depth++; ECHO; }
<CONSTANT_EXPRESSION>"]"	{ if (!(--ce_depth)) {
							    BEGIN S_FUNCTION_ARG;
							    }
							  ECHO;
							  }
<S_FUNCTION_ARG>";"			{ might_need_comma = 1; }
<S_FUNCTION_ARG>"#"" "*"else"		{ BEGIN 0;
							  if (might_need_comma&&vargs) fprintf(yyout,",");
							  might_need_comma = 0;
							  if (prefer_newline) {
							    fprintf(yyout,"\n");
							    prefer_newline = 0;
							    }
							  if (vargs) fprintf(yyout,"...");
							  else if (noargs) fprintf(yyout,"void");
							  fprintf(yyout,"\n#endif\n");
							  fprintf(yyout,");\n"); }
<S_FUNCTION_ARG>"{"			{ BEGIN 0;
							  if (might_need_comma&&vargs) fprintf(yyout,",");
							  might_need_comma = 0;
							  if (prefer_newline) {
							    fprintf(yyout,"\n");
							    prefer_newline = 0;
							    }
							  if (vargs) fprintf(yyout,"...");
							  else if (noargs) fprintf(yyout,"void");
							  fprintf(yyout,"\n#endif\n");
							  fprintf(yyout,");\n"); }
<S_FUNCTION_ARG>.			{ if (might_need_comma) {
							    fprintf(yyout,",");
							    might_need_comma = 0;
							    if (prefer_newline) {
                                  fprintf(yyout,"\n");
							      prefer_newline = 0;
							      }
							    }
							  ECHO; }
<CONSTANT_EXPRESSION>.		{ ECHO; }
.							;
<S_FUNCTION_DEC>"\n"		{ ECHO; }
<S_FUNCTION_ARG>"\n"		{ prefer_newline = 1; }
<S_FUNCTION_DEC>[ \t]		{ ECHO; }
<S_FUNCTION_ARG>[ \t]		{ ECHO; }
"\n"						;
%%


int
yywrap()
{
  return 1;
}
