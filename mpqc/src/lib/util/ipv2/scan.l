%{
#include <tmpl.h>
#if !defined(SUN4)
#include <string.h>
#endif
#include "ip_types.h"
#include "ip_global.h"
#include "ip_error.gbl"
#include "scan.gbl"
#include "scan.lcl"
#include "parse.h"
#ifdef __STDC__
#define YY_USE_PROTOS
#endif
#ifndef FLEX
#undef input
#define input() (((yytchar=yysptr>yysbuf?U(*--yysptr):getc(yyin))==10?(advline(),yylineno++,yytchar):yytchar)==EOF?0:(advchar(),yytchar))
#undef unput
#define unput(c) {bckchar();yytchar= (c);if(yytchar=='\n'){yylineno--;bckline();};*yysptr++=yytchar;}
#endif
%}
string  [A-Za-z0-9_\.*+-/]*
qstring \"[^"\n]+\"
%%
{string}		{ int strlenyytext = strlen(yytext);
				  if (strlenyytext==1) {
				    if (yytext[0]=='*') return '*';
				    if (yytext[0]=='/') return '/';
				    if (yytext[0]=='-') return '-';
				    if (yytext[0]=='+') return '+';
				    }
				  yylval.str = (char *)malloc(strlenyytext+1);
				  if (!yylval.str) {
				    perror("{string} rule");
				    ip_error("{string} rule: malloc failed");
				    }
				  strcpy(yylval.str,yytext);
                  if (ip_uppercase) cvs_toupper(yylval.str);
				  return(T_STRING);
				  }
{qstring}		{ yylval.str = (char *)malloc(strlen(yytext));
				  if (!yylval.str) {
				    perror("{qstring} rule");
				    ip_error("{qstring} rule: malloc failed");
				    }
				  strcpy(yylval.str,&yytext[1]);
				  yylval.str[strlen(yylval.str)-1] = '\0';
				  return(T_STRING);
				  }
[ \t\n]+		;
%.*$			;
#               { return(0); }
"("				{ return(T_KEYWORD_LEFT); }
")"				{ return(T_KEYWORD_RIGHT); }
"["				{ return(T_ARRAY_LEFT); }
"]"				{ return(T_ARRAY_RIGHT); }
"{"				{ return(T_TABLE_LEFT); }
"}"				{ return(T_TABLE_RIGHT); }
[,<>;=:\$]		{ return((int) yytext[0]); }
.				{ ip_error("Illegal character"); }
%%

/* Convert a string to uppercase. */
LOCAL_FUNCTION VOID
cvs_toupper(s)
char *s;
{
  for (; *s!='\0'; s++) {
    if (*s>='a' && *s <='z') *s = *s + 'A' - 'a';
    }
  }

#ifndef FLEX

#define MAX_LL 80
static char cur_line[MAX_LL] = {'\0'};
static int eline = 0;
static int lostline = 0;

/* Advance the input line. */
GLOBAL_FUNCTION VOID
advline()
{
  lostline = 0;
  eline = 0;
  cur_line[0] = '\0';
  }

/* Move the input line back. */
GLOBAL_FUNCTION VOID
bckline()
{
  lostline = 1;
  eline = 0;
  cur_line[0] = '\0';
  }

/* Advance the input char. */
GLOBAL_FUNCTION VOID
advchar()
{
  int i;
  /* Shift the line if it is too long. */
  if (eline==MAX_LL-1) {
    eline += -MAX_LL/2;
    for (i=0; i<eline; i++) {
      cur_line[i] = cur_line[i+MAX_LL/2];
      }
    }
  cur_line[eline++] = (char) yytchar;
  cur_line[eline] = '\0';
  }

GLOBAL_FUNCTION VOID
bckchar()
{
  if (eline==0) {
    bckline();
    return;
    }
  eline--;
  cur_line[eline] = '\0';
  }
#endif
  

/* Show position. */
GLOBAL_FUNCTION VOID
showpos()
{
#ifndef FLEX
  int i;

  fprintf(stderr,"the error occurred at roughly line number %d\n",yylineno+1);

  if (lostline) return;

  /* The first charater is a newline--don't print it. */
  fprintf(stderr,"%s\n",&cur_line[1]);
  for (i=0; i<eline-2-yyleng+1; i++) fprintf(stderr,"_");
  for (i=0; i<yyleng; i++) fprintf(stderr,"^");
  fprintf(stderr,"\n");
#else
  printf("flex has been used to generate the scanner so I don't\n");
  printf("know where the error is\n");
#endif
  }

