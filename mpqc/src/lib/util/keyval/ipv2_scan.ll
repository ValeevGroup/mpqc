%{
#define yyin IPV2::yin
#define yyout IPV2::yout
#undef YY_DECL
#define YY_DECL IPV2::ylex()
#if !defined(SUN4)
#include <string.h>
#endif
#include "ipv2.h"
#include "ipv2_parse.h"
#define YY_USE_PROTOS
#define yyrestart IPV2::yrestart
#ifndef OLD_FLEX_SCANNER
#  define yy_get_next_buffer IPV2::y_get_next_buffer
#  define yyunput IPV2::yunput
#  define yyinput IPV2::yinput
   FILE* IPV2::yin=0;
   FILE* IPV2::yout=0;
#endif
%}
string  [A-Za-z0-9_\.*+-/]*
qstring \"[^"\n]+\"
%%
{string}        { int strlenyytext = strlen(yytext);
                  if (strlenyytext==1) {
                    if (yytext[0]=='*') return '*';
                    if (yytext[0]=='/') return '/';
                    if (yytext[0]=='-') return '-';
                    if (yytext[0]=='+') return '+';
                    }
                  yylval.str = (char *)malloc(strlenyytext+1);
                  if (!yylval.str) {
                    perror("{string} rule");
                    error("{string} rule: malloc failed");
                    }
                  strcpy(yylval.str,yytext);
                  if (ip_uppercase) cvs_toupper(yylval.str);
                  return(T_STRING);
                  }
{qstring}       { yylval.str = (char *)malloc(strlen(yytext));
                  if (!yylval.str) {
                    perror("{qstring} rule");
                    error("{qstring} rule: malloc failed");
                    }
                  strcpy(yylval.str,&yytext[1]);
                  yylval.str[strlen(yylval.str)-1] = '\0';
                  return(T_STRING);
                  }
[ \t]+          ; 
"\n"            lineno++;
%.*$            ;
"("             { return(T_KEYWORD_LEFT); }
")"             { return(T_KEYWORD_RIGHT); }
"["             { return(T_ARRAY_LEFT); }
"]"             { return(T_ARRAY_RIGHT); }
"{"             { return(T_TABLE_LEFT); }
"}"             { return(T_TABLE_RIGHT); }
[,<>;=:\$]      { return((int) yytext[0]); }
.               { error("Illegal character"); }
%%

/* Convert a string to uppercase. */
void
IPV2::cvs_toupper(char*s)
{
  for (; *s!='\0'; s++) {
    if (*s>='a' && *s <='z') *s = *s + 'A' - 'a';
    }
  }

/* Show position. */
void
IPV2::showpos()
{
  printf("error occurred at line number %d\n",lineno);
  }

