%option c++ prefix="IPV2" yylineno

%{

#if !defined(SUN4)
#include <string.h>
#endif

#include <util/keyval/ipv2_scan.h>
#include <util/keyval/ipv2_parse.h>

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
                    cerr << "IPV2: {string} rule: malloc failed" << endl;
                    abort();
                    }
                  strcpy(yylval.str,yytext);
                  return(T_STRING);
                  }
{qstring}       { yylval.str = (char *)malloc(strlen(yytext));
                  if (!yylval.str) {
                    cerr << "IPV2: {qstring} rule: malloc failed" << endl;
                    abort();
                    }
                  strcpy(yylval.str,&yytext[1]);
                  yylval.str[strlen(yylval.str)-1] = '\0';
                  return(T_STRING);
                  }
[ \t]+          ; 
%.*$            ;
"\n"            ;
"("             { return(T_KEYWORD_LEFT); }
")"             { return(T_KEYWORD_RIGHT); }
"["             { return(T_ARRAY_LEFT); }
"]"             { return(T_ARRAY_RIGHT); }
"{"             { return(T_TABLE_LEFT); }
"}"             { return(T_TABLE_RIGHT); }
[,<>;=:\$]      { return((int) yytext[0]); }
.               { cerr<<"IPV2: Illegal character: \""<<yytext[0]<<"\"\n"; }
%%

int
IPV2wrap()
{
  return 1;
}
