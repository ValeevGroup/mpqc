/*
 * ipv2_scan.ll
 *
 * Copyright (C) 1996 Limit Point Systems, Inc.
 *
 * Author: Curtis Janssen <cljanss@ca.sandia.gov>
 * Maintainer: LPS
 *
 * This file is part of the SC Toolkit.
 *
 * The SC Toolkit is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * The SC Toolkit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public License
 * along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
 * the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * The U.S. Government is granted a limited license as per AL 91-7.
 */
%option c++ prefix="IPV2" yylineno

%{

#if !defined(SUN4)
#include <string.h>
#endif

#include <util/misc/exenv.h>
#include <util/keyval/ipv2_scan.h>
#include <util/keyval/ipv2_parse.h>

using namespace sc;

#define YY_NO_UNISTD_H
extern "C" int IPV2wrap();

#ifndef yywrap
#  define yywrap IPV2wrap
#endif

%}
string  [A-Za-z0-9_\.*+-/]*
qstring \"[^"\n]+\"
%%
"##"            { return T_CONCAT; }
{string}        { int strlenyytext = strlen(yytext);
                  if (strlenyytext==1) {
                    if (yytext[0]=='*') return '*';
                    if (yytext[0]=='/') return '/';
                    if (yytext[0]=='-') return '-';
                    if (yytext[0]=='+') return '+';
                    }
                  yylval.str = (char *)malloc(strlenyytext+1);
                  if (!yylval.str) {
                    ExEnv::errn() << "IPV2: {string} rule: malloc failed" << endl;
                    abort();
                    }
                  strcpy(yylval.str,yytext);
                  return(T_STRING);
                  }
{qstring}       { yylval.str = (char *)malloc(strlen(yytext));
                  if (!yylval.str) {
                    ExEnv::errn() << "IPV2: {qstring} rule: malloc failed" << endl;
                    abort();
                    }
                  strcpy(yylval.str,&yytext[1]);
                  yylval.str[strlen(yylval.str)-1] = '\0';
                  return(T_QUOTED_STRING);
                  }
[ \t]+          ; 
%.*$            ;
[\n\r\f]        ;
"("             { return(T_KEYWORD_LEFT); }
")"             { return(T_KEYWORD_RIGHT); }
"["             { return(T_ARRAY_LEFT); }
"]"             { return(T_ARRAY_RIGHT); }
"{"             { return(T_TABLE_LEFT); }
"}"             { return(T_TABLE_RIGHT); }
[,<>;=:\$]      { return((int) yytext[0]); }
.               { ExEnv::errn()<<"IPV2: Illegal character: \""<<yytext[0]<<"\"\n"; }
%%

int
IPV2wrap()
{
  return 1;
}
