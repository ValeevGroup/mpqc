%option c++ prefix="MkClasses" yylineno

%{

#define in_mkclasses_scanner
#include "mkclasses.h"
#include "parse.h"

%}
integer [0-9]+
string  [A-Za-z_\./][A-Za-z0-9_\./]*
qstring \"[^"\n]+\"
%%
"##"                    { return T_CONCAT; }
"classes"               { return T_CLASSES; }
"files"                 { return T_FILES; }
"unique_name"           { return T_UNIQUE_NAME; }
"header_includes"       { return T_HEADER_INCLUDES; }
"source_includes"       { return T_SOURCE_INCLUDES; }
"member_data"           { return T_MEMBER_DATA; }
"name"                  { return T_NAME; }
"type"                  { return T_TYPE; }
"long_name"             { return T_LONG_NAME; }
"parents"               { return T_PARENTS; }
"version"               { return T_VERSION; }
"description"           { return T_DESCRIPTION; }
"options"               { return T_OPTIONS; }
"castdown"              { return T_CASTDOWN; }
"require_castdown"      { return T_REQUIRE_CASTDOWN; }
"static_class_desc"     { return T_STATIC_CLASS_DESC; }
"ctors"                 { return T_CTORS; }
"yes"                   { return T_YES; }
"no"                    { return T_NO; }
"true"                  { return T_TRUE; }
"false"                 { return T_FALSE; }
"statein"               { return T_STATEIN; }
"keyval"                { return T_KEYVAL; }
"void"                  { return T_VOID; }
{integer}       { yylval.i = atoi(yytext);
                  return T_INTEGER;
                  }
{string}        { int strlenyytext = strlen(yytext);
                  yylval.charstr = new char[strlenyytext+1];
                  if (!yylval.charstr) {
                    cerr << "mkclasses: {string} rule: malloc failed" << endl;
                    abort();
                    }
                  strcpy(yylval.charstr,yytext);
                  return T_STRING;
                  }
{qstring}       { yylval.charstr = new char[strlen(yytext)];
                  if (!yylval.charstr) {
                    cerr << "mkclasses: {qstring} rule: malloc failed" << endl;
                    abort();
                    }
                  strcpy(yylval.charstr,&yytext[1]);
                  yylval.charstr[strlen(yylval.charstr)-1] = '\0';
                  return(T_STRING);
                  }
[ \t]+          ; 
%.*$            ;
"\n"            ;
[\(\)\[\]=:\$]  { return((int) yytext[0]); }
.               { cerr<<"mkclasses: Illegal character: \""<<yytext[0]<<"\"\n"; }
%%
