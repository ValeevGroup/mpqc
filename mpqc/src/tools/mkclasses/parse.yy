
%{
#include <stdlib.h>
#include "mkclasses.h"

#define yyerror MkClasses::yerror
#define yyparse MkClasses::yparse
#define yylex MkClasses::ylex
#define yywrap MkClasses::ywrap
%}

%union {
  int i;
  char *charstr;
  string *str;
  vector<string> *vecstr;
  }

%token T_CLASSES T_FILES
%token T_UNIQUE_NAME T_HEADER_INCLUDES T_SOURCE_INCLUDES
%token T_NAME T_PARENTS T_VERSION T_DESCRIPTION
%token T_LONG_NAME T_TYPE T_MEMBER_DATA
%token T_OPTIONS
%token T_CASTDOWN T_REQUIRE_CASTDOWN T_STATIC_CLASS_DESC
%token T_CTORS
%token T_CONCAT
%token T_VOID T_KEYVAL T_STATEIN
%token T_YES T_NO T_TRUE T_FALSE
%token <charstr> T_STRING
%token <i> T_INTEGER


%type <str> string substring partialstring
%type <vecstr> string_array strings
%type <i> boolean integer

%start input
%%

input:          group_defs
            ;

group_defs:     group_defs group_def
            |
            |   error
            ;

group_def:      classes_group
            |   files_group
        ;

files_group: T_FILES ':' '(' file_defs ')'
        ;

file_defs:     file_defs file_def
            |
            |  error
            ;

file_def:       T_UNIQUE_NAME '=' string
                        {set_unique_name(*$3);
                         delete $3;}
            |   T_HEADER_INCLUDES '=' string_array
                        {set_header_inc(*$3);
                         delete $3;}
            |   T_SOURCE_INCLUDES '=' string_array
                        {set_source_inc(*$3);
                         delete $3;}
        ;

classes_group:  T_CLASSES ':' '[' class_defs ']' 
        ;

class_defs:     class_defs class_def
            |
            |   error
            ;

class_def:      '(' { start_class(); } class_items { end_class(); } ')'
        ;

class_items:    class_items class_item
            |
            |   error
            ;

class_item: T_NAME '=' string
                        { current_class().set_name(*$3);
                          delete $3; }
            |   T_PARENTS '=' string_array
                        { current_class().set_parents(*$3);
                          delete $3; }
            |   T_VERSION '=' integer
                        { current_class().set_version($3); }
            |   T_DESCRIPTION '=' string
                        { current_class().set_description(*$3);
                          delete $3; }
            |   class_options_group
            |   class_ctors_group
            |   member_data_group
        ;

member_data_group: T_MEMBER_DATA ':' '[' member_data ']'
        ;

member_data:    member_data '('    { current_class().start_member_datum(); }
                member_datum_items { current_class().end_member_datum(); }
                ')'
            |
        ;

member_datum_items: member_datum_items member_datum_item
            |
            | error
        ;

member_datum_item:
                T_NAME '=' string
     { current_class().current_member_datum().set_name(*$3);
       delete $3; }
            |   T_TYPE '=' string
     { current_class().current_member_datum().set_type(*$3);
       delete $3; }
            |   T_LONG_NAME '=' string
     { current_class().current_member_datum().set_long_name(*$3);
       delete $3; }
            |   T_DESCRIPTION '=' string
     { current_class().current_member_datum().set_description(*$3);
       delete $3; }
        ;

class_options_group: T_OPTIONS ':' '(' class_options ')'
        ;

class_options: class_options class_option
            |
            | error
        ;

class_option:   T_CASTDOWN '=' boolean
               { current_class().set_castdown($3); }
            |   T_REQUIRE_CASTDOWN '=' boolean
               { current_class().set_require_castdown($3); }
            |   T_STATIC_CLASS_DESC '=' boolean
               { current_class().set_static_class_desc($3); }
        ;

class_ctors_group: T_CTORS ':' '(' class_ctors ')'
        ;

class_ctors: class_ctors class_ctor
            |
            | error
        ;

class_ctor:     T_VOID '=' boolean   
                        { current_class().set_ctor_void($3); }
            |   T_KEYVAL '=' boolean
                        { current_class().set_ctor_keyval($3); }
            |   T_STATEIN '=' boolean
                        { current_class().set_ctor_statein($3); }
        ;

boolean:        integer { $$ = $1; }
            |   T_YES   { $$ = 1; }
            |   T_NO    { $$ = 0; }
            |   T_TRUE  { $$ = 1; }
            |   T_FALSE { $$ = 0; }
        ;

string_array:   '[' strings ']' { $$ = $2; }
        ;

strings:        strings string  { $$ = $1; $$->insert($$->end(), *($2)); }
        |                       { $$ = new vector<string>; }
        |       error           { $$ = new vector<string>; }
        ;

string:         partialstring   { $$ = $1; }
        ;

partialstring: partialstring T_CONCAT substring { $$ = new string(*$1 + *$3);
                                    delete $1; delete $3; }
        | substring               { $$ = $1; }
        ;

substring:      T_STRING  { $$ = new string($1); }
        ;

integer:        T_INTEGER { $$ = $1; }
        ;

%%
