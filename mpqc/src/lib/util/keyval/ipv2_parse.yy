/*
 * ipv2_parse.yy
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
%{
#ifdef DEC
#include <math.h>
#else
#include <stdlib.h>
#endif
#include <string.h>
#ifdef BISON
#define YYDEBUG 0
#if YYDEBUG != 0
int yydebug =1;
#endif /* YYDEBUG != 0 */
#endif /* BISON */
#if defined(SABER)
#define xmalloc malloc
#endif
#if defined(SGI)
#include <alloca.h>
#endif
#include <util/keyval/ipv2.h>
#define yyerror sc::IPV2::yerror
#define yyparse sc::IPV2::yparse
#define yylex sc::IPV2::ylex
#define yywrap sc::IPV2::ywrap
%}

%union {
  char *str;
  sc::ip_string_list_t *sl;
  double dbl;
  }

%token T_KEYWORD_LEFT T_ARRAY_LEFT T_TABLE_LEFT
%token T_KEYWORD_RIGHT T_ARRAY_RIGHT T_TABLE_RIGHT
%token T_CONCAT
%token <str> T_STRING T_QUOTED_STRING
%type <str> string quoted_string var_key var_sub table_key polymorph
%type <sl> stringlist parentlist commalist
%type <str> expression
%type <dbl> subexpression expvalue

%start input
%%

input:          group_defs
            ;

group_defs:     group_defs group_def
            |
            ;

group_def:      keyword ':' T_KEYWORD_LEFT group_defs T_KEYWORD_RIGHT
                                                { ip_pop_keyword(); }
            |   keyword ':' karray_defs
                                                { ip_pop_keyword(); }
            |   keyword ':' group_def
                                                { ip_pop_keyword(); }
            |   keyword '=' value
                                                { ip_pop_keyword(); }
            |   T_TABLE_LEFT stringlist T_TABLE_RIGHT
                                                { ip_begin_table($2); }
                 '='
                T_TABLE_LEFT tablevalues T_TABLE_RIGHT
                                                { ip_done_table(); }
            |   implicit_keyword ':' T_KEYWORD_LEFT group_defs T_KEYWORD_RIGHT
            |   implicit_keyword ':' karray_defs
            |   implicit_keyword ':' group_def
            |   implicit_keyword '=' value
            ;
/* old table construction stuff
            |   T_TABLE_LEFT stringlist T_TABLE_RIGHT
                 '='
                T_TABLE_LEFT stringlist T_TABLE_RIGHT
                                                { ip_construct_table($2,$6); }
            ;
 */

/* One or more strings in a stringlist.
 * This is needed to distinguish array construction from table contruction
 * (with no columns-which is silly anyway).  If we allow zero or more,
 * then YACC complains about the reduce/reduce conflict. */
stringlist:     stringlist table_key
                            { $$ = ip_add_string_list($1,$2); }
            |   table_key
                            { $$ = ip_string_to_string_list($1); }
            ;


table_key:      table_key ':' string    { $$ = ip_append_keystrings($1,$3); }
            |   string                  { $$ = $1; }
            ;

tablevalues:    tablevalues value
            |
            ;

karray_defs:    array_left karray_elems array_right
            ;

karray_elems:   karray_elems { ip_incr_karray(); } karray_elem
            |                                   { ip_init_karray(); }
            ;

karray_elem:    karray_defs
            |   T_KEYWORD_LEFT group_defs T_KEYWORD_RIGHT
            |   group_def
            ;

array_left:     T_ARRAY_LEFT                    { ip_start_karray(); }
            ;

array_right:    T_ARRAY_RIGHT                   { ip_pop_karray(); }
            ;

/*
keyword:        T_STRING
                                                { ip_push_keyword($1); }
            |   T_STRING '<' T_STRING '>'       { ip_push_keyclass($1,$3); }
            ;
*/

keyword:        string                        { ip_push_keyword($1); }
            |   string polymorph              { ip_push_keyclass($1,$2,0); }
            |   string parentlist             { ip_push_keyclass($1,0,$2); }
            |   string parentlist polymorph   { ip_push_keyclass($1,$3,$2); }
            |   string polymorph parentlist   { ip_push_keyclass($1,$2,$3); }
            ;
implicit_keyword:
                polymorph                       { ip_push_keyclass(0,$1,0); }
            ;

polymorph:      '<' string '>'                { $$ = $2; }
            ;

parentlist:     '<' '<' commalist '>' '>'       { $$ = $3; }
            ;

commalist:      commalist ',' string  { $$ = ip_add_string_list($1,$3); }
            |   string                { $$ = ip_string_to_string_list($1); }
            ;

value:          array
            |   string
                                    { ip_assign_value($1); }
            |   var_sub             { ip_assign_variable($1); }
            |   expression          { ip_assign_value($1); }
            ;

expression:     T_KEYWORD_LEFT subexpression T_KEYWORD_RIGHT
                                    { $$ = ip_double_to_string($2); }
            ;

subexpression:  expvalue '*' expvalue   { $$ = $1 * $3; }
            |   expvalue '-' expvalue   { $$ = $1 - $3; }
            |   expvalue '+' expvalue   { $$ = $1 + $3; }
            |   expvalue '/' expvalue   { $$ = $1 / $3; }
            ;


expvalue:       var_sub             { $$ = ip_get_variable_double($1); }
            |   string              { $$ = atof($1); free($1); }
            ; 

var_sub:        '$' var_key         { $$ = $2; }
            ;


var_key:        var_key ':' string    { $$ = ip_append_keystrings($1,$3); }
            |   ':' string            { $$ = ip_append_keystrings(NULL,$2); }
            |   string                { $$ = $1; }
            ;

array:          array_left values array_right
            ;

values:         values { ip_incr_karray(); } value
            |
                                                { ip_init_karray(); }
            ;

quoted_string: quoted_string T_CONCAT T_QUOTED_STRING
                                { $$ = (char*) malloc(strlen($1)+strlen($3)+1);
                                  strcpy($$, $1);
                                  strcat($$, $3);
                                  free($1);
                                  free($3);
                                }
            |  T_QUOTED_STRING
                                                { $$ = $1; }
            ;

string:        T_STRING                         { $$ = $1; }
            |  quoted_string                    { $$ = $1; }
            ;

%%
