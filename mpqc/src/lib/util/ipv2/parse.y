%{
#include <stdio.h>
#ifdef DEC
#include <math.h>
#else
#include <stdlib.h>
#endif
#include <tmpl.h>
#ifdef BISON
#define YYDEBUG 0
#if YYDEBUG != 0
int yydebug =1;
#endif /* YYDEBUG != 0 */
#endif /* BISON */
#if defined(SABER)
#define xmalloc malloc
#endif
#if defined(SGI) || defined(sparc)
#include <alloca.h>
#endif
#include "ip_types.h"
#include "ip_read.gbl"
#include "ip_error.gbl"
%}

%union {
  char *str;
  ip_string_list_t *sl;
  double dbl;
  }

%token T_KEYWORD_LEFT T_ARRAY_LEFT T_TABLE_LEFT
%token T_KEYWORD_RIGHT T_ARRAY_RIGHT T_TABLE_RIGHT
%token <str> T_STRING
%type <str> scalar var_key var_sub table_key polymorph
%type <sl> stringlist parentlist commalist
%type <str> expression
%type <dbl> subexpression expvalue

%start input
%%

input:			group_defs
			;

group_defs:		group_defs group_def
			|
			;

group_def:		keyword ':' T_KEYWORD_LEFT group_defs T_KEYWORD_RIGHT
												{ ip_pop_keyword(); }
			|	keyword ':' karray_defs
												{ ip_pop_keyword(); }
			|	keyword ':' group_def
												{ ip_pop_keyword(); }
			|	keyword '=' value
												{ ip_pop_keyword(); }
			|	T_TABLE_LEFT stringlist T_TABLE_RIGHT
												{ ip_begin_table($2); }
				 '='
				T_TABLE_LEFT tablevalues T_TABLE_RIGHT
												{ ip_done_table(); }
			;
/* old table construction stuff
			|	T_TABLE_LEFT stringlist T_TABLE_RIGHT
				 '='
				T_TABLE_LEFT stringlist T_TABLE_RIGHT
												{ ip_construct_table($2,$6); }
			;
 */

/* One or more strings in a stringlist.
 * This is needed to distinguish array construction from table contruction
 * (with no columns-which is silly anyway).  If we allow zero or more,
 * then YACC complains about the reduce/reduce conflict. */
stringlist:		stringlist table_key
							{ $$ = ip_add_string_list($1,$2); }
			|	table_key
							{ $$ = ip_string_to_string_list($1); }
			;


table_key:		table_key ':' T_STRING	{ $$ = ip_append_keystrings($1,$3); }
			|	T_STRING				{ $$ = $1; }
			;

tablevalues:	tablevalues value
			|
			;

karray_defs:	array_left karray_elems array_right
			;

karray_elems:	karray_elems { ip_incr_karray(); } karray_elem
			|									{ ip_init_karray(); }
			;

karray_elem:	karray_defs
			|	T_KEYWORD_LEFT group_defs T_KEYWORD_RIGHT
			|	group_def
			;

array_left:		T_ARRAY_LEFT					{ ip_start_karray(); }
			;

array_right:	T_ARRAY_RIGHT					{ ip_pop_karray(); }
			;

/*
keyword:		T_STRING
												{ ip_push_keyword($1); }
			|	T_STRING '<' T_STRING '>'		{ ip_push_keyclass($1,$3); }
			;
*/

keyword:		T_STRING 						{ ip_push_keyword($1); }
			|	T_STRING polymorph				{ ip_push_keyclass($1,$2,0); }
			|	T_STRING parentlist				{ ip_push_keyclass($1,0,$2); }
			|	T_STRING parentlist polymorph	{ ip_push_keyclass($1,$3,$2); }
			|	T_STRING polymorph parentlist	{ ip_push_keyclass($1,$2,$3); }
			;

polymorph:		'<' T_STRING '>'				{ $$ = $2; }
			;

parentlist:		'<' '<' commalist '>' '>'		{ $$ = $3; }
			;

commalist:		commalist ',' T_STRING	{ $$ = ip_add_string_list($1,$3); }
			|	T_STRING				{ $$ = ip_string_to_string_list($1); }
			;

value:			array
			|	scalar
									{ ip_assign_value($1); }
			|	var_sub				{ ip_assign_variable($1); }
			|	expression			{ ip_assign_value($1); }
			;

expression:		T_KEYWORD_LEFT subexpression T_KEYWORD_RIGHT
									{ $$ = ip_double_to_string($2); }
			;

subexpression:	expvalue '*' expvalue	{ $$ = $1 * $3; }
			|	expvalue '-' expvalue	{ $$ = $1 - $3; }
			|	expvalue '+' expvalue	{ $$ = $1 + $3; }
			|	expvalue '/' expvalue	{ $$ = $1 / $3; }
			;


expvalue:		var_sub				{ $$ = ip_get_variable_double($1); }
			|	scalar				{ $$ = atof($1); free($1); }
			; 

var_sub:		'$' var_key			{ $$ = $2; }
			;


var_key:		var_key ':' T_STRING	{ $$ = ip_append_keystrings($1,$3); }
			|	':' T_STRING			{ $$ = ip_append_keystrings(NULL,$2); }
			|	T_STRING				{ $$ = $1; }
			;

array:			array_left values array_right
			;

values:			values { ip_incr_karray(); } value
			|
												{ ip_init_karray(); }
			;

scalar:			T_STRING
												{ $$ = $1; }
			;

%%

int
yywrap()
{return 1;}

int
yyerror(s)
char *s;
{ip_error(s);}
