%{

/* $Log$
 * Revision 1.2  1994/02/23 02:05:49  cljanss
 * Added a yyerror declaration.
 *
 * Revision 1.1.1.1  1993/12/29  12:53:58  etseidl
 * SC source tree 0.1
 *
 * Revision 1.3  1992/04/01  17:28:35  seidl
 * remove token names
 *
 * Revision 1.2  1992/03/30  23:08:34  seidl
 * merge in sandia changes
 *
 * Revision 1.4  91/09/28  16:40:32  cljanss
 * new naming convention is used to keep names <= 14 characters.
 * 
 * Revision 1.3  91/07/22  14:43:17  cljanss
 * ';' is required after default modules for consistency.
 * 
 * Revision 1.2  1991/07/19  14:42:33  cljanss
 * Changed the way that generation modules are selected.
 *
 * Revision 1.1  1991/06/15  21:14:16  janssen
 * Initial revision
 * */
static char *rcsid = "$Id$";

#include <stdio.h>
#include <tmpl.h>
#include "types.h"
#include "global.h"
#include "sgen_read.gbl"

void yyerror(char*);
%}

%union {
  index_list_t *il;
  member_list_t *ml;
  name_list_t *nl;
  member_t *m;
  char *str;
  int i;
  }

%token T_UNION T_STRUCT T_SIGNED T_UNSIGNED T_EQEQ T_DEFAULT_MODULES
%token <str> T_ID
%type <str> exclude add
%type <nl> mod_opt mod_list
%type <il> indices
%type <str>  index
%type <ml> umembers members union
%type <m>  umember member
%type <i>  pointer qualifier

%start input
%%

input:			decs
			;

decs:			decs dec
			|	decs mod_def
			|
			;

mod_def:		T_DEFAULT_MODULES '(' mod_list ')' ';'
										{ set_default_modules($3); }
			;

dec:			T_ID mod_opt '{' members '}' ';'
										{ declaration($1, $2, $4); }
			;

mod_opt:		'(' mod_list ')'		{ $$ = $2; }
			|							{ $$ = NULL; }
			;

mod_list:		mod_list add		{ $$ = add_module($1,$2); }
			|	mod_list exclude	{ $$ = exclude_module($1,$2); }
			|							{ $$ = NULL; }
			;

add:				T_ID				{ $$ = $1; }
			;

exclude:		'!' T_ID				{ $$ = $2; }
			;

members:		members member
										{ $$ = add_member($1,$2); }
			|	members union
										{ $$ = merge_members($1,$2); }
			|
										{ $$ = NULL; }
			;

member:			qualifier T_ID pointer T_ID indices ';'
					{ $$ = create_member($1,$2,$3,$4,$5,NULL,NULL,NULL); }
			/* Old style unions: */
			|	T_UNION '(' T_ID T_EQEQ T_ID ')'
				        qualifier T_ID pointer T_ID indices T_ID ';'
					{ $$ = create_member($7,$8,$9,$10,$11,$3,$5,$12); }
			;

/* New style union specification: */
union:			T_UNION '(' T_ID ')' '{' umembers '}' T_ID ';'
							{ $$ = setup_union($3,$8,$6); }
			;

umembers:		umembers umember
										{ $$ = add_member($1,$2); }
			|
										{ $$ = NULL; }
			;

umember:		T_ID ':' qualifier T_ID pointer T_ID indices ';'
					{ $$ = create_member($3,$4,$5,$6,$7,NULL,$1,NULL); }
			;

qualifier:		T_STRUCT				{ $$ = Q_STRUCT; }
			|	T_SIGNED				{ $$ = Q_SIGNED; }
			|	T_UNSIGNED				{ $$ = Q_UNSIGNED; }
			|							{ $$ = Q_NONE; }
			;

pointer		:	pointer '*' 			{ $$ = $1 + 1; }
			|							{ $$ = 0; }
			;

indices:		indices index
										{ $$ = add_index($1,$2); }
			|
										{ $$ = NULL; }
			;

index:			'[' T_ID ']'
										{ $$ = $2; }
			;

%%

int
yywrap()
{return 1;}

void
yyerror(s)
char *s;
{error(s);}
