/* A Bison parser, made by GNU Bison 1.875.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     T_KEYWORD_LEFT = 258,
     T_ARRAY_LEFT = 259,
     T_TABLE_LEFT = 260,
     T_KEYWORD_RIGHT = 261,
     T_ARRAY_RIGHT = 262,
     T_TABLE_RIGHT = 263,
     T_CONCAT = 264,
     T_STRING = 265,
     T_QUOTED_STRING = 266
   };
#endif
#define T_KEYWORD_LEFT 258
#define T_ARRAY_LEFT 259
#define T_TABLE_LEFT 260
#define T_KEYWORD_RIGHT 261
#define T_ARRAY_RIGHT 262
#define T_TABLE_RIGHT 263
#define T_CONCAT 264
#define T_STRING 265
#define T_QUOTED_STRING 266




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 53 "ipv2_parse.yy"
typedef union YYSTYPE {
  char *str;
  sc::ip_string_list_t *sl;
  double dbl;
  } YYSTYPE;
/* Line 1240 of yacc.c.  */
#line 64 "ipv2_parse.tmp.hh"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;



