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
     T_MOLECULE = 258,
     T_MULTIPLICITY = 259,
     T_CHARGE = 260,
     T_METHOD = 261,
     T_BASIS = 262,
     T_AUXBASIS = 263,
     T_EQUALS = 264,
     T_OPTIMIZE = 265,
     T_GRADIENT = 266,
     T_BEG_OPT = 267,
     T_END_OPT = 268,
     T_CARTESIAN = 269,
     T_INTERNAL = 270,
     T_REDUNDANT = 271,
     T_RESTART = 272,
     T_CHECKPOINT = 273,
     T_COLON = 274,
     T_XC = 275,
     T_SYMMETRY = 276,
     T_MEMORY = 277,
     T_BOHR = 278,
     T_ANGSTROM = 279,
     T_GRID = 280,
     T_FREQUENCIES = 281,
     T_DOCC = 282,
     T_SOCC = 283,
     T_FROZEN_DOCC = 284,
     T_FROZEN_UOCC = 285,
     T_ALPHA = 286,
     T_BETA = 287,
     T_OO_INPUT_KEYWORD = 288,
     T_STRING = 289,
     T_BOOL = 290
   };
#endif
#define T_MOLECULE 258
#define T_MULTIPLICITY 259
#define T_CHARGE 260
#define T_METHOD 261
#define T_BASIS 262
#define T_AUXBASIS 263
#define T_EQUALS 264
#define T_OPTIMIZE 265
#define T_GRADIENT 266
#define T_BEG_OPT 267
#define T_END_OPT 268
#define T_CARTESIAN 269
#define T_INTERNAL 270
#define T_REDUNDANT 271
#define T_RESTART 272
#define T_CHECKPOINT 273
#define T_COLON 274
#define T_XC 275
#define T_SYMMETRY 276
#define T_MEMORY 277
#define T_BOHR 278
#define T_ANGSTROM 279
#define T_GRID 280
#define T_FREQUENCIES 281
#define T_DOCC 282
#define T_SOCC 283
#define T_FROZEN_DOCC 284
#define T_FROZEN_UOCC 285
#define T_ALPHA 286
#define T_BETA 287
#define T_OO_INPUT_KEYWORD 288
#define T_STRING 289
#define T_BOOL 290




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 28 "parse.yy"
typedef union YYSTYPE {
  char *str;
  int i;
  std::vector<int> *nniv;
  } YYSTYPE;
/* Line 1240 of yacc.c.  */
#line 112 "parse.tmp.hh"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE MPQCInylval;



