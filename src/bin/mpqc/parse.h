/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

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
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     T_NOT = 258,
     T_MOLECULE = 259,
     T_MULTIPLICITY = 260,
     T_CHARGE = 261,
     T_METHOD = 262,
     T_BASIS = 263,
     T_AUXBASIS = 264,
     T_DFBASIS = 265,
     T_EQUALS = 266,
     T_OPTIMIZE = 267,
     T_GRADIENT = 268,
     T_BEG_OPT = 269,
     T_END_OPT = 270,
     T_CARTESIAN = 271,
     T_INTERNAL = 272,
     T_REDUNDANT = 273,
     T_RESTART = 274,
     T_CHECKPOINT = 275,
     T_COLON = 276,
     T_SYMMETRY = 277,
     T_MEMORY = 278,
     T_TMPDIR = 279,
     T_TMPSTORE = 280,
     T_DEBUG = 281,
     T_ACCURACY = 282,
     T_BOHR = 283,
     T_ANGSTROM = 284,
     T_FREQUENCIES = 285,
     T_LINDEP = 286,
     T_MAXITER = 287,
     T_SCF = 288,
     T_UC = 289,
     T_DOCC = 290,
     T_SOCC = 291,
     T_FROZEN_DOCC = 292,
     T_FROZEN_UOCC = 293,
     T_ALPHA = 294,
     T_BETA = 295,
     T_XC = 296,
     T_GRID = 297,
     T_RI = 298,
     T_F12 = 299,
     T_APP = 300,
     T_ANSATZ = 301,
     T_OO_INPUT_KEYWORD = 302,
     T_STRING = 303,
     T_BOOL = 304
   };
#endif
/* Tokens.  */
#define T_NOT 258
#define T_MOLECULE 259
#define T_MULTIPLICITY 260
#define T_CHARGE 261
#define T_METHOD 262
#define T_BASIS 263
#define T_AUXBASIS 264
#define T_DFBASIS 265
#define T_EQUALS 266
#define T_OPTIMIZE 267
#define T_GRADIENT 268
#define T_BEG_OPT 269
#define T_END_OPT 270
#define T_CARTESIAN 271
#define T_INTERNAL 272
#define T_REDUNDANT 273
#define T_RESTART 274
#define T_CHECKPOINT 275
#define T_COLON 276
#define T_SYMMETRY 277
#define T_MEMORY 278
#define T_TMPDIR 279
#define T_TMPSTORE 280
#define T_DEBUG 281
#define T_ACCURACY 282
#define T_BOHR 283
#define T_ANGSTROM 284
#define T_FREQUENCIES 285
#define T_LINDEP 286
#define T_MAXITER 287
#define T_SCF 288
#define T_UC 289
#define T_DOCC 290
#define T_SOCC 291
#define T_FROZEN_DOCC 292
#define T_FROZEN_UOCC 293
#define T_ALPHA 294
#define T_BETA 295
#define T_XC 296
#define T_GRID 297
#define T_RI 298
#define T_F12 299
#define T_APP 300
#define T_ANSATZ 301
#define T_OO_INPUT_KEYWORD 302
#define T_STRING 303
#define T_BOOL 304




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 28 "parse.yy"
{
  char *str;
  int i;
  std::vector<int> *nniv;
  }
/* Line 1529 of yacc.c.  */
#line 153 "parse.tmp.hh"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE MPQCInylval;

