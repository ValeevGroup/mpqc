/* A Bison parser, made by GNU Bison 2.7.12-4996.  */

/* Bison interface for Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2013 Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

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

#ifndef YY_YY_PARSE_TMP_HH_INCLUDED
# define YY_YY_PARSE_TMP_HH_INCLUDED
/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

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
     T_CONVERGENCE = 273,
     T_REDUNDANT = 274,
     T_RESTART = 275,
     T_CHECKPOINT = 276,
     T_COLON = 277,
     T_SYMMETRY = 278,
     T_MEMORY = 279,
     T_TMPDIR = 280,
     T_TMPSTORE = 281,
     T_DEBUG = 282,
     T_ACCURACY = 283,
     T_BOHR = 284,
     T_ANGSTROM = 285,
     T_FREQUENCIES = 286,
     T_PRECISE_FINDIF = 287,
     T_LINDEP = 288,
     T_MAXITER = 289,
     T_SCF = 290,
     T_UC = 291,
     T_PUREAM = 292,
     T_SPLIT = 293,
     T_DOCC = 294,
     T_SOCC = 295,
     T_FROZEN_DOCC = 296,
     T_FROZEN_UOCC = 297,
     T_ALPHA = 298,
     T_BETA = 299,
     T_PCCSD = 300,
     T_XC = 301,
     T_GRID = 302,
     T_RI = 303,
     T_F12 = 304,
     T_APP = 305,
     T_ANSATZ = 306,
     T_OO_INPUT_KEYWORD = 307,
     T_STRING = 308,
     T_BOOL = 309
   };
#endif


#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{
/* Line 2053 of yacc.c  */
#line 28 "parse.yy"

  char *str;
  int i;
  std::vector<int> *nniv;
  

/* Line 2053 of yacc.c  */
#line 118 "parse.tmp.hh"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE MPQCInylval;

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */

#endif /* !YY_YY_PARSE_TMP_HH_INCLUDED  */
