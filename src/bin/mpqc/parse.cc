/* A Bison parser, made by GNU Bison 2.7.12-4996.  */

/* Bison implementation for Yacc-like parsers in C
   
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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.7.12-4996"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
/* Line 371 of yacc.c  */
#line 1 "parse.yy"

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
#include "mpqcin.h"
#define yyerror sc::MPQCIn::yerror
#define yyparse sc::MPQCIn::yparse
#define yylex sc::MPQCIn::ylex
#define yynerrs MPQCInyynerrs
#define yychar MPQCInyychar

/* Line 371 of yacc.c  */
#line 95 "parse.tmp.cc"

# ifndef YY_NULL
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULL nullptr
#  else
#   define YY_NULL 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "parse.tmp.hh".  */
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
/* Line 387 of yacc.c  */
#line 28 "parse.yy"

  char *str;
  int i;
  std::vector<int> *nniv;
  

/* Line 387 of yacc.c  */
#line 199 "parse.tmp.cc"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE MPQCInylval;

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus

#else

#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus

#else

#endif
#endif /* ! YYPARSE_PARAM */

#endif /* !YY_YY_PARSE_TMP_HH_INCLUDED  */

/* Copy the second part of user declarations.  */

/* Line 390 of yacc.c  */
#line 227 "parse.tmp.cc"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef __attribute__
/* This feature is available in gcc versions 2.5 and later.  */
# if (! defined __GNUC__ || __GNUC__ < 2 \
      || (__GNUC__ == 2 && __GNUC_MINOR__ < 5))
#  define __attribute__(Spec) /* empty */
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif


/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(N) (N)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (YYID (0))
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  3
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   157

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  55
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  39
/* YYNRULES -- Number of rules.  */
#define YYNRULES  103
/* YYNRULES -- Number of states.  */
#define YYNSTATES  193

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   309

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     5,     8,     9,    10,    15,    19,    23,
      27,    32,    37,    42,    47,    52,    56,    61,    65,    69,
      73,    77,    81,    85,    89,    93,    97,   101,   105,   109,
     113,   117,   123,   127,   131,   133,   137,   140,   141,   145,
     146,   149,   150,   152,   154,   156,   160,   164,   165,   168,
     169,   173,   176,   179,   180,   186,   190,   191,   194,   195,
     199,   203,   204,   207,   208,   210,   212,   216,   217,   220,
     221,   225,   229,   233,   237,   241,   245,   249,   250,   253,
     254,   256,   258,   260,   264,   265,   268,   269,   271,   273,
     275,   279,   280,   283,   284,   286,   288,   290,   294,   295,
     298,   299,   303,   305
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      56,     0,    -1,    57,    -1,    57,    58,    -1,    -1,    -1,
       4,    22,    59,    68,    -1,    23,    22,    92,    -1,     5,
      22,    92,    -1,     6,    22,    92,    -1,     7,    22,    92,
      77,    -1,     8,    22,    92,    80,    -1,     9,    22,    92,
      83,    -1,    10,    22,    92,    86,    -1,    12,    22,    93,
      62,    -1,    13,    22,    93,    -1,    31,    22,    93,    65,
      -1,    20,    22,    93,    -1,    21,    22,    93,    -1,    32,
      22,    93,    -1,    24,    22,    92,    -1,    26,    22,    92,
      -1,    25,    22,    92,    -1,    28,    22,    92,    -1,    33,
      22,    92,    -1,    39,    22,    60,    -1,    40,    22,    60,
      -1,    43,    22,    60,    -1,    44,    22,    60,    -1,    41,
      22,    60,    -1,    42,    22,    60,    -1,    45,    22,    92,
      92,    92,    -1,    35,    22,    89,    -1,    27,    22,    92,
      -1,    92,    -1,    14,    61,    15,    -1,    61,    92,    -1,
      -1,    14,    63,    15,    -1,    -1,    63,    64,    -1,    -1,
      16,    -1,    17,    -1,    19,    -1,    18,    11,    92,    -1,
      14,    66,    15,    -1,    -1,    66,    67,    -1,    -1,    28,
      11,    92,    -1,    74,    69,    -1,    69,    70,    -1,    -1,
      92,    92,    92,    92,    71,    -1,    14,    72,    15,    -1,
      -1,    72,    73,    -1,    -1,     6,    11,    92,    -1,    14,
      75,    15,    -1,    -1,    75,    76,    -1,    -1,    29,    -1,
      30,    -1,    14,    78,    15,    -1,    -1,    78,    79,    -1,
      -1,    46,    11,    92,    -1,    47,    11,    92,    -1,    49,
      11,    92,    -1,    50,    11,    92,    -1,    48,    11,    92,
      -1,    51,    11,    92,    -1,    14,    81,    15,    -1,    -1,
      81,    82,    -1,    -1,    36,    -1,    37,    -1,    38,    -1,
      14,    84,    15,    -1,    -1,    84,    85,    -1,    -1,    36,
      -1,    37,    -1,    38,    -1,    14,    87,    15,    -1,    -1,
      87,    88,    -1,    -1,    36,    -1,    37,    -1,    38,    -1,
      14,    90,    15,    -1,    -1,    90,    91,    -1,    -1,    34,
      11,    92,    -1,    53,    -1,    54,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,    54,    54,    57,    58,    61,    61,    63,    65,    67,
      69,    71,    73,    75,    77,    79,    81,    83,    85,    87,
      89,    91,    93,    95,    97,    99,   101,   103,   105,   107,
     109,   111,   113,   114,   119,   120,   125,   126,   130,   131,
     135,   136,   140,   141,   142,   143,   147,   148,   152,   153,
     157,   160,   163,   164,   167,   172,   173,   177,   178,   182,
     186,   187,   191,   192,   196,   197,   201,   202,   206,   207,
     211,   212,   213,   214,   215,   216,   220,   221,   225,   226,
     230,   231,   232,   236,   237,   241,   242,   246,   247,   248,
     252,   253,   257,   258,   262,   263,   264,   268,   269,   273,
     274,   278,   281,   284
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "T_NOT", "T_MOLECULE", "T_MULTIPLICITY",
  "T_CHARGE", "T_METHOD", "T_BASIS", "T_AUXBASIS", "T_DFBASIS", "T_EQUALS",
  "T_OPTIMIZE", "T_GRADIENT", "T_BEG_OPT", "T_END_OPT", "T_CARTESIAN",
  "T_INTERNAL", "T_CONVERGENCE", "T_REDUNDANT", "T_RESTART",
  "T_CHECKPOINT", "T_COLON", "T_SYMMETRY", "T_MEMORY", "T_TMPDIR",
  "T_TMPSTORE", "T_DEBUG", "T_ACCURACY", "T_BOHR", "T_ANGSTROM",
  "T_FREQUENCIES", "T_PRECISE_FINDIF", "T_LINDEP", "T_MAXITER", "T_SCF",
  "T_UC", "T_PUREAM", "T_SPLIT", "T_DOCC", "T_SOCC", "T_FROZEN_DOCC",
  "T_FROZEN_UOCC", "T_ALPHA", "T_BETA", "T_PCCSD", "T_XC", "T_GRID",
  "T_RI", "T_F12", "T_APP", "T_ANSATZ", "T_OO_INPUT_KEYWORD", "T_STRING",
  "T_BOOL", "$accept", "input", "assignments", "assignment", "$@1",
  "nonnegative_int_vector", "nonnegative_int_sequence",
  "optimize_options_list", "optimize_options", "optimize_option",
  "freq_options_list", "freq_options", "freq_option", "molecule", "atoms",
  "atom", "atom_options_list", "atom_options", "atom_option",
  "molecule_options_list", "molecule_options", "molecule_option",
  "method_options_list", "method_options", "method_option",
  "basis_options_list", "basis_options", "basis_option",
  "abasis_options_list", "abasis_options", "abasis_option",
  "dbasis_options_list", "dbasis_options", "dbasis_option",
  "scf_options_list", "scf_options", "scf_option", "string", "bool", YY_NULL
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    55,    56,    57,    57,    59,    58,    58,    58,    58,
      58,    58,    58,    58,    58,    58,    58,    58,    58,    58,
      58,    58,    58,    58,    58,    58,    58,    58,    58,    58,
      58,    58,    58,    58,    60,    60,    61,    61,    62,    62,
      63,    63,    64,    64,    64,    64,    65,    65,    66,    66,
      67,    68,    69,    69,    70,    71,    71,    72,    72,    73,
      74,    74,    75,    75,    76,    76,    77,    77,    78,    78,
      79,    79,    79,    79,    79,    79,    80,    80,    81,    81,
      82,    82,    82,    83,    83,    84,    84,    85,    85,    85,
      86,    86,    87,    87,    88,    88,    88,    89,    89,    90,
      90,    91,    92,    93
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     2,     0,     0,     4,     3,     3,     3,
       4,     4,     4,     4,     4,     3,     4,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     5,     3,     3,     1,     3,     2,     0,     3,     0,
       2,     0,     1,     1,     1,     3,     3,     0,     2,     0,
       3,     2,     2,     0,     5,     3,     0,     2,     0,     3,
       3,     0,     2,     0,     1,     1,     3,     0,     2,     0,
       3,     3,     3,     3,     3,     3,     3,     0,     2,     0,
       1,     1,     1,     3,     0,     2,     0,     1,     1,     1,
       3,     0,     2,     0,     1,     1,     1,     3,     0,     2,
       0,     3,     1,     1
};

/* YYDEFACT[STATE-NAME] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       4,     0,     2,     1,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     3,     5,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    98,     0,     0,     0,     0,     0,     0,
       0,    61,   102,     8,     9,    67,    77,    84,    91,   103,
      39,    15,    17,    18,     7,    20,    22,    21,    33,    23,
      47,    19,    24,   100,    32,    37,    25,    34,    26,    29,
      30,    27,    28,     0,    63,     6,    53,    69,    10,    79,
      11,    86,    12,    93,    13,    41,    14,    49,    16,     0,
       0,     0,     0,    51,     0,     0,     0,     0,     0,     0,
      97,     0,    99,    35,    36,    31,    60,    64,    65,    62,
      52,     0,    66,     0,     0,     0,     0,     0,     0,    68,
      76,    80,    81,    82,    78,    83,    87,    88,    89,    85,
      90,    94,    95,    96,    92,    38,    42,    43,     0,    44,
      40,    46,     0,    48,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   101,     0,    70,    71,    74,    72,
      73,    75,    45,    50,    56,    58,    54,     0,     0,    55,
      57,     0,    59
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     1,     2,    32,    61,    86,   110,   106,   118,   160,
     108,   119,   163,    95,   113,   130,   186,   187,   190,    96,
     112,   129,    98,   114,   139,   100,   115,   144,   102,   116,
     149,   104,   117,   154,    84,   109,   122,    87,    70
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -35
static const yytype_int16 yypact[] =
{
     -35,    17,    23,   -35,   -13,    -2,     3,    12,    16,    19,
      20,    30,    31,    39,    58,    65,    66,    67,    73,    74,
      76,    90,    91,    92,    93,    94,    95,    96,    97,    98,
      99,   100,   -35,   -35,    70,    70,    70,    70,    70,    70,
      71,    71,    71,    71,    70,    70,    70,    70,    70,    70,
      71,    71,    70,   110,    -8,    -8,    -8,    -8,    -8,    -8,
      70,   112,   -35,   -35,   -35,   113,   114,   115,   126,   -35,
     128,   -35,   -35,   -35,   -35,   -35,   -35,   -35,   -35,   -35,
     129,   -35,   -35,   -35,   -35,   -35,   -35,   -35,   -35,   -35,
     -35,   -35,   -35,    70,   -35,   -35,   -35,   -35,   -35,   -35,
     -35,   -35,   -35,   -35,   -35,   -35,   -35,   -35,   -35,     6,
       4,    70,    -7,    70,    24,    45,    54,    63,    87,     9,
     -35,   133,   -35,   -35,   -35,   -35,   -35,   -35,   -35,   -35,
     -35,    70,   -35,   134,   135,   136,   137,   138,   139,   -35,
     -35,   -35,   -35,   -35,   -35,   -35,   -35,   -35,   -35,   -35,
     -35,   -35,   -35,   -35,   -35,   -35,   -35,   -35,   140,   -35,
     -35,   -35,   141,   -35,    70,    70,    70,    70,    70,    70,
      70,    70,    70,    70,   -35,    70,   -35,   -35,   -35,   -35,
     -35,   -35,   -35,   -35,   142,   -35,   -35,     1,   143,   -35,
     -35,    70,   -35
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -35,   -35,   -35,   -35,   -35,    52,   -35,   -35,   -35,   -35,
     -35,   -35,   -35,   -35,   -35,   -35,   -35,   -35,   -35,   -35,
     -35,   -35,   -35,   -35,   -35,   -35,   -35,   -35,   -35,   -35,
     -35,   -35,   -35,   -35,   -35,   -35,   -35,   -34,    43
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
      63,    64,    65,    66,    67,    68,    85,   188,   126,    33,
      74,    75,    76,    77,    78,    79,   189,     3,    82,   123,
      34,   120,   127,   128,   161,    35,    93,     4,     5,     6,
       7,     8,     9,    10,    36,    11,    12,   162,    37,   132,
     121,    38,    39,    13,    14,    62,    15,    16,    17,    18,
      19,    20,    40,    41,    21,    22,    23,    62,    24,   111,
     140,    42,    25,    26,    27,    28,    29,    30,    31,   145,
     133,   134,   135,   136,   137,   138,   124,   125,   150,   131,
      43,   141,   142,   143,    71,    72,    73,    44,    45,    46,
     146,   147,   148,    80,    81,    47,    48,   165,    49,   151,
     152,   153,   155,   156,   157,   158,   159,    88,    89,    90,
      91,    92,    50,    51,    52,    53,    54,    55,    56,    57,
      58,    59,    60,    62,    83,    69,    94,    97,    99,   101,
     174,   175,   176,   177,   178,   179,   180,   181,   182,   183,
     103,   184,   105,   107,   164,   166,   167,   168,   169,   170,
     171,   172,   173,     0,   191,     0,   185,   192
};

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-35)))

#define yytable_value_is_error(Yytable_value) \
  YYID (0)

static const yytype_int16 yycheck[] =
{
      34,    35,    36,    37,    38,    39,    14,     6,    15,    22,
      44,    45,    46,    47,    48,    49,    15,     0,    52,    15,
      22,    15,    29,    30,    15,    22,    60,     4,     5,     6,
       7,     8,     9,    10,    22,    12,    13,    28,    22,    15,
      34,    22,    22,    20,    21,    53,    23,    24,    25,    26,
      27,    28,    22,    22,    31,    32,    33,    53,    35,    93,
      15,    22,    39,    40,    41,    42,    43,    44,    45,    15,
      46,    47,    48,    49,    50,    51,   110,   111,    15,   113,
      22,    36,    37,    38,    41,    42,    43,    22,    22,    22,
      36,    37,    38,    50,    51,    22,    22,   131,    22,    36,
      37,    38,    15,    16,    17,    18,    19,    55,    56,    57,
      58,    59,    22,    22,    22,    22,    22,    22,    22,    22,
      22,    22,    22,    53,    14,    54,    14,    14,    14,    14,
     164,   165,   166,   167,   168,   169,   170,   171,   172,   173,
      14,   175,    14,    14,    11,    11,    11,    11,    11,    11,
      11,    11,    11,    -1,    11,    -1,    14,   191
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    56,    57,     0,     4,     5,     6,     7,     8,     9,
      10,    12,    13,    20,    21,    23,    24,    25,    26,    27,
      28,    31,    32,    33,    35,    39,    40,    41,    42,    43,
      44,    45,    58,    22,    22,    22,    22,    22,    22,    22,
      22,    22,    22,    22,    22,    22,    22,    22,    22,    22,
      22,    22,    22,    22,    22,    22,    22,    22,    22,    22,
      22,    59,    53,    92,    92,    92,    92,    92,    92,    54,
      93,    93,    93,    93,    92,    92,    92,    92,    92,    92,
      93,    93,    92,    14,    89,    14,    60,    92,    60,    60,
      60,    60,    60,    92,    14,    68,    74,    14,    77,    14,
      80,    14,    83,    14,    86,    14,    62,    14,    65,    90,
      61,    92,    75,    69,    78,    81,    84,    87,    63,    66,
      15,    34,    91,    15,    92,    92,    15,    29,    30,    76,
      70,    92,    15,    46,    47,    48,    49,    50,    51,    79,
      15,    36,    37,    38,    82,    15,    36,    37,    38,    85,
      15,    36,    37,    38,    88,    15,    16,    17,    18,    19,
      64,    15,    28,    67,    11,    92,    11,    11,    11,    11,
      11,    11,    11,    11,    92,    92,    92,    92,    92,    92,
      92,    92,    92,    92,    92,    14,    71,    72,     6,    15,
      73,    11,    92
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  However,
   YYFAIL appears to be in use.  Nevertheless, it is formally deprecated
   in Bison 2.4.2's NEWS entry, where a plan to phase it out is
   discussed.  */

#define YYFAIL		goto yyerrlab
#if defined YYFAIL
  /* This is here to suppress warnings from the GCC cpp's
     -Wunused-macros.  Normally we don't worry about that warning, but
     some users do, and we want to make it easy for users to remove
     YYFAIL uses, which will produce warnings from Bison 2.5.  */
#endif

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      MPQCInylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))

/* Error token number */
#define YYTERROR	1
#define YYERRCODE	256


/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */
#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULL, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULL;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - Assume YYFAIL is not used.  It's too flawed to consider.  See
       <http://lists.gnu.org/archive/html/bison-patches/2009-12/msg00024.html>
       for details.  YYERROR is fine as it does not invoke this
       function.
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULL, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YYUSE (yytype);
}




/* The lookahead symbol.  */
int yychar;


#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

/* The semantic value of the lookahead symbol.  */
YYSTYPE MPQCInylval YY_INITIAL_VALUE(yyval_default);

/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &MPQCInylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &MPQCInylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = MPQCInylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 5:
/* Line 1787 of yacc.c  */
#line 61 "parse.yy"
    { begin_molecule(); }
    break;

  case 6:
/* Line 1787 of yacc.c  */
#line 62 "parse.yy"
    { end_molecule(); }
    break;

  case 7:
/* Line 1787 of yacc.c  */
#line 64 "parse.yy"
    { set_symmetry((yyvsp[(3) - (3)].str)); }
    break;

  case 8:
/* Line 1787 of yacc.c  */
#line 66 "parse.yy"
    { set_multiplicity((yyvsp[(3) - (3)].str)); }
    break;

  case 9:
/* Line 1787 of yacc.c  */
#line 68 "parse.yy"
    { set_charge((yyvsp[(3) - (3)].str)); }
    break;

  case 10:
/* Line 1787 of yacc.c  */
#line 70 "parse.yy"
    { set_method((yyvsp[(3) - (4)].str)); }
    break;

  case 11:
/* Line 1787 of yacc.c  */
#line 72 "parse.yy"
    { basis_.set_name((yyvsp[(3) - (4)].str)); }
    break;

  case 12:
/* Line 1787 of yacc.c  */
#line 74 "parse.yy"
    { auxbasis_.set_name((yyvsp[(3) - (4)].str)); }
    break;

  case 13:
/* Line 1787 of yacc.c  */
#line 76 "parse.yy"
    { dfbasis_.set_name((yyvsp[(3) - (4)].str)); }
    break;

  case 14:
/* Line 1787 of yacc.c  */
#line 78 "parse.yy"
    { set_optimize((yyvsp[(3) - (4)].i)); }
    break;

  case 15:
/* Line 1787 of yacc.c  */
#line 80 "parse.yy"
    { set_gradient((yyvsp[(3) - (3)].i)); }
    break;

  case 16:
/* Line 1787 of yacc.c  */
#line 82 "parse.yy"
    { set_frequencies((yyvsp[(3) - (4)].i)); }
    break;

  case 17:
/* Line 1787 of yacc.c  */
#line 84 "parse.yy"
    { set_restart((yyvsp[(3) - (3)].i)); }
    break;

  case 18:
/* Line 1787 of yacc.c  */
#line 86 "parse.yy"
    { set_checkpoint((yyvsp[(3) - (3)].i)); }
    break;

  case 19:
/* Line 1787 of yacc.c  */
#line 88 "parse.yy"
    { set_precise_findif((yyvsp[(3) - (3)].i)); }
    break;

  case 20:
/* Line 1787 of yacc.c  */
#line 90 "parse.yy"
    { set_memory((yyvsp[(3) - (3)].str)); }
    break;

  case 21:
/* Line 1787 of yacc.c  */
#line 92 "parse.yy"
    { set_tmpstore((yyvsp[(3) - (3)].str)); }
    break;

  case 22:
/* Line 1787 of yacc.c  */
#line 94 "parse.yy"
    { set_tmpdir((yyvsp[(3) - (3)].str)); }
    break;

  case 23:
/* Line 1787 of yacc.c  */
#line 96 "parse.yy"
    { set_accuracy((yyvsp[(3) - (3)].str)); }
    break;

  case 24:
/* Line 1787 of yacc.c  */
#line 98 "parse.yy"
    { set_lindep((yyvsp[(3) - (3)].str)); }
    break;

  case 25:
/* Line 1787 of yacc.c  */
#line 100 "parse.yy"
    { set_docc((yyvsp[(3) - (3)].nniv)); }
    break;

  case 26:
/* Line 1787 of yacc.c  */
#line 102 "parse.yy"
    { set_socc((yyvsp[(3) - (3)].nniv)); }
    break;

  case 27:
/* Line 1787 of yacc.c  */
#line 104 "parse.yy"
    { set_alpha((yyvsp[(3) - (3)].nniv)); }
    break;

  case 28:
/* Line 1787 of yacc.c  */
#line 106 "parse.yy"
    { set_beta((yyvsp[(3) - (3)].nniv)); }
    break;

  case 29:
/* Line 1787 of yacc.c  */
#line 108 "parse.yy"
    { set_frozen_docc((yyvsp[(3) - (3)].nniv)); }
    break;

  case 30:
/* Line 1787 of yacc.c  */
#line 110 "parse.yy"
    { set_frozen_uocc((yyvsp[(3) - (3)].nniv)); }
    break;

  case 31:
/* Line 1787 of yacc.c  */
#line 112 "parse.yy"
    { set_pccsd((yyvsp[(3) - (5)].str),(yyvsp[(4) - (5)].str),(yyvsp[(5) - (5)].str)); }
    break;

  case 33:
/* Line 1787 of yacc.c  */
#line 115 "parse.yy"
    { set_debug((yyvsp[(3) - (3)].str)); }
    break;

  case 34:
/* Line 1787 of yacc.c  */
#line 119 "parse.yy"
    { (yyval.nniv) = make_nnivec(0,(yyvsp[(1) - (1)].str)); }
    break;

  case 35:
/* Line 1787 of yacc.c  */
#line 121 "parse.yy"
    { (yyval.nniv) = (yyvsp[(2) - (3)].nniv); }
    break;

  case 36:
/* Line 1787 of yacc.c  */
#line 125 "parse.yy"
    { (yyval.nniv) = make_nnivec((yyvsp[(1) - (2)].nniv),(yyvsp[(2) - (2)].str)); }
    break;

  case 37:
/* Line 1787 of yacc.c  */
#line 126 "parse.yy"
    { (yyval.nniv) = make_nnivec(0,0); }
    break;

  case 42:
/* Line 1787 of yacc.c  */
#line 140 "parse.yy"
    { set_opt_type(T_CARTESIAN); }
    break;

  case 43:
/* Line 1787 of yacc.c  */
#line 141 "parse.yy"
    { set_opt_type(T_INTERNAL); }
    break;

  case 44:
/* Line 1787 of yacc.c  */
#line 142 "parse.yy"
    { set_redund_coor(1); }
    break;

  case 45:
/* Line 1787 of yacc.c  */
#line 143 "parse.yy"
    { set_opt_convergence((yyvsp[(3) - (3)].str)); }
    break;

  case 50:
/* Line 1787 of yacc.c  */
#line 157 "parse.yy"
    { set_freq_accuracy((yyvsp[(3) - (3)].str)); }
    break;

  case 54:
/* Line 1787 of yacc.c  */
#line 168 "parse.yy"
    { add_atom((yyvsp[(1) - (5)].str),(yyvsp[(2) - (5)].str),(yyvsp[(3) - (5)].str),(yyvsp[(4) - (5)].str)); }
    break;

  case 59:
/* Line 1787 of yacc.c  */
#line 182 "parse.yy"
    { set_atom_charge((yyvsp[(3) - (3)].str)); }
    break;

  case 64:
/* Line 1787 of yacc.c  */
#line 196 "parse.yy"
    { set_molecule_bohr(1); }
    break;

  case 65:
/* Line 1787 of yacc.c  */
#line 197 "parse.yy"
    { set_molecule_bohr(0); }
    break;

  case 70:
/* Line 1787 of yacc.c  */
#line 211 "parse.yy"
    { set_dftmethod_xc((yyvsp[(3) - (3)].str)); }
    break;

  case 71:
/* Line 1787 of yacc.c  */
#line 212 "parse.yy"
    { set_dftmethod_grid((yyvsp[(3) - (3)].str)); }
    break;

  case 72:
/* Line 1787 of yacc.c  */
#line 213 "parse.yy"
    { set_r12method_f12((yyvsp[(3) - (3)].str)); }
    break;

  case 73:
/* Line 1787 of yacc.c  */
#line 214 "parse.yy"
    { set_r12method_app((yyvsp[(3) - (3)].str)); }
    break;

  case 74:
/* Line 1787 of yacc.c  */
#line 215 "parse.yy"
    { set_r12method_ri((yyvsp[(3) - (3)].str)); }
    break;

  case 75:
/* Line 1787 of yacc.c  */
#line 216 "parse.yy"
    { set_r12method_ansatz((yyvsp[(3) - (3)].str)); }
    break;

  case 80:
/* Line 1787 of yacc.c  */
#line 230 "parse.yy"
    { basis_.set_uc(true); }
    break;

  case 81:
/* Line 1787 of yacc.c  */
#line 231 "parse.yy"
    { basis_.set_puream(true); }
    break;

  case 82:
/* Line 1787 of yacc.c  */
#line 232 "parse.yy"
    { basis_.set_split(true); }
    break;

  case 87:
/* Line 1787 of yacc.c  */
#line 246 "parse.yy"
    { auxbasis_.set_uc(true); }
    break;

  case 88:
/* Line 1787 of yacc.c  */
#line 247 "parse.yy"
    { auxbasis_.set_puream(true); }
    break;

  case 89:
/* Line 1787 of yacc.c  */
#line 248 "parse.yy"
    { auxbasis_.set_split(true); }
    break;

  case 94:
/* Line 1787 of yacc.c  */
#line 262 "parse.yy"
    { dfbasis_.set_uc(true); }
    break;

  case 95:
/* Line 1787 of yacc.c  */
#line 263 "parse.yy"
    { dfbasis_.set_puream(true); }
    break;

  case 96:
/* Line 1787 of yacc.c  */
#line 264 "parse.yy"
    { dfbasis_.set_split(true); }
    break;

  case 101:
/* Line 1787 of yacc.c  */
#line 278 "parse.yy"
    { set_scf_maxiter((yyvsp[(3) - (3)].str)); }
    break;

  case 102:
/* Line 1787 of yacc.c  */
#line 281 "parse.yy"
    { (yyval.str) = (yyvsp[(1) - (1)].str); }
    break;

  case 103:
/* Line 1787 of yacc.c  */
#line 284 "parse.yy"
    { (yyval.i) = (yyvsp[(1) - (1)].i); }
    break;


/* Line 1787 of yacc.c  */
#line 1948 "parse.tmp.cc"
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &MPQCInylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = MPQCInylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &MPQCInylval);
    }
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


/* Line 2050 of yacc.c  */
#line 287 "parse.yy"

