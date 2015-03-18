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
#include "mpqcin.h"
#define yyerror sc::MPQCIn::yerror
#define yyparse sc::MPQCIn::yparse
#define yylex sc::MPQCIn::ylex
#define yynerrs MPQCInyynerrs
#define yychar MPQCInyychar
%}

%union {
  char *str;
  int i;
  std::vector<int> *nniv;
  }

%token T_NOT
%token T_MOLECULE T_MULTIPLICITY T_CHARGE T_METHOD T_BASIS T_AUXBASIS T_DFBASIS T_EQUALS
%token T_OPTIMIZE T_GRADIENT T_BEG_OPT T_END_OPT T_CARTESIAN T_INTERNAL T_CONVERGENCE
%token T_REDUNDANT T_RESTART T_CHECKPOINT T_COLON T_SYMMETRY T_MEMORY T_TMPDIR T_TMPSTORE
%token T_DEBUG T_ACCURACY T_BOHR T_ANGSTROM T_FREQUENCIES T_PRECISE_FINDIF T_LINDEP T_MAXITER
%token T_SCF T_UC
%token T_DOCC T_SOCC T_FROZEN_DOCC T_FROZEN_UOCC T_ALPHA T_BETA
%token T_PCCSD
%token T_XC T_GRID
%token T_RI T_F12 T_APP T_ANSATZ
%token T_OO_INPUT_KEYWORD
%token <str> T_STRING
%token <i> T_BOOL
%type <str> string
%type <i> bool
%type <nniv> nonnegative_int_vector nonnegative_int_sequence

%start input
%%

input:          assignments
            ;

assignments:    assignments assignment
            |
            ;

assignment:     T_MOLECULE T_COLON              { begin_molecule(); }
                 molecule                       { end_molecule(); }
            |   T_SYMMETRY T_COLON string
                                                { set_symmetry($3); }
            |   T_MULTIPLICITY T_COLON string
                                                { set_multiplicity($3); }
            |   T_CHARGE T_COLON string
                                                { set_charge($3); }
            |   T_METHOD T_COLON string method_options_list
                                                { set_method($3); }
            |   T_BASIS T_COLON string basis_options_list
                                                { basis_.set_name($3); }
            |   T_AUXBASIS T_COLON string abasis_options_list
                                                { auxbasis_.set_name($3); }
            |   T_DFBASIS T_COLON string dbasis_options_list
                                                { dfbasis_.set_name($3); }
            |   T_OPTIMIZE T_COLON bool optimize_options_list
                                                { set_optimize($3); }
            |   T_GRADIENT T_COLON bool
                                                { set_gradient($3); }
            |   T_FREQUENCIES T_COLON bool freq_options_list
                                                { set_frequencies($3); }
            |   T_RESTART T_COLON bool
                                                { set_restart($3); }
            |   T_CHECKPOINT T_COLON bool
                                                { set_checkpoint($3); }
            |   T_PRECISE_FINDIF T_COLON bool
                                                { set_precise_findif($3); }
            |   T_MEMORY T_COLON string
                                                { set_memory($3); }
            |   T_TMPSTORE T_COLON string
                                                { set_tmpstore($3); }
            |   T_TMPDIR T_COLON string
                                                { set_tmpdir($3); }
            |   T_ACCURACY T_COLON string
                                                { set_accuracy($3); }
            |   T_LINDEP T_COLON string
                                                { set_lindep($3); }
            |   T_DOCC T_COLON nonnegative_int_vector
                                                { set_docc($3); }
            |   T_SOCC T_COLON nonnegative_int_vector
                                                { set_socc($3); }
            |   T_ALPHA T_COLON nonnegative_int_vector
                                                { set_alpha($3); }
            |   T_BETA T_COLON nonnegative_int_vector
                                                { set_beta($3); }
            |   T_FROZEN_DOCC T_COLON nonnegative_int_vector
                                                { set_frozen_docc($3); }
            |   T_FROZEN_UOCC T_COLON nonnegative_int_vector
                                                { set_frozen_uocc($3); }
            |   T_PCCSD T_COLON string string string
                                                { set_pccsd($3,$4,$5); }
            |   T_SCF T_COLON scf_options_list
            |   T_DEBUG T_COLON string
                                                { set_debug($3); }
            ;

nonnegative_int_vector:
                string                          { $$ = make_nnivec(0,$1); }
            |   T_BEG_OPT nonnegative_int_sequence T_END_OPT
                                                { $$ = $2; }
            ;

nonnegative_int_sequence:
                nonnegative_int_sequence string { $$ = make_nnivec($1,$2); }
            |                                   { $$ = make_nnivec(0,0); }
            ;

optimize_options_list:
                T_BEG_OPT optimize_options T_END_OPT
            |
            ;

optimize_options:
                optimize_options optimize_option
            |
            ;

optimize_option:
                T_CARTESIAN                     { set_opt_type(T_CARTESIAN); }
            |   T_INTERNAL                      { set_opt_type(T_INTERNAL); }
            |   T_REDUNDANT                     { set_redund_coor(1); }
            |   T_CONVERGENCE T_EQUALS string   { set_opt_convergence($3); }
            ;

freq_options_list:
                T_BEG_OPT freq_options T_END_OPT
            |
            ;

freq_options:
                freq_options freq_option
            |
            ;

freq_option:
                T_ACCURACY T_EQUALS string   { set_freq_accuracy($3); }
            ;

molecule:       molecule_options_list atoms
            ;

atoms:          atoms atom
            |
            ;

atom:           string string string string atom_options_list
                                                { add_atom($1,$2,$3,$4); }
            ;

atom_options_list:
                T_BEG_OPT atom_options T_END_OPT
            |
            ;

atom_options:
                atom_options atom_option
            |
            ;

atom_option:
                T_CHARGE T_EQUALS string        { set_atom_charge($3); }
            ;

molecule_options_list:
                T_BEG_OPT molecule_options T_END_OPT
            |
            ;

molecule_options:
                molecule_options molecule_option
            |
            ;

molecule_option:
                T_BOHR            { set_molecule_bohr(1); }
            |   T_ANGSTROM        { set_molecule_bohr(0); }
            ;

method_options_list:
                T_BEG_OPT method_options T_END_OPT
            |
            ;

method_options:
                method_options method_option
            |
            ;

method_option:
                T_XC T_EQUALS string        { set_dftmethod_xc($3); }
            |   T_GRID T_EQUALS string      { set_dftmethod_grid($3); }
            |   T_F12 T_EQUALS string       { set_r12method_f12($3); }
            |   T_APP T_EQUALS string       { set_r12method_app($3); }
            |   T_RI T_EQUALS string        { set_r12method_ri($3); }
            |   T_ANSATZ T_EQUALS string    { set_r12method_ansatz($3); }
            ;

basis_options_list:
                T_BEG_OPT basis_options T_END_OPT
            |
            ;

basis_options:
                basis_options basis_option
            |
            ;

basis_option:
                T_UC              { basis_.set_uc(true); }
            ;

abasis_options_list:
                T_BEG_OPT abasis_options T_END_OPT
            |
            ;

abasis_options:
                abasis_options abasis_option
            |
            ;

abasis_option:
                T_UC              { auxbasis_.set_uc(true); }
            ;

dbasis_options_list:
                T_BEG_OPT dbasis_options T_END_OPT
            |
            ;

dbasis_options:
                dbasis_options dbasis_option
            |
            ;

dbasis_option:
                T_UC              { dfbasis_.set_uc(true); }
            ;

scf_options_list:
                T_BEG_OPT scf_options T_END_OPT
            |
            ;

scf_options:
                scf_options scf_option
            |
            ;

scf_option:
                T_MAXITER T_EQUALS string   { set_scf_maxiter($3); }
            ;

string:         T_STRING                        { $$ = $1; }
            ;

bool:           T_BOOL                          { $$ = $1; }
            ;

%%
