
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tmpl.h>
#include <ctype.h>

#include <chemistry/qc/intv3/atominfo.gbl>
#include <chemistry/qc/intv3/atominfo.lcl>

#define N_ATOMS 110

struct {
  char *atom;
  char *symbol;
  int an;
  } atominfo[N_ATOMS] = {
   {"hydrogen",    "h",   1},
   {"helium",      "he",  2},
   {"lithium",     "li",  3},
   {"beryllium",   "be",  4},
   {"boron",       "b",   5},
   {"carbon",      "c",   6},
   {"nitrogen",    "n",   7},
   {"oxygen",      "o",   8},
   {"fluorine",    "f",   9},
   {"neon",        "ne", 10},
   {"sodium",      "na", 11},
   {"magnesium",   "mg", 12},
   {"aluminum",    "al", 13},
   {"silicon",     "si", 14},
   {"phosphorus",  "p",  15},
   {"sulfur",      "s",  16},
   {"chlorine",    "cl", 17},
   {"argon",       "ar", 18},
   {"xenon",       "xe", 54},
   {NULL,          NULL,  0}
  };

/* Convert an atomic number to a symbol.  The returned character pointer is
 * malloced. */
char *
IntV3::an_to_sym(int an)
{
  int i;
  char *result;

  for (i=0; atominfo[i].an != 0; i++) {
    if (atominfo[i].an == an) {
      result = (char *)malloc(strlen(atominfo[i].symbol)+1);
      strcpy(result,atominfo[i].symbol);
      result[0] += 'A' - 'a';
      return result;
      }
    }
  return NULL;
  }

/* Converts a symbol to an atom name.  If the symbol name is unknown
 * then the symbol name is returned. */
char *
IntV3::sym_to_atom(char *sym)
{
  int i;
  char tmpsym[10];

  if (!sym) return NULL;

  /* Convert the passed name to lowercase. */
  strcpy(tmpsym,sym);
  for (i=0; i<strlen(sym); i++) {
    if (isupper(tmpsym[i])) tmpsym[i] += 'a' - 'A';
    }

  for (i=0; atominfo[i].atom != NULL; i++) { 
    if (!strcmp(tmpsym,atominfo[i].atom) || !strcmp(tmpsym,atominfo[i].symbol)) { 
      return atominfo[i].atom;
      }
    }
  return NULL;
  }

int
IntV3::atom_to_an(char *atom)
{
  int i;
  char tmpatom[50];

  if (!atom) return 0;

  /* Convert the passed name to lowercase. */
  strcpy(tmpatom,atom);
  for (i=0; i<strlen(atom); i++) {
    if (isupper(tmpatom[i])) tmpatom[i] += 'a' - 'A';
    }

  for (i=0; atominfo[i].atom != NULL; i++) { 
    if (!strcmp(tmpatom,atominfo[i].atom) || !strcmp(tmpatom,atominfo[i].symbol)) { 
      return atominfo[i].an;
      }
    }

  return 0;
  }

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
