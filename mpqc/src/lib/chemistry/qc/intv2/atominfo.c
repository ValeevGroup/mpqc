
/* $Log$
 * Revision 1.3  1994/08/26 22:45:07  etseidl
 * fix a bunch of warnings, get rid of rcs id's, get rid of bread/bwrite and
 * fread/fwrite modules
 *
 * Revision 1.2  1994/01/14  10:50:51  seidl
 * add xenon to atominfo struct
 *
 * Revision 1.1.1.1  1993/12/29  12:53:02  etseidl
 * SC source tree 0.1
 *
 * Revision 1.3  1992/06/17  22:04:21  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/06/02  10:57:35  seidl
 * change "phosphorous" to "phosphorus"
 *
 * Revision 1.1.1.1  1992/03/17  16:32:01  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:32:00  seidl
 * Initial revision
 *
 * Revision 1.1  1992/02/03  15:49:49  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/14  00:19:36  cljanss
 * Initial revision
 *
 * Revision 1.5  91/10/31  14:43:30  cljanss
 * added ctype.h
 * 
 * Revision 1.4  91/09/28  19:26:41  cljanss
 * Switch to new file naming scheme
 * 
 * Revision 1.3  91/09/10  19:32:02  cljanss
 * sym_to_atom now returns NULL when it can't find the sym
 * 
 * Revision 1.2  1991/08/13  17:47:43  cljanss
 * routines are now case independent
 *
 * Revision 1.1  1991/06/16  16:40:07  janssen
 * Initial revision
 * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tmpl.h>
#include <ctype.h>

#include "atominfo.gbl"
#include "atominfo.lcl"

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
GLOBAL_FUNCTION char *
an_to_sym(an)
int an;
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
GLOBAL_FUNCTION char *
sym_to_atom(sym)
char *sym;
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

GLOBAL_FUNCTION int
atom_to_an(atom)
char *atom;
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
