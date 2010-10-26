//
// atominfo.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <chemistry/qc/intv3/atominfo.gbl>
#include <chemistry/qc/intv3/atominfo.lcl>

using namespace sc;

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
   {0,             0,    0}
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
  return 0;
  }

/* Converts a symbol to an atom name.  If the symbol name is unknown
 * then the symbol name is returned. */
char *
IntV3::sym_to_atom(char *sym)
{
  int i;
  char tmpsym[10];

  if (!sym) return 0;

  /* Convert the passed name to lowercase. */
  strcpy(tmpsym,sym);
  for (i=0; i<strlen(sym); i++) {
    if (isupper(tmpsym[i])) tmpsym[i] += 'a' - 'A';
    }

  for (i=0; atominfo[i].atom != 0; i++) { 
    if (!strcmp(tmpsym,atominfo[i].atom) || !strcmp(tmpsym,atominfo[i].symbol)) { 
      return atominfo[i].atom;
      }
    }
  return 0;
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

  for (i=0; atominfo[i].atom != 0; i++) { 
    if (!strcmp(tmpatom,atominfo[i].atom) || !strcmp(tmpatom,atominfo[i].symbol)) { 
      return atominfo[i].an;
      }
    }

  return 0;
  }

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
