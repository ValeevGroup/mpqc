//
// basisk.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <util/keyval/keyval.h>

extern "C" {
#include <util/sgen/sgen.h>
#include <chemistry/qc/intv2/atoms.h>
#include <chemistry/qc/intv2/atomsinit.h>

#include <chemistry/qc/intv2/atominfo.gbl>
#include <chemistry/qc/intv2/normalize.gbl>
}

#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/basis/gaussshell.h>

enum whats { READ, COUNT };

static int parse_am(char*, shell_type_t&);
static int read_shells(const RefKeyVal&, const char*, const char*,
                       basis_t&, int&,
                       enum whats);

/////////////////////////////////////////////////////////////////////////
//
// the keyval version of the above
// this one actually does the work
// use malloc since the centers struct may be freed by C routines
//

int
int_read_basis(const RefKeyVal& topkeyval,char *atom, const char *bname,
               basis_t &basis)
{
  init_basis(&basis);

  basis.name = strdup(bname);
  if (!basis.name) {
    fprintf(stderr,"read_basis: could not strdup(%s)\n",bname);
    return -1;
  }

  // construct a keyval that contains the basis library

  // this ParsedKeyVal CTOR looks at the basisdir and basisfiles
  // variables to find out what basis set files are to be read in
  RefKeyVal libkeyval = new ParsedKeyVal("basis",topkeyval);
  RefKeyVal aggkeyval = new AggregateKeyVal(topkeyval,libkeyval);
  RefKeyVal keyval = new PrefixKeyVal(":basis",aggkeyval);

  int tmp=0;
  basis.n = read_shells(keyval,atom,bname,basis,tmp,COUNT);
  if (basis.n < 0) {
    fprintf(stderr,"read_basis:  could not count number of shells %s %s\n",
                   atom,bname);
    return -1;
  }

  basis.shell = (shell_t*) malloc(sizeof(shell_t)*basis.n);
  if (!basis.shell) {
    fprintf(stderr,"read_basis:  could not malloc %d shells\n",basis.n);
    return -1;
  }

  tmp=0;
  if (read_shells(keyval,atom,bname,basis,tmp,READ) < 0) {
    fprintf(stderr,"read_basis: trouble in read_shells\n");
    return -1;
  }

  return 0;
}

/////////////////////////////////////////////////////////////////////////
//
// given a keyval, fill in centers.  this is really obsolete and will
// disappear eventually.
//

int
int_read_centers(const RefKeyVal& keyval, centers_t& centers)
{
  char *basis = keyval->pcharvalue("basis");
  if (!basis) {
    fprintf(stderr,"int_read_centers: can't read basis\n");
    return -1;
  }

 // initialize some things
  centers.shell_offset = 0;
  centers.prim_offset = 0;
  centers.func_offset = 0;
  centers.nshell = 0;
  centers.nprim = 0;
  centers.nfunc = 0;
  centers.center_num = 0;
  centers.shell_num = 0;
  centers.func_num = 0;

 // find out how many atoms there are
  centers.n = keyval->count("atoms");
  if (keyval->error() != KeyVal::OK) {
    fprintf(stderr,"int_read_centers: can't count atoms\n");
    return -1;
  }

 // and alloc memory for centers
  centers.center = (center_t*) malloc(sizeof(center_t)*centers.n);
  if (!centers.center) {
    fprintf(stderr,"int_read_centers: could not malloc center\n");
    return -1;
  }

  for (int i=0; i < centers.n; i++) {
    char *atom = keyval->pcharvalue("atoms",i);
    if (keyval->error() != KeyVal::OK) {
      fprintf(stderr,"int_read_centers: can't read atoms[%d]\n",i);
      return -1;
    }

    centers.center[i].atom = strdup(atom);
    if (!centers.center[i].atom) {
      fprintf(stderr,"int_read_centers: can't strdup atoms[%d]\n",i);
      return -1;
    }
    delete[] atom;

    centers.center[i].charge = keyval->doublevalue("charges",i);
    if (keyval->error() != KeyVal::OK)
      centers.center[i].charge = atom_to_an(centers.center[i].atom);

    centers.center[i].r = (double*) malloc(sizeof(double)*3);
    if (!centers.center[i].atom) {
      fprintf(stderr,"int_read_centers: can't malloc r for atom %d\n",i);
      return -1;
    }

    for (int j=0; j<3; j++) {
      centers.center[i].r[j] = keyval->doublevalue("geometry",i,j);
      if (keyval->error() != KeyVal::OK) {
        fprintf(stderr,"int_read_centers: can't read r(%d,%d)\n",i,j);
        return -1;
      }
    }

    if (int_read_basis(keyval, sym_to_atom(centers.center[i].atom), basis,
                       centers.center[i].basis) < 0) {
      fprintf(stderr,"int_read_centers: can't read basis %d\n",i);
      return -1;
    }
  }

  delete[] basis;

  int_initialize_centers(&centers);

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

static int
read_shells(const RefKeyVal& keyval, const char *atom, const char *bname,
            basis_t& basis, int& nsh, enum whats what)
{

  char key[512];

 // first count what's at :basis:"atom":"bname", this will give us the
 // number of shells + get statements

  sprintf(key,"%s:%s",atom,bname);
  int nelem = keyval->count(key);
  if (keyval->error() != KeyVal::OK) {
    fprintf(stderr,"read_shells: could not count %s\n",key);
    keyval->dump();
    return -1;
  }

 // now loop over the number of elements in :basis:atom:bname

  for (int i=0; i < nelem; i++) {

   // see if :basis:atom:bname:i:type exists, if so, then read in the
   // shell info (if what==READ), and increment nsh by one.

    sprintf(key,"%s:%s:%d:type",atom,bname,i);
    if (keyval->exists(key)) {
      if (what==READ) {
        int j;

       // zero out the bits we aren't responsible for
        basis.shell[nsh].nfunc = 0;
        
       // count type to get the number of general contractions
        int ncon = keyval->count(key);
        if (keyval->error() != KeyVal::OK) {
          fprintf(stderr,"read_shells: could not count %s\n",key);
          return -1;
        }

       // count exp to get the number of primitives
        sprintf(key,"%s:%s:%d:exp",atom,bname,i);
        int nprim = keyval->count(key);
        if (keyval->error() != KeyVal::OK) {
          fprintf(stderr,"read_shells: could not count %s\n",key);
          return -1;
        }

        basis.shell[nsh].nprim = nprim;
        basis.shell[nsh].ncon = ncon;

       // alloc memory for the exp vector, and read it in
        basis.shell[nsh].exp = (double*) malloc(sizeof(double)*nprim);
        if (!basis.shell[nsh].exp) {
          fprintf(stderr,"read_shells: could not malloc exp\n");
          return -1;
        }

        for (j=0; j < nprim; j++) {
          basis.shell[nsh].exp[j] = keyval->doublevalue(key,j);
          if (keyval->error() != KeyVal::OK) {
            fprintf(stderr,"read_shells: could not read %s:%d\n",key,j);
            return -1;
          }
        }

       // alloc memory for the coef matrix, and read it in
        basis.shell[nsh].coef = (double**) malloc(sizeof(double*)*ncon);
        if (!basis.shell[nsh].coef) {
          fprintf(stderr,"read_shells: could not malloc coef\n");
          return -1;
        }

        sprintf(key,"%s:%s:%d:coef",atom,bname,i);
        for (j=0; j < ncon; j++) {
          basis.shell[nsh].coef[j] = (double *) malloc(sizeof(double)*nprim);
          if (!basis.shell[nsh].coef) {
            fprintf(stderr,"read_shells: could not malloc coef[%d]\n",j);
            return -1;
          }

          for (int k=0; k < nprim; k++) {
            basis.shell[nsh].coef[j][k] = keyval->doublevalue(key,j,k);
            if (keyval->error() != KeyVal::OK) {
              fprintf(stderr,"read_shells: could not read %s:%d:%d\n",key,j,k);
              return -1;
            }
          }
        }

       // alloc memory for the type vector and read it in
        basis.shell[nsh].type =
                           (shell_type_t*) malloc(sizeof(shell_type_t)*ncon);
        if (!basis.shell[nsh].type) {
          fprintf(stderr,"read_shells: could not malloc type\n");
          return -1;
        }

        for (j=0; j < ncon; j++) {
          sprintf(key,"%s:%s:%d:type:%d:am",atom,bname,i,j);
          char *am = keyval->pcharvalue(key);
          if (keyval->error() != KeyVal::OK) {
            fprintf(stderr,"read_shells: could not read %s\n",key);
            return -1;
          }

          if (parse_am(am,basis.shell[nsh].type[j]) < 0) {
            fprintf(stderr,"read_shells: could not parse %s\n",am);
            return -1;
          }

          delete[] am;
        }
      }

      nsh++;

    } else {
      sprintf(key,"%s:%s:%d:get",atom,bname,i);
      char *nbas = keyval->pcharvalue(key);
      if (keyval->error() != KeyVal::OK) {
        fprintf(stderr,"read_shells: "
                       "there is neither a shell nor a get here\n");
        return -1;
      }

      if (read_shells(keyval,atom,nbas,basis,nsh,what) < 0) {
        return -1;
      }

      delete[] nbas;
    }
  }

  return nsh;
}

///////////////////////////////////////////////////////////////////////////

static int
parse_am(char *am, shell_type_t& st)
{
  if (!am) {
    fprintf(stderr,"parse_am: am is null\n");
    return -1;
  }

  if (am[0] == 's' || am[0] == 'S') {
    st.am = 0;
  } else if (am[0] == 'p' || am[0] == 'P') {
    st.am = 1;
  } else if (am[0] == 'd' || am[0] == 'D') {
    st.am = 2;
  } else if (am[0] == 'f' || am[0] == 'F') {
    st.am = 3;
  } else if (am[0] == 'g' || am[0] == 'G') {
    st.am = 4;
  } else if (am[0] == 'h' || am[0] == 'H') {
    st.am = 5;
  } else if (am[0] == 'i' || am[0] == 'I') {
    st.am = 6;
  } else if (am[0] == 'j' || am[0] == 'J') {
    st.am = 7;
  } else if (am[0] == 'k' || am[0] == 'K') {
    st.am = 8;
  } else if (am[0] == 'l' || am[0] == 'L') {
    st.am = 9;
  } else if (am[0] == 'm' || am[0] == 'M') {
    st.am = 10;
  } else if (am[0] == 'n' || am[0] == 'N') {
    st.am = 11;
  } else if (am[0] == 'o' || am[0] == 'O') {
    st.am = 12;
  } else if (am[0] == 'q' || am[0] == 'Q') {
    st.am = 13;
  } else if (am[0] == 'r' || am[0] == 'R') {
    st.am = 14;
  } else if (am[0] == 't' || am[0] == 'T') {
    st.am = 15;
  } else if (am[0] == 'u' || am[0] == 'U') {
    st.am = 16;
  } else if (am[0] == 'v' || am[0] == 'V') {
    st.am = 17;
  } else if (am[0] == 'w' || am[0] == 'W') {
    st.am = 18;
  } else if (am[0] == 'x' || am[0] == 'X') {
    st.am = 19;
  } else if (am[0] == 'y' || am[0] == 'Y') {
    st.am = 20;
  } else if (am[0] == 'z' || am[0] == 'Z') {
    st.am = 21;
  }

 // if s, p, or just a single letter, then use cartesian functions
  if (st.am < 2 || am[1] == '\0') {
    st.puream = 0;
  } else {
    int npure = atoi(&am[1]);

    if (npure != 2*st.am + 1) {
      fprintf(stderr,"parse_am:  illegal am %s\n",am);
      return -1;
    }

    st.puream = 1;
  }

  return 0;
}

/////////////////////////////////////////////////////////////////////////
//
// returns a centers struct given a GaussianBasisSet

centers_t*
int_centers_from_gbs(const RefGaussianBasisSet& gbs_)
{
  GaussianBasisSet& gbs = *gbs_.pointer();
  
  centers_t* c;
  c = (centers_t*)malloc(sizeof(centers_t));
  c->n = gbs.ncenter();
  c->center = (center_t*) malloc(sizeof(center_t)*gbs.ncenter());
  c->center_num = 0;
  c->shell_num = 0;
  c->func_num = 0;
  c->nfunc = 0;
  c->nshell = 0;
  c->nprim = 0;
  c->func_offset = 0;
  c->prim_offset = 0;
  c->shell_offset = 0;

  Molecule *mol = gbs.molecule().pointer();
  
  for (int icenter = 0; icenter < gbs.ncenter(); icenter++) {
    if (mol) {
      const Molecule& molr = *mol;
      c->center[icenter].atom = strdup(molr[icenter].element().symbol());
      c->center[icenter].charge = molr[icenter].element().charge();
    } else {
      c->center[icenter].atom = strdup("dummy");
      c->center[icenter].charge = 0.0;
    }

    c->center[icenter].r = (double*)malloc(sizeof(double)*3);
    for (int xyz=0; xyz < 3; xyz++)
      c->center[icenter].r[xyz] = gbs.r(icenter,xyz);

    c->center[icenter].basis.n = gbs.nshell_on_center(icenter);
    c->center[icenter].basis.name = strdup(gbs.name());

    c->center[icenter].basis.shell =
      (shell_t*)malloc(sizeof(shell_t)*gbs.nshell_on_center(icenter));

    shell_t*shell = c->center[icenter].basis.shell;
    int ishell;
    for (ishell = 0; ishell < gbs.nshell_on_center(icenter); ishell++) {
      const GaussianShell &igshell = gbs(icenter, ishell);

      int nprim = shell[ishell].nprim = igshell.nprimitive();
      int ncon = shell[ishell].ncon = igshell.ncontraction();

      shell[ishell].exp = (double *) malloc(sizeof(double)*nprim);
      shell[ishell].type = (shell_type_t *) malloc(sizeof(shell_type_t)*ncon);
      shell[ishell].coef = (double **) malloc(sizeof(double *)*ncon);

      int i;
      for (i=0; i < ncon; i++) {
        shell[ishell].coef[i] = (double *) malloc(sizeof(double)*nprim);
        shell[ishell].type[i].am = igshell.am(i);
        shell[ishell].type[i].puream = igshell.is_pure(i);
      }

      for (i=0; i < nprim; i++) {
        shell[ishell].exp[i] = igshell.exponent(i);
        for (int j=0; j < ncon; j++)
          shell[ishell].coef[j][i] = igshell.coefficient_norm(j,i);
      }
    }
  }

  int_initialize_centers(c);
  return c;
}
