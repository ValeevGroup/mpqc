//
// pointgrp.cc
//
// Copyright (C) 1997 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <janssen@limitpt.com>
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/misc/formio.h>
#include <math/symmetry/corrtab.h>

////////////////////////////////////////////////////////////////////////

CorrelationTable::CorrelationTable():
  n_(0),
  ngamma_(0),
  gamma_(0)
{
}

CorrelationTable::CorrelationTable(const RefPointGroup& group,
                                   const RefPointGroup& subgroup):
  n_(0),
  ngamma_(0),
  gamma_(0)
{
  int rc = initialize_table(group,subgroup);
  if (rc == -1) {
      cerr << node0
           << "ERROR: CorrelationTable: too many symop matches" << endl;
      abort();
    }
  else if (rc == -2) {
      cerr << node0
           << "ERROR: CorrelationTable: not a subgroup or wrong ref frame"
           << endl;
      abort();
    }
  else if (rc == -3) {
      cerr << node0 << "ERROR: CorrelationTable: degeneracies don't add up:"
           << endl
           << "  only nondegenerate groups supported"
           << endl;
      abort();
    }
}

CorrelationTable::~CorrelationTable()
{
  clear();
}

int
CorrelationTable::initialize_table(const RefPointGroup& group,
                                   const RefPointGroup& subgroup)
{
  clear();

  group_ = group;
  subgroup_ = subgroup;

  int i, j, k, l;

  CharacterTable ct = group_->char_table();
  CharacterTable subct = subgroup_->char_table();

  n_ = ct.nirrep();
  ngamma_ = new int[n_];
  gamma_ = new int*[n_];

  // CAN ONLY HANDLE NONDEGENERATE POINT GROUPS

  for (i=0; i<n_; i++) {
    ngamma_[i] = 0;
    gamma_[i] = 0;
    }

  // map the ops in the high order to low order groups
  int *so_to_subso = new int[ct.order()];
  int *subso_to_so = new int[subct.order()];
  for (i=0; i<subct.order(); i++) subso_to_so[i] = -1;
  for (i=0; i<ct.order(); i++) {
    SymmetryOperation so = ct.symm_operation(i);
    int found = 0;
    so_to_subso[i] = -1;
    for (j=0; j<subct.order(); j++) {
      SymmetryOperation subso = subct.symm_operation(j);
      double sumsquare = 0.0;
      for (k=0; k<3; k++) {
        for (l=0; l<3; l++) {
          double diff = so(k,l)-subso(k,l);
          sumsquare += diff*diff;
        }
      }
      if (sumsquare < 1.0e-12) {
        found++;
        //cout << scprintf("symmop %d in %s is %d in %s",
        //                 i,ct.symbol(),j,subct.symbol()) << endl;
        so_to_subso[i] = j;
        subso_to_so[j] = i;
        }
      }
    if (found > 1) {
      delete[] so_to_subso;
      delete[] subso_to_so;
      return -1;
      }
    }
  for (i=0; i<subct.order(); i++) {
    if (subso_to_so[i] == -1) {
      delete[] so_to_subso;
      delete[] subso_to_so;
      return -2;
      }
    }

  for (i=0; i<ct.nirrep(); i++) {
    int match;
    for (j=0; j<subct.nirrep(); j++) {
      match = j;
      for (k=0; k<ct.order(); k++) {
        double chr = ct.gamma(i).character(k);
        if (so_to_subso[k] >= 0) {
          double subchr = subct.gamma(j).character(so_to_subso[k]);
          if (fabs(subchr - chr) > 1.0e-6) {
            match = -1;
            break;
            }
          }
        }
      if (match >= 0) {
        int *newgamma = new int[ngamma_[i] + 1];
        memcpy(newgamma,gamma_[i],ngamma_[i]*sizeof(int));
        newgamma[ngamma_[i]] = j;
        ngamma_[i]++;
        delete[] gamma_[i];
        gamma_[i] = newgamma;
        }
      }
    }

  delete[] so_to_subso;
  delete[] subso_to_so;

  for (i=0; i<n(); i++) {
    int degen = ct.gamma(i).degeneracy();
    int subdegen = 0;
    for (j=0; j<ngamma(i); j++) {
      subdegen += subct.gamma(gamma(i,j)).degeneracy();
      }
    if (degen != subdegen) {
      return -3;
      }
    }
}

void
CorrelationTable::clear()
{
  for (int i=0; i<n_; i++) {
    delete[] gamma_[i];
    }
  delete[] ngamma_;
  delete[] gamma_;
}

void
CorrelationTable::print(ostream &o) const
{
  o << node0 << indent
    << "Correlation Table from "
    << group_->symbol() << " to " << subgroup_->symbol() << ":" << endl;

  CharacterTable ct = group_->char_table();
  CharacterTable subct = subgroup_->char_table();

  o << incindent;
  for (int i=0; i<n(); i++) {
    o << node0 << indent
      << ct.gamma(i).symbol() << ":";
    for (int j=0; j<ngamma(i); j++) {
      o << node0 << indent
        << " " << subct.gamma(gamma(i,j)).symbol();
      }
    o << endl;
    }
  o << decindent;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
