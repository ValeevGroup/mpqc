//
// sobasis.cc --- implementation of the Integral class
//
// Copyright (C) 1998 Limit Point Systems, Inc.
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <string.h>

#include <util/misc/formio.h>
#include <math/symmetry/pointgrp.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/sobasis.h>

using namespace std;
using namespace sc;

SOBasis::SOBasis(const Ref<GaussianBasisSet> &basis, const Ref<Integral>&integral)
{
  int i,j,k;

  basis_ = basis;

  Ref<Molecule> mol = basis_->molecule();

  CharacterTable ct = mol->point_group()->char_table();
  nirrep_ = ct.nirrep();

  // count the number of so shells
  nshell_ = 0;
  for (i=0; i<mol->nunique(); i++) {
    nshell_ += basis_->nshell_on_center(mol->unique(i));
    }

  // map each ao shell to an so shell
  int *aoshell_to_soshell = new int[basis_->nshell()];
  int soshell = 0;
  for (i=0; i<mol->nunique(); i++) {
    for (j=0; j<basis_->nshell_on_center(mol->unique(i)); j++) {
      for (k=0; k<mol->nequivalent(i); k++) {
        int aoshell = basis_->shell_on_center(mol->equivalent(i,k),j);
        aoshell_to_soshell[aoshell] = soshell;
        }
      soshell++;
      }
    }

  ncomp_ = new int[nirrep_];
  for (i=0; i<nirrep_; i++) {
    ncomp_[i] = ct.gamma(i).degeneracy();
    if (ncomp_[i] != 1) {
      ExEnv::out0()
           << "WARNING: SOBasis not tested for degenerate point groups"
           << endl;
      }
    }

  naofunc_ = new int[nshell_];
  memset(naofunc_, 0, sizeof(int)*nshell_);

  nfunc_ = new int*[nshell_];
  funcoff_ = new int*[nshell_];
  for (i=0; i<nshell_; i++) {
    nfunc_[i] = new int[nirrep_];
    funcoff_[i] = new int[nirrep_];
    for (j=0; j<nirrep_; j++) {
      nfunc_[i][j] = 0;
      }
    }

  Ref<PetiteList> petite = new PetiteList(basis_, integral);

  int nblocks = petite->nblocks();
  SO_block *soblocks = petite->aotoso_info();

  trans_ = new SOTransform[nshell_];
  for (i=0; i<nblocks; i++) {
    for (j=0; j<soblocks[i].len; j++) {
      if (soblocks[i].so[j].length == 0) continue;
      int bfn0 = soblocks[i].so[j].cont[0].bfn;
      int aoshell0 = basis_->function_to_shell(bfn0);
      int soshell0 = aoshell_to_soshell[aoshell0];
      int atom0 = basis_->shell_to_center(aoshell0);
      int nequiv0 = mol->nequivalent(mol->atom_to_unique(atom0));
      trans_[soshell0].set_naoshell(nequiv0);
      }
    }

  int nfuncall = 0;
  for (i=0; i<nblocks; i++) {
    int irrep = ct.which_irrep(i);
    for (j=0; j<soblocks[i].len; j++) {
      if (soblocks[i].so[j].length == 0) continue;
      int bfn0 = soblocks[i].so[j].cont[0].bfn;
      int aoshell0 = basis_->function_to_shell(bfn0);
      int soshell0 = aoshell_to_soshell[aoshell0];
      int sofunc = nfunc_[soshell0][irrep];

      int naofunc = basis_->shell(aoshell0).nfunction();
      if (naofunc_[soshell0] && (naofunc_[soshell0] != naofunc)) {
        ExEnv::errn() << "ERROR: SOBasis: mismatch in naofunc" << endl;
        abort();
        }
      naofunc_[soshell0] = naofunc;

      nfunc_[soshell0][irrep]++;
      nfuncall++;

      for (k=0; k<soblocks[i].so[j].length; k++) {
        int bfn = soblocks[i].so[j].cont[k].bfn;
        double coef = soblocks[i].so[j].cont[k].coef;
        int aoshell = basis_->function_to_shell(bfn);
        int aoshellfunc = bfn - basis_->shell_to_function(aoshell);
        int soshell = aoshell_to_soshell[aoshell];

        if (soshell != soshell0) {
          ExEnv::outn() << "ERROR: SOBasis: shell changed" << endl;
          abort();
          }

        trans_[soshell].add_transform(aoshell,irrep, coef,aoshellfunc,sofunc);
        }
      }
    }

  if (nfuncall != basis_->nbasis()) {
    ExEnv::out0() << "ERROR: SOBasis: miscounted number of functions"
         << endl;
    print();
    abort();
    }

  delete[] soblocks;
  delete[] aoshell_to_soshell;

  for (i=0; i<nshell_; i++) {
    funcoff_[i][0] = 0;
    for (j=1; j<nirrep_; j++) {
      funcoff_[i][j] = funcoff_[i][j-1] + nfunc_[i][j-1];
      }
    }

  func_ = new int[nshell_];
  irrep_ = new int[basis_->nbasis()];
  func_within_irrep_ = new int[basis_->nbasis()];
  nfunc_in_irrep_ = new int[nirrep_];

  for (i=0; i<nirrep_; i++) nfunc_in_irrep_[i] = 0;

  if (nshell_) {
    func_[0] = 0;
    for (i=1; i<nshell_; i++) {
      func_[i] = func_[i-1] + nfunction(i-1);
      }
    int ibasis_ = 0;
    for (i=0; i<nshell_; i++) {
      for (j=0; j<nirrep_; j++) {
        for (k=0; k<nfunc_[i][j]; k++,ibasis_++) {
          irrep_[ibasis_] = j;
          func_within_irrep_[ibasis_] = nfunc_in_irrep_[j]++;
          }
        }
      }
    }
}

SOBasis::~SOBasis()
{
  for (int i=0; i<nshell_; i++) {
    delete[] nfunc_[i];
    delete[] funcoff_[i];
    }
  delete[] nfunc_;
  delete[] funcoff_;
  delete[] naofunc_;
  delete[] ncomp_;
  delete[] trans_;
  delete[] func_;
  delete[] irrep_;
  delete[] func_within_irrep_;
  delete[] nfunc_in_irrep_;
}

int
SOBasis::max_nfunction_in_shell() const
{
  int maxn = 0;
  for (int i=0; i<nshell_; i++) {
    int n = nfunction(i);
    if (n > maxn) maxn = n;
    }
  return maxn;
}

int
SOBasis::nfunction(int ishell) const
{
  int n=0;
  for (int i=0; i<nirrep_; i++) {
    n += nfunc_[ishell][i];
    }
  return n;
}

void
SOBasis::print(ostream &o) const
{
  int i,j,k;

  ExEnv::out0()
       << indent << "SOBasis:" << endl
       << incindent
       << basis_
       << indent << "nshell(SO) = " << nshell_ << endl
       << indent << "nirrep = " << nirrep_ << endl;

  ExEnv::out0() << indent << "ncomp = [";
  for (i=0; i<nirrep_; i++) ExEnv::out0() << " " << ncomp_[i];
  ExEnv::out0() << " ]" << endl;

  ExEnv::out0() << indent << "nfunc:" << endl;
  for (i=0; i<nshell_; i++) {
    ExEnv::out0() << indent << "  " << i << ":";
    for (j=0; j<nirrep_; j++) ExEnv::out0() << " " << nfunc_[i][j];
    ExEnv::out0() << endl;
    }

  ExEnv::out0() << indent << "transform:" << endl;
  ExEnv::out0() << incindent;
  for (i=0; i<nshell_; i++) {
    if (i>0) ExEnv::out0() << endl;
    for (j=0; j<trans_[i].naoshell; j++) {
      for (k=0; k<trans_[i].aoshell[j].nfunc; k++) {
        ExEnv::out0() << indent
             << scprintf("SO(%3d %2d %d [%2d]) += % 12.8f * AO(%3d %2d)",
                         i,
                         trans_[i].aoshell[j].func[k].sofunc,
                         trans_[i].aoshell[j].func[k].irrep,
                         function_offset_within_shell(
                           i, trans_[i].aoshell[j].func[k].irrep)
                         + trans_[i].aoshell[j].func[k].sofunc,
                         trans_[i].aoshell[j].func[k].coef,
                         trans_[i].aoshell[j].aoshell,
                         trans_[i].aoshell[j].func[k].aofunc
                         )
             << endl;
        }
      }
    }
  ExEnv::out0() << decindent;

  ExEnv::out0() << decindent;
}

/////////////////////////////////////////////////////////////////////////////

SOTransform::SOTransform()
{
  naoshell_allocated = 0;
  naoshell = 0;
  aoshell = 0;
}

SOTransform::~SOTransform()
{
  delete[] aoshell;
}

void
SOTransform::set_naoshell(int n)
{
  naoshell = 0;
  delete[] aoshell;
  naoshell_allocated = n;
  aoshell = new SOTransformShell[n];
}

void
SOTransform::add_transform(int aoshellnum, int irrep,
                           double coef, int aofunc, int sofunc)
{
  int i;
  for (i=0; i<naoshell; i++) {
    if (aoshell[i].aoshell == aoshellnum) break;
    }
  if (i>=naoshell_allocated) {
    ExEnv::outn() << "ERROR: SOTransform: add_transform allocation too small"
         << endl;
    abort();
    }
  aoshell[i].add_func(irrep,coef,aofunc,sofunc);
  aoshell[i].aoshell = aoshellnum;
  if (i==naoshell) naoshell++;
}

/////////////////////////////////////////////////////////////////////////////

SOTransformShell::SOTransformShell()
{
  nfunc = 0;
  func = 0;
}

SOTransformShell::~SOTransformShell()
{
  delete[] func;
}

void
SOTransformShell::add_func(int irrep, double coef, int aofunc, int sofunc)
{
  SOTransformFunction *newfunc = new SOTransformFunction[nfunc+1];
  for (int i=0; i<nfunc; i++) newfunc[i] = func[i];
  delete[] func;
  func = newfunc;
  func[nfunc].irrep = irrep;
  func[nfunc].coef = coef;
  func[nfunc].aofunc = aofunc;
  func[nfunc].sofunc = sofunc;
  nfunc++;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
