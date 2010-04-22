//
// union.cc
//
// Copyright (C) 2009 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <sstream>
#include <stdexcept>
#include <numeric>
#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <math/scmat/blocked.h>
#include <chemistry/qc/basis/union.h>

using namespace sc;
using namespace std;

ClassDesc
UnionBasisSet::class_desc_(typeid(UnionBasisSet),"UnionBasisSet",2,"public GaussianBasisSet",
  0, create<UnionBasisSet>, create<UnionBasisSet>);

UnionBasisSet::UnionBasisSet(const Ref<KeyVal>& keyval) :
  basis1_(require_dynamic_cast<GaussianBasisSet*>(keyval->describedclassvalue("basis1"),"UnionBasisSet:basis1")),
  basis2_(require_dynamic_cast<GaussianBasisSet*>(keyval->describedclassvalue("basis2"),"UnionBasisSet:basis1"))
{
  sum(basis1_,basis2_);
}

UnionBasisSet::UnionBasisSet(StateIn& si) :
  GaussianBasisSet(si)
{
  basis1_ << SavableState::restore_state(si);
  basis2_ << SavableState::restore_state(si);

  // maps
  si.get(shell_to_basis_);
  si.get(function_to_basis_);
}

void
UnionBasisSet::save_data_state(StateOut& so)
{
  GaussianBasisSet::save_data_state(so);
  SavableState::save_state(basis1_.pointer(),so);
  SavableState::save_state(basis2_.pointer(),so);

  // maps
  so.put(shell_to_basis_);
  so.put(function_to_basis_);
}

UnionBasisSet::UnionBasisSet(const Ref<GaussianBasisSet>& bs1,
                             const Ref<GaussianBasisSet>& bs2) :
  basis1_(bs1), basis2_(bs2)
{
  sum(basis1_,basis2_);
}

void
UnionBasisSet::sum(const Ref<GaussianBasisSet>& A,
                   const Ref<GaussianBasisSet>& B)
{
  GaussianBasisSet* a = A.pointer();
  GaussianBasisSet* b = B.pointer();
  if (a->molecule_.pointer() != b->molecule_.pointer())
    throw std::runtime_error("GaussianBasisSetSum::sum -- cannot sum basis sets, molecules are different");

  Ref<SCMatrixKit> matrixkit = a->matrixkit();
  Ref<Molecule> molecule = a->molecule();
  const int ncenter = a->ncenter();
  std::vector<int> center_to_nshell(ncenter);
  std::vector<int> center_to_shell(ncenter);

  // eliminate duplicate shells to compute the number of shells in the composite basis set
  const int nshell1 = a->nshell();
  const int nshell2 = b->nshell();
  std::vector<bool> shell1_in_bs2(nshell1, false);
  std::vector<bool> shell2_in_bs1(nshell2, false);
  for(int s1=0; s1<nshell1; ++s1) {
    const int c1 = a->shell_to_center(s1);
    const int s1_in_b2 = b->find(c1, a->shell(s1));
    if (s1_in_b2 != -1) // s1 is found in bs2
      shell1_in_bs2[s1] = true;
  }
  int nshell = nshell1;
  for(int s2=0; s2<nshell2; ++s2) {
    const int c2 = b->shell_to_center(s2);
    const int s2_in_b1 = a->find(c2, b->shell(s2));
    if (s2_in_b1 == -1) // s2 is not found in bs1
      ++nshell;
    else
      shell2_in_bs1[s2] = true;
  }

  GaussianShell** shell = new GaussianShell*[nshell];
  int* func_per_shell = new int[nshell];

  // now compute how many shells per each center we have
  // and compute maps from bs12 to bs1 and bs2
  shell_to_basis_.resize(nshell);
  std::vector<unsigned int> shell12_to_shell(nshell);
  for(int ss=0, c=0; c<ncenter; ++c) { // ss is the absolute index in bs12
    const int ns1 = a->nshell_on_center(c);
    const int ns2 = b->nshell_on_center(c);
    const int s1off = a->shell_on_center(c,0);
    const int s2off = b->shell_on_center(c,0);

    // map all shells from bs1 to bs12
    for(int s1=0; s1<ns1; ++s1) {
      shell_to_basis_[ss + s1] = shell1_in_bs2[s1] ? Basis1_and_Basis2 : Basis1;
      shell12_to_shell[ss + s1] = s1+s1off;
    }

    // map shells from bs2 to bs12, skipping shells from bs2 that as also in bs1
    int ns = ns1;
    ss += ns1;
    for(int s2=0; s2<ns2; ++s2)
      if (shell2_in_bs1[s2+s2off] == false) {
        shell_to_basis_[ss] = Basis2;
        shell12_to_shell[ss] = s2 + s2off;
        ++ns;
        ++ss;
      }

    center_to_nshell[c] = ns;
  }

  // compute shell offsets for each center
  center_to_shell[0] = 0;
  for(int c=1; c<ncenter; c++) {
    center_to_shell[c] = center_to_shell[c-1] + center_to_nshell[c-1];
  }

  // copy shells
  std::vector<int> shell_to_center(nshell, -1);
  for(int c=0; c<ncenter; c++) {

    const int ns = center_to_nshell[c];
    const int soff = center_to_shell[c];

    // i is the shell index on this center in the composite basis
    // ii is the absolute index
    for (int i=0; i<ns; i++) {
      const int ii = soff + i;
      const GaussianShell* gsi;
      if (shell_to_basis_[ii] == Basis2)
        gsi = &b->shell( shell12_to_shell[ii] );
      else
        gsi = &a->shell( shell12_to_shell[ii] );

      shell_to_center[ii] = c;

      int nc=gsi->ncontraction();
      int np=gsi->nprimitive();
      func_per_shell[ii] = gsi->nfunction();

      int *ams = new int[nc];
      int *pure = new int[nc];
      double *exps = new double[np];
      double **coefs = new double*[nc];

      for (int j=0; j < nc; j++) {
        ams[j] = gsi->am(j);
        pure[j] = gsi->is_pure(j);
        coefs[j] = new double[np];
        for (int k=0; k < np; k++)
        coefs[j][k] = gsi->coefficient_unnorm(j,k);
      }

      for (int j=0; j < np; j++)
      exps[j] = gsi->exponent(j);

      shell[ii] = new GaussianShell(nc, np, exps, ams, pure, coefs,
          GaussianShell::Unnormalized);
    }
  }

  const int nbasis = std::accumulate(func_per_shell, func_per_shell+nshell, 0);
  RefSCDimension basisdim = new SCDimension(nbasis, nshell, func_per_shell, "basis set dimension");

  const char* A_name = a->name();
  const char* B_name = b->name();
  char* AplusB_name = 0;
  if (A_name && B_name) {
    ostringstream oss;
    oss << "[" << A_name << "]+[" << B_name << "]";
    std::string tmpname = oss.str();
    AplusB_name = strcpy(new char[tmpname.size()+1],tmpname.c_str());
  }
  char* AplusB_label = 0;
  if (AplusB_name) {
    AplusB_label = AplusB_name;
  }
  else {
    ostringstream oss;
    const char* A_label = a->label();
    const char* B_label = b->label();
    oss << "[" << A_label << "]+[" << B_label << "]";
    std::string tmpname = oss.str();
    AplusB_label = strcpy(new char[tmpname.size()+1],tmpname.c_str());
  }

  this->init(AplusB_name, AplusB_label, molecule,
             matrixkit, Ref<SCMatrixKit>(new BlockedSCMatrixKit(matrixkit)),
             shell, shell_to_center);

  delete[] func_per_shell;
  delete[] AplusB_name;
  if (AplusB_name != AplusB_label)
  delete[] AplusB_label;

  //
  // compute the rest of the maps using shell_to_basis_
  //
  function_to_basis_.resize(nbasis);
  for(int s=0; s<nshell; ++s) {
    const int nf = this->shell(s).nfunction();
    const int foff = this->shell_to_function(s);
    int ff = foff;
    for(int f=0; f<nf; ++f, ++ff) {
      function_to_basis_[ff] = shell_to_basis_[s];
    }
  }
}

UnionBasisSet::Basis12
UnionBasisSet::shell_to_basis(int s) const {
  return static_cast<UnionBasisSet::Basis12>(shell_to_basis_[s]);
}

UnionBasisSet::Basis12
UnionBasisSet::function_to_basis(int f) const {
  return static_cast<UnionBasisSet::Basis12>(function_to_basis_[f]);
}


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
