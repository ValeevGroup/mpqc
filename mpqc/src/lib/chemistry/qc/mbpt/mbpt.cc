//
// mbpt.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Ida Nielsen <ida@kemi.aau.dk>
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
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/mbpt/mbpt.h>

/////////////////////////////////////////////////////////////////
// Function dquicksort performs a quick sort (smaller -> larger) 
// of the double data in item by the integer indices in index;
// data in item remain unchanged
/////////////////////////////////////////////////////////////////
static void
dqs(double *item,int *index,int left,int right)
{
  register int i,j;
  double x;
  int y;
 
  i=left; j=right;
  x=item[index[(left+right)/2]];
 
  do {
    while(item[index[i]]<x && i<right) i++;
    while(x<item[index[j]] && j>left) j--;
 
    if (i<=j) {
      if (item[index[i]] != item[index[j]]) {
        y=index[i];
        index[i]=index[j];
        index[j]=y;
        }
      i++; j--;
      }
    } while(i<=j);
       
  if (left<j) dqs(item,index,left,j);
  if (i<right) dqs(item,index,i,right);
}
static void
dquicksort(double *item,int *index,int n)
{
  int i;
  if (n<=0) return;
  for (i=0; i<n; i++) {
    index[i] = i;
    }
  dqs(item,index,0,n-1);
  }

///////////////////////////////////////////////////////////////////////////
// MBPT2

#define CLASSNAME MBPT2
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public Wavefunction
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
MBPT2::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Wavefunction::_castdown(cd);
  return do_castdowns(casts,cd);
}

MBPT2::MBPT2(StateIn& s):
  Wavefunction(s)
  maybe_SavableState(s)
{
  reference_.restore_state(s);
  s.get(nfzc);
  s.get(nfzv);
  s.get(mem_alloc);
  s.getstring(method_);
  s.getstring(algorithm_);
  s.get(debug_);

  eliminate_in_gmat_ = 1;

  mem = MemoryGrp::initial_memorygrp();
  msg_ = MessageGrp::get_default_messagegrp();
}

MBPT2::MBPT2(const RefKeyVal& keyval):
  Wavefunction(keyval)
{
  reference_ = keyval->describedclassvalue("reference");
  if (reference_.null()) {
      cerr << "MBPT2::MBPT2: no reference wavefunction" << endl;
      abort();
    }
  nfzc = keyval->intvalue("nfzc");
  nfzv = keyval->intvalue("nfzv");
  mem_alloc = keyval->intvalue("memory");
  if (keyval->error() != KeyVal::OK) {
      mem_alloc = 8000000;
    }
  mem = keyval->describedclassvalue("memorygrp");
  if (mem.null()) {
      mem = MemoryGrp::initial_memorygrp();
    }
  msg_ = MessageGrp::get_default_messagegrp();

  method_ = keyval->pcharvalue("method");

  algorithm_ = keyval->pcharvalue("algorithm");

  debug_ = keyval->booleanvalue("debug");
  if (keyval->error() == KeyVal::WrongType)
      debug_ = keyval->intvalue("debug");

  eliminate_in_gmat_ = 1;
}

MBPT2::~MBPT2()
{
  delete[] method_;
  delete[] algorithm_;
}

void
MBPT2::save_data_state(StateOut& s)
{
  Wavefunction::save_data_state(s);
  reference_.save_state(s);
  s.put(nfzc);
  s.put(nfzv);
  s.put(mem_alloc);
  s.putstring(method_);
  s.putstring(algorithm_);
  s.put(debug_);
}

void
MBPT2::print(ostream&o)
{
  o << node0 << indent << "MBPT2:" << endl;
  o << incindent;
  Wavefunction::print(o);
  o << node0 << indent << "Reference Wavefunction:" << endl;
  o << incindent; reference_->print(o); o << decindent << endl;
  o << decindent;
}

//////////////////////////////////////////////////////////////////////////////

RefSymmSCMatrix
MBPT2::density()
{
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

void
MBPT2::compute()
{
  init_variables();

  reference_->set_desired_value_accuracy(desired_value_accuracy()/100.);
  if (do_gradient()) {
      if (nsocc) {
          cerr << "MBPT2: cannot compute open shell gradients" << endl;
          abort();
        }
      compute_cs_grad();
    }
  else {
      if (nsocc && algorithm_ && !strcmp(algorithm_,"memgrp")) {
          cerr << "MBPT2: memgrp algorithm cannot compute open shell energy"
               << endl;
          abort();
        }
      if (!nsocc && (!algorithm_ || !strcmp(algorithm_,"memgrp"))) {
          compute_cs_grad();
        }
      else if ((!algorithm_ && msg_->n() <= 32)
               || (algorithm_ && !strcmp(algorithm_,"v1"))) {
          compute_hsos_v1();
        }
      else if (!algorithm_ || !strcmp(algorithm_,"v2")) {
          compute_hsos_v2();
        }
      else if (!strcmp(algorithm_,"v2lb")) {
          compute_hsos_v2_lb();
        }
      else {
          cerr << "MBPT2: unknown algorithm: " << algorithm_ << endl;
          abort();
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

void
MBPT2::obsolete()
{
  reference_->obsolete();
  Compute::obsolete();
}

//////////////////////////////////////////////////////////////////////////////

int
MBPT2::gradient_implemented()
{
  init_variables();

  if (nsocc) return 0;
  return 1;
}

//////////////////////////////////////////////////////////////////////////////

int
MBPT2::value_implemented()
{
  return 1;
}

//////////////////////////////////////////////////////////////////////////////

void
MBPT2::eigen(RefDiagSCMatrix &vals, RefSCMatrix &vecs, RefDiagSCMatrix &occs)
{
  int i;
  if (nsocc) {
      if (reference_->n_fock_matrices() != 2) {
          cerr << "MBPT2: given open reference with"
               << " wrong number of Fock matrices" << endl;
          abort();
        }

      // Notation: oo = orthonormal symmetry orbital basis
      //           ao = atomic orbital basis
      //           so = symmetrized atomic orbital basis
      //           mo1 = SCF molecular orbital basis
      //           mo2 = MBPT molecular orbital basis

      // get the closed shell and open shell AO fock matrices
      RefSymmSCMatrix fock_c_so = reference_->fock(0);
      RefSymmSCMatrix fock_o_so = reference_->fock(1);

      // transform the AO fock matrices to the MO basis
      RefSymmSCMatrix fock_c_mo1 = fock_c_so.clone();
      RefSymmSCMatrix fock_o_mo1 = fock_o_so.clone();
      RefSCMatrix vecs_so_mo1 = reference_->eigenvectors();

      fock_c_mo1.assign(0.0);
      fock_o_mo1.assign(0.0);
      fock_c_mo1.accumulate_transform(vecs_so_mo1.t(), fock_c_so);
      fock_o_mo1.accumulate_transform(vecs_so_mo1.t(), fock_o_so);
      fock_c_so = 0;
      fock_o_so = 0;

     /* Convert to the Guest & Saunders general form.
        This is the form used for an OPT2 calculation.

            C        O         V
        ----------
        |        |
     C  |   fc   |
        |        |
        -------------------
        |        |        |
     O  | 2fc-fo |   fc   |
        |        |        |
        ----------------------------
        |        |        |        |
     V  |   fc   |   fo   |   fc   |
        |        |        |        |
        ----------------------------
      */
      RefSymmSCMatrix fock_eff_mo1 = fock_c_mo1.clone();
      fock_eff_mo1.assign(fock_c_mo1);
      for (i=0; i<nbasis; i++) {
          double occi = reference_->occupation(i);
          for (int j=0; j<=i; j++) {
              double occj = reference_->occupation(j);
              if (occi == 2.0 && occj == 1.0
                  || occi == 1.0 && occj == 2.0) {
                  fock_eff_mo1.accumulate_element(i,j,
                                                  fock_c_mo1(i,j)-fock_o_mo1(i,j));
                }
              else if (occi == 0.0 && occj == 1.0
                  || occi == 1.0 && occj == 0.0) {
                  fock_eff_mo1.accumulate_element(i,j,
                                                  fock_o_mo1(i,j)-fock_c_mo1(i,j));
                }
            }
        }

      // diagonalize the new fock matrix
      RefDiagSCMatrix vals_mo2(fock_eff_mo1.dim(), fock_eff_mo1.kit());
      RefSCMatrix vecs_mo1_mo2(fock_eff_mo1.dim(), fock_eff_mo1.dim(),
                               fock_eff_mo1.kit());
      fock_eff_mo1.diagonalize(vals_mo2, vecs_mo1_mo2);
      vals = vals_mo2;

      // compute the AO to new MO scf vector
      RefSCMatrix so_ao = reference_->integral()->petite_list()->sotoao();
      vecs = vecs_mo1_mo2.t() * vecs_so_mo1.t() * so_ao;
    }
  else {
      // get the closed shell AO fock matrices
      RefSymmSCMatrix fock_c_so = reference_->fock(0);

      // transform the AO fock matrices to the MO basis
      RefSymmSCMatrix fock_c_mo1 = fock_c_so.clone();
      RefSCMatrix vecs_so_mo1 = reference_->eigenvectors();

      fock_c_mo1.assign(0.0);
      fock_c_mo1.accumulate_transform(vecs_so_mo1.t(), fock_c_so);
      fock_c_so = 0;

      // diagonalize the fock matrix
      vals = fock_c_mo1.eigvals();

      // compute the AO to new MO scf vector
      RefSCMatrix so_ao = reference_->integral()->petite_list()->sotoao();
      vecs = vecs_so_mo1.t() * so_ao;
    }
  // fill in the occupations
  occs = matrixkit()->diagmatrix(vals.dim());
  for (i=0; i<nbasis; i++) {
      occs(i) = reference_->occupation(i);
    }
  // sort the eigenvectors and values if symmetry is not c1
  if (molecule()->point_group().char_table().order() != 1) {
      int n = vals.n();
      double *evals = new double[n];
      vals->convert(evals);
      int *indices = new int[n];
      dquicksort(evals,indices,n);
      delete[] evals;
      // make sure all nodes see the same indices and evals
      msg_->bcast(indices,n);
      RefSCMatrix newvecs(vecs.rowdim(), vecs.coldim(), matrixkit());
      RefDiagSCMatrix newvals(vals.dim(), matrixkit());
      RefDiagSCMatrix newoccs(vals.dim(), matrixkit());
      for (i=0; i<n; i++) {
          newvals(i) = vals(indices[i]);
          newoccs(i) = occs(indices[i]);
          for (int j=0; j<n; j++) {
              newvecs(i,j) = vecs(indices[i],j);
            }
        }
      occs = newoccs;
      vecs = newvecs;
      vals = newvals;
      delete[] indices;
    }
}

/////////////////////////////////////////////////////////////////////////////

void
MBPT2::init_variables()
{
  nbasis = basis()->nbasis();
  nocc = nvir = nsocc = 0;
  for (int i=0; i<nbasis; i++) {
    if (reference_->occupation(i) == 2.0) nocc++;
    else if (reference_->occupation(i) == 1.0) nsocc++;
    else nvir++;
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
