//
// mp2r12_energy_util.cc
//
// Copyright (C) 2006 Edward Valeev
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

#include <cassert>
#include <chemistry/qc/mbptr12/mp2r12_energy_util.h>
#include <math/mmisc/pairiter.h>

using namespace std;
using namespace sc;

/*******************************
 * class MP2R12EnergyUtil_base *
 *******************************/

MP2R12EnergyUtil_base::MP2R12EnergyUtil_base() {
}

MP2R12EnergyUtil_base::~MP2R12EnergyUtil_base() {
}

MP2R12EnergyUtil_base::MP2R12EnergyUtil_base(const RefSCDimension& oodim,
                     const RefSCDimension& xydim,
                     const RefSCDimension& f12dim,
                     const unsigned int nocc_act) :
oodim_(oodim), xydim_(xydim), f12dim_(f12dim), nf12_(f12dim.n()/xydim.n()), nocc_act_(nocc_act)
{
  gdim_ = new SCDimension(nf12_);
  if (f12dim_.n()%xydim_.n())
    throw ProgrammingError("MP2R12EnergyUtil::MP2R12EnergyUtil -- rank of f12dim must be divisible by rank of xydim",__FILE__,__LINE__);
}

void MP2R12EnergyUtil_base::check_dims(const RefSCMatrix& A) const
{
  const int nrow = A.rowdim().n();
  const int ncol = A.coldim().n();
  if (nrow != f12dim_.n() && nrow != oodim_.n())
    throw ProgrammingError("MP2R12EnergyUtil::check_dims -- row dimension does not match",__FILE__,__LINE__);
  if (ncol != f12dim_.n() && ncol != oodim_.n())
    throw ProgrammingError("MP2R12EnergyUtil::check_dims -- column dimension does not match",__FILE__,__LINE__);
}

void MP2R12EnergyUtil_base::check_dims(const RefSymmSCMatrix& A) const
{
  const int n = A.dim().n();
  if (n != f12dim_.n())
    throw ProgrammingError("MP2R12EnergyUtil::check_dims -- dimension does not match",__FILE__,__LINE__);
}

////

/*******************************
 * class MP2R12EnergyUtil_Diag *
 *******************************/

MP2R12EnergyUtil_Diag::MP2R12EnergyUtil_Diag() {
}

MP2R12EnergyUtil_Diag::~MP2R12EnergyUtil_Diag() {
}

MP2R12EnergyUtil_Diag::MP2R12EnergyUtil_Diag(const RefSCDimension& oodim,
                     const RefSCDimension& xydim,
                     const RefSCDimension& f12dim,
                     const unsigned int nocc_act) :
oodim_(oodim), xydim_(xydim), f12dim_(f12dim), nf12_(f12dim.n()/xydim.n()), nocc_act_(nocc_act)
{
  gdim_ = new SCDimension(nf12_);
  if (f12dim_.n()%xydim_.n())
    throw ProgrammingError("MP2R12EnergyUtil::MP2R12EnergyUtil -- rank of f12dim must be divisible by rank of xydim",__FILE__,__LINE__);
}

void MP2R12EnergyUtil_Diag::check_dims(const RefSCMatrix& A) const
{
  const int nrow = A.rowdim().n();
  const int ncol = A.coldim().n();
  if (nrow != f12dim_.n() && nrow != oodim_.n())
    throw ProgrammingError("MP2R12EnergyUtil::check_dims -- row dimension does not match",__FILE__,__LINE__);
  if (ncol != f12dim_.n() && ncol != oodim_.n())
    throw ProgrammingError("MP2R12EnergyUtil::check_dims -- column dimension does not match",__FILE__,__LINE__);
}

void MP2R12EnergyUtil_Diag::check_dims(const RefSymmSCMatrix& A) const
{
  const int n = A.dim().n();
  if (n != f12dim_.n())
    throw ProgrammingError("MP2R12EnergyUtil::check_dims -- dimension does not match",__FILE__,__LINE__);
}


namespace sc {

  /**********************************
   * class MP2R12EnergyUtil_Nondiag *
   **********************************/

  MP2R12EnergyUtil_Nondiag::MP2R12EnergyUtil_Nondiag(const RefSCDimension& oodim,
                                                     const RefSCDimension& xydim,
                                                     const RefSCDimension& f12dim,
                                                     const unsigned int nocc_act) :
                                                       MP2R12EnergyUtil_base(oodim,xydim,f12dim,nocc_act) {}

  unsigned int MP2R12EnergyUtil_Nondiag::nrowblks(const RefSCMatrix& A) const {
    throw ProgrammingError("MP2R12EnergyUtil_Nondiag::nrowblks -- should not be used when Diag=false",__FILE__,__LINE__);
  }


  unsigned int MP2R12EnergyUtil_Nondiag::ncolblks(const RefSCMatrix& A) const {
    throw ProgrammingError("MP2R12EnergyUtil_Nondiag::ncolblks -- should not be used when Diag=false",__FILE__,__LINE__);
  }

  unsigned int MP2R12EnergyUtil_Nondiag::nblks(
                                                         const RefSymmSCMatrix& A) const {
    throw ProgrammingError("MP2R12EnergyUtil_Nondiag::nblks -- should not be used when Diag=false",__FILE__,__LINE__);
  }

  void MP2R12EnergyUtil_Nondiag::invert(RefSymmSCMatrix& A) const {
    A->gen_invert_this();
  }

  RefDiagSCMatrix MP2R12EnergyUtil_Nondiag::eigenvalues(
                                              const RefSymmSCMatrix& A) const {
    return A.eigvals();
  }

  void MP2R12EnergyUtil_Nondiag::diagonalize(
                                             const RefSymmSCMatrix& A,
                                             RefDiagSCMatrix& evals,
                                             RefSCMatrix& evecs) const {
    evals = A.kit()->diagmatrix(A.dim());
    evecs = A.kit()->matrix(A.dim(), A.dim());
    A.diagonalize(evals, evecs);
  }

  void MP2R12EnergyUtil_Nondiag::transform(const RefSymmSCMatrix& B,
                                                     const RefDiagSCMatrix& A,
                                                     const RefSCMatrix& U) const {
    B.assign(0.0);
    B.accumulate_transform(U, A);
  }

  // Solves A*X = B
  void MP2R12EnergyUtil_Nondiag::solve_linear_system(
                                                     const RefSymmSCMatrix& A,
                                                     RefSCMatrix& X,
                                                     const RefSCMatrix& B) const {
    sc::lapack_linsolv_symmnondef(A, X, B);
  }

  void MP2R12EnergyUtil_Nondiag::solve_linear_system(unsigned int ij,
                                                     const RefSymmSCMatrix& A,
                                                     RefSCMatrix& X,
                                                     const RefSCMatrix& B) const {
    RefSCVector Xij = X.get_column(ij);
    RefSCVector Bij = B.get_column(ij);
    sc::lapack_linsolv_symmnondef(A,Xij,Bij);
    X.assign_column(Xij,ij);
  }

  // computes y = A x
  void MP2R12EnergyUtil_Nondiag::times(const RefSymmSCMatrix& A,
                                       const RefSCMatrix& x,
                                       RefSCMatrix& y) const {
    y = A*x;
  }

  void MP2R12EnergyUtil_Nondiag::times(unsigned int ij,
                                       const RefSymmSCMatrix& A,
                                       const RefSCMatrix& x,
                                       RefSCMatrix& y) const {
    RefSCVector xij = x.get_column(ij);
    RefSCVector yij = A*xij;
    y.assign_column(yij,ij);
  }

  RefSCVector MP2R12EnergyUtil_Nondiag::dot_product(
                                                    const RefSCMatrix& A,
                                                    const RefSCMatrix& B) const {
    check_dims(A);
    check_dims(B);
    if (A.coldim().n() != B.coldim().n() ||A.rowdim().n() != B.rowdim().n() ||A.coldim().n() != oodim_.n() ||A.rowdim().n() != f12dim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::dot_product -- dimensions do not match",__FILE__,__LINE__);
    RefSCMatrix AB = A.t() * B;
    const int noo = oodim_.n();
    RefSCVector result = AB.kit()->vector(oodim_);
    for (int ij=0; ij<noo; ij++)
      result(ij) = AB.get_element(ij,ij);
    return result;
  }

  void MP2R12EnergyUtil_Nondiag::print(const char* label,
                                                 const RefSCMatrix& A,
                                                 std::ostream& os) const {
    A.print(label, os);
  }

  void MP2R12EnergyUtil_Nondiag::print(const char* label,
                                                 const RefSymmSCMatrix& A,
                                                 std::ostream& os) const {
    A.print(label, os);
  }

  void MP2R12EnergyUtil_Nondiag::print(const char* label,
                                                 const RefDiagSCMatrix& A,
                                                 std::ostream& os) const {
    A.print(label, os);
  }

  /*********************************************
   * class MP2R12EnergyUtil_Diag_DifferentSpin *
   *********************************************/

  MP2R12EnergyUtil_Diag_DifferentSpin::MP2R12EnergyUtil_Diag_DifferentSpin(const RefSCDimension& oodim,
                                                                             const RefSCDimension& xydim,
                                                                             const RefSCDimension& f12dim,
                                                                             const unsigned int nocc_act) :
    MP2R12EnergyUtil_Diag(oodim,xydim,f12dim,nocc_act) {
    gdim_ = new SCDimension(nf12_);
    if (f12dim_.n()%xydim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::MP2R12EnergyUtil -- rank of f12dim must be divisible by rank of xydim",__FILE__,__LINE__);
    if (oodim_.n() != xydim_.n())
      throw ProgrammingError("MP2R12EnergyUtil_Diag_DifferentSpin::MP2R12EnergyUtil_Diag_DifferentSpin -- number of generating pairs must be as nij if diagonal ansatz is chosen",__FILE__,__LINE__);
  }

  // number of blocks should only be needed for diagonal ansatze
  unsigned int MP2R12EnergyUtil_Diag_DifferentSpin::nrowblks(const RefSCMatrix& A) const {
    check_dims(A);
    return A.rowdim().n()/oodim_.n();
  }

  // number of blocks should only be needed for diagonal ansatze
  unsigned int MP2R12EnergyUtil_Diag_DifferentSpin::ncolblks(const RefSCMatrix& A) const {
    check_dims(A);
    return A.coldim().n()/oodim_.n();
  }

  // number of blocks should only be needed for diagonal ansatze
  unsigned int MP2R12EnergyUtil_Diag_DifferentSpin::nblks(const RefSymmSCMatrix& A) const {
    check_dims(A);
    return A.dim().n()/oodim_.n();
  }

  void MP2R12EnergyUtil_Diag_DifferentSpin::get(unsigned int ij, const RefSCMatrix& A, const RefSCVector& Aij) const {
    const unsigned int nij = oodim_.n();
    const int nrow = A.rowdim().n();
    const int ncol = A.coldim().n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::get -- ij >= nij",__FILE__,__LINE__);
    if (nrow != f12dim_.n() || ncol != oodim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::get -- incorrect dimensions of A",__FILE__,__LINE__);
    int i=ij/nocc_act_;
    int j=ij%nocc_act_;
    int ji;
    if(i!=j){
      ji=j*nocc_act_+i;
      for (unsigned int g=0; g<nf12_; g++) {
        int goffset=2*g;
        Aij.set_element(goffset  ,A.get_element(g*nij+ij, ij));
        Aij.set_element(goffset+1,A.get_element(g*nij+ji, ij));
      }
    }  // i!=j
    else {  // i==j
      for (unsigned int g=0; g<nf12_; g++) {
        Aij.set_element(g, A.get_element(g*nij+ij, ij));
      }
    }  // i==j
  }

  void MP2R12EnergyUtil_Diag_DifferentSpin::get(unsigned int ij, const RefSCMatrix& A, const RefSCMatrix& Aij) const {
    const unsigned int nij = oodim_.n();
    const int nrow = A.rowdim().n();
    const int ncol = A.coldim().n();
    const int nrb = nrowblks(A);
    const int ncb = ncolblks(A);
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::get -- ij >= nij",__FILE__,__LINE__);
    if (nrow != f12dim_.n() || ncol != oodim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::get -- incorrect dimensions of A",__FILE__,__LINE__);
    if(nrb!=ncb){
      throw ProgrammingError("MP2R12EnergyUtil_Diag::get -- not implemented for case of different number of column and row blocks",__FILE__,__LINE__);
    }
    int i=ij/nocc_act_;
    int j=ij%nocc_act_;
    int ji;
    if(i!=j) {
      ji=j*nocc_act_+i;
      unsigned int gij = 0;
      for (unsigned int g=0; g<nrb; ++g, gij+=nij) {
        int goffset=2*g;
        unsigned int fij = 0;
        for (unsigned int f=0; f<ncb; ++f, fij+=nij) {
          int foffset=2*f;
          Aij.set_element(goffset  ,foffset  ,A.get_element(gij+ij, fij+ij));
          Aij.set_element(goffset  ,foffset+1,A.get_element(gij+ij, fij+ji));
          Aij.set_element(goffset+1,foffset  ,A.get_element(gij+ji, fij+ij));
          Aij.set_element(goffset+1,foffset+1,A.get_element(gij+ji, fij+ji));
        }
      }
    }  // i!=j
    else {  // i==j
      unsigned int gij = 0;
      for (unsigned int g=0; g<nrb; ++g, gij+=nij) {
        unsigned int fij = 0;
        for (unsigned int f=0; f<ncb; ++f, fij+=nij) {
          //Aij.set_element(g, f, A.get_element(gij, fij));
          Aij.set_element(g, f, A.get_element(gij+ij, fij+ij));
        }
      }
    }  // i==j
  }

  void MP2R12EnergyUtil_Diag_DifferentSpin::get(unsigned int ij, const RefSymmSCMatrix& A, const RefSymmSCMatrix& Aij) const {
    const unsigned int nij = oodim_.n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::get -- ij >= nij",__FILE__,__LINE__);
    if (A.dim().n() != f12dim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::get -- dimension of A does not match",__FILE__,__LINE__);
    int i,j,ji;
    i=ij/nocc_act_;
    j=ij%nocc_act_;
    if(i!=j) {
      ji=j*nocc_act_+i;
      for (unsigned int g=0; g<nf12_; g++) {
        int goffset=2*g;
        int gij=g*nij;
        for (unsigned int f=0; f<=g; f++) {
          int foffset=2*f;
          int fij=f*nij;
          Aij.set_element(goffset  ,foffset  ,A.get_element(gij+ij, fij+ij));
          Aij.set_element(goffset  ,foffset+1,A.get_element(gij+ij, fij+ji));
          Aij.set_element(goffset+1,foffset  ,A.get_element(gij+ji, fij+ij));
          Aij.set_element(goffset+1,foffset+1,A.get_element(gij+ji, fij+ji));
        }
      }
    }  // i!=j
    else{  // i==j
      for (unsigned int g=0; g<nf12_; g++) {
        for (unsigned int f=0; f<=g; f++) {
          Aij.set_element(g, f, A.get_element(g*nij+ij, f*nij+ij));
        }
      }
    }  // i==j
  }

  void MP2R12EnergyUtil_Diag_DifferentSpin::get(unsigned int ij, const RefDiagSCMatrix& A, const RefDiagSCMatrix& Aij) const {
    const unsigned int nij = oodim_.n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::get -- ij >= nij",__FILE__,__LINE__);
    if (A.dim().n() != f12dim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::get -- dimension of A does not match",__FILE__,__LINE__);
    unsigned int gij;
    int i,j,ji;
    i=ij/nocc_act_;
    j=ij%nocc_act_;
    if(i!=j) {
      ji=j*nocc_act_+i;
      gij=0;
      for (unsigned int g=0; g<nf12_; g++, gij+=nij) {
        int goffset=2*g;
        Aij.set_element(goffset  ,A.get_element(gij+ij));
        Aij.set_element(goffset+1,A.get_element(gij+ji));
      }
    }  // i!=j
    else {  // i==j
      gij=0;
      for (unsigned int g=0; g<nf12_; g++, gij+=nij) {
        Aij.set_element(g, A.get_element(gij+ij));
      }
    }  // i==j
  }

  void MP2R12EnergyUtil_Diag_DifferentSpin::put(unsigned int ij, const RefSCMatrix& A, const RefSCVector& Aij) const {
    const unsigned int nij = oodim_.n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::put -- ij >= nij",__FILE__,__LINE__);
    const int nrow = A.rowdim().n();
    const int ncol = A.coldim().n();
    if (nrow != f12dim_.n() || ncol != oodim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::put -- incorrect dimensions of A",__FILE__,__LINE__);
    int i,j,ji;
    i=ij/nocc_act_;
    j=ij%nocc_act_;
    ji=j*nocc_act_+i;
    if(i!=j) {
      for (unsigned int g=0; g<nf12_; g++) {
        int goffset=2*g;
        A.set_element(g*nij+ij, ij, Aij.get_element(goffset));
        A.set_element(g*nij+ji, ij, Aij.get_element(goffset+1));
      }
    }  // i!=j
    else {  // i==j
      for (unsigned int g=0; g<nf12_; g++) {
        A.set_element(g*nij+ij, ij, Aij.get_element(g));
      }
    }  // i==j
  }

  void MP2R12EnergyUtil_Diag_DifferentSpin::put(unsigned int ij, const RefSCMatrix& A, const RefSCMatrix& Aij) const {
    const unsigned int nij = oodim_.n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::put -- ij >= nij",__FILE__,__LINE__);
    const unsigned int nrb = nrowblks(A);
    const unsigned int ncb = ncolblks(A);
    if(nrb!=ncb){
      throw ProgrammingError("MP2R12EnergyUtil_Diag::put -- only implemented for the case of same number of row and column blocks.",__FILE__,__LINE__);
    }
    int i,j,ji;
    i=ij/nocc_act_;
    j=ij%nocc_act_;
    if(i!=j) {
      ji=j*nocc_act_+i;
      unsigned int gij = 0;
      for (unsigned int g=0; g<nrb; ++g, gij+=nij) {
        int goffset=2*g;
        unsigned int fij = 0;
        for (unsigned int f=0; f<ncb; ++f, fij+=nij) {
          int foffset=2*f;
          A.set_element(gij+ij, fij+ij, Aij.get_element(goffset,foffset));
          A.set_element(gij+ij, fij+ji, Aij.get_element(goffset,foffset+1));
          A.set_element(gij+ji, fij+ij, Aij.get_element(goffset+1,foffset));
          A.set_element(gij+ji, fij+ji, Aij.get_element(goffset+1,foffset+1));
        }
      }
    }
    else {
      unsigned int gij = 0;
      for (unsigned int g=0; g<nrb; ++g, gij+=nij) {
        unsigned int fij = 0;
        for (unsigned int f=0; f<ncb; ++f, fij+=nij) {
          A.set_element(gij+ij, fij+ij, Aij.get_element(g, f));
        }
      }
    }
  }

  void MP2R12EnergyUtil_Diag_DifferentSpin::put(unsigned int ij, const RefSymmSCMatrix& A, const RefSymmSCMatrix& Aij) const {
    const unsigned int nij = oodim_.n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::put -- ij >= nij",__FILE__,__LINE__);
    if (A.dim().n() != f12dim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::put -- dimension of A does not match",__FILE__,__LINE__);
    int i,j,ji;
    i=ij/nocc_act_;
    j=ij%nocc_act_;
    if(i!=j) {  // i!=j
      ji=j*nocc_act_+i;
      for (unsigned int g=0; g<nf12_; g++) {
        int gij=g*nij;
        int goffset=2*g;
        for (unsigned int f=0; f<=g; f++) {
          int fij=f*nij;
          int foffset=2*f;
          A.set_element(gij+ij,fij+ij,Aij.get_element(goffset  ,foffset));
          A.set_element(gij+ij,fij+ji,Aij.get_element(goffset  ,foffset+1));
          A.set_element(gij+ji,fij+ij,Aij.get_element(goffset+1,foffset));
          A.set_element(gij+ji,fij+ji,Aij.get_element(goffset+1,foffset+1));
        }
      }
    }  // i!=j
    else {  // i==j
      for (unsigned int g=0; g<nf12_; g++) {
        for (unsigned int f=0; f<=g; f++) {
          A.set_element(g*nij+ij, f*nij+ij, Aij.get_element(g, f));
        }
      }
    }  // i==j
  }

  void MP2R12EnergyUtil_Diag_DifferentSpin::put(unsigned int ij, const RefDiagSCMatrix& A, const RefDiagSCMatrix& Aij) const {
    const unsigned int nij = oodim_.n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::get -- ij >= nij",__FILE__,__LINE__);
    if (A.dim().n() != f12dim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::get -- dimension of A does not match",__FILE__,__LINE__);
    int i,j,ji;
    i=ij/nocc_act_;
    j=ij%nocc_act_;
    if(i!=j){
      ji=j*nocc_act_+i;
      unsigned int gij = 0;
      for (unsigned int g=0; g<nf12_; g++, gij+=nij) {
        int goffset=2*g;
        A.set_element(gij+ij, Aij.get_element(goffset));
        A.set_element(gij+ji, Aij.get_element(goffset+1));
      }
    } // i!=j
    else {  // i==j
      unsigned int gij = 0;
      for (unsigned int g=0; g<nf12_; g++, gij+=nij) {
        A.set_element(gij+ij, Aij.get_element(g));
      }
    }  // i==j
  }

  void MP2R12EnergyUtil_Diag_DifferentSpin::invert(RefSymmSCMatrix& A) const {
    check_dims(A);
    int i,j;
    RefSymmSCMatrix Aij;
    //RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    const int noo = oodim_.n();
    for (int ij=0; ij<noo; ij++) {
      i=ij/nocc_act_;
      j=ij%nocc_act_;
      if(i!=j){
        Aij = A.kit()->symmmatrix(two_gdim);
      }
      else {
        Aij = A.kit()->symmmatrix(gdim_);
      }

      if(i<=j){
        get(ij, A, Aij);
        //Aij.print("B_pair");
        Aij->gen_invert_this();
        put(ij, A, Aij);
      }
    }
  }

  RefDiagSCMatrix MP2R12EnergyUtil_Diag_DifferentSpin::eigenvalues(const RefSymmSCMatrix& A) const {
    check_dims(A);
    RefDiagSCMatrix evals = A.kit()->diagmatrix(f12dim_);
    //RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    int i,j;
    RefSymmSCMatrix Aij;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    const int noo = oodim_.n();
    for (int ij=0; ij<noo; ij++) {
      i=ij/nocc_act_;
      j=ij%nocc_act_;
      if(i!=j){
        Aij = A.kit()->symmmatrix(two_gdim);
      }
      else {
        Aij = A.kit()->symmmatrix(gdim_);
      }

      get(ij, A, Aij);
      RefDiagSCMatrix evals_ij = Aij.eigvals();
      put(ij, evals, evals_ij);
    }

    return evals;
  }

  void MP2R12EnergyUtil_Diag_DifferentSpin::diagonalize(const RefSymmSCMatrix& A, RefDiagSCMatrix& evals, RefSCMatrix& evecs) const {
    check_dims(A);
    evals = A.kit()->diagmatrix(A.dim());
    evecs = A.kit()->matrix(A.dim(), A.dim());

    //RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    RefSymmSCMatrix Aij;
    //RefDiagSCMatrix evals_ij = A.kit()->diagmatrix(gdim_);
    RefDiagSCMatrix evals_ij;
    //RefSCMatrix evecs_ij = A.kit()->matrix(gdim_, gdim_);
    RefSCMatrix evecs_ij;
    int i,j;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    const unsigned int nij = oodim_.n();
    for (unsigned int ij=0; ij<nij; ++ij) {
      i=ij/nocc_act_;
      j=ij%nocc_act_;
      if(i!=j){
        Aij = A.kit()->symmmatrix(two_gdim);
        evals_ij = A.kit()->diagmatrix(two_gdim);
        evecs_ij = A.kit()->matrix(two_gdim, two_gdim);
      }
      else {
        Aij = A.kit()->symmmatrix(gdim_);
        evals_ij = A.kit()->diagmatrix(gdim_);
        evecs_ij = A.kit()->matrix(gdim_, gdim_);
      }

      get(ij, A, Aij);
      Aij.diagonalize(evals_ij, evecs_ij);
      put(ij, evals, evals_ij);
      put(ij, evecs, evecs_ij);
    }
  }

  void MP2R12EnergyUtil_Diag_DifferentSpin::transform(const RefSymmSCMatrix& B, const RefDiagSCMatrix& A, const RefSCMatrix& U) const {
    check_dims(B);
    check_dims(U);
    B.assign(0.0);
    int i,j;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    //RefSymmSCMatrix Bij = B.kit()->symmmatrix(gdim_);
    RefSymmSCMatrix Bij;
    //RefDiagSCMatrix Aij = A.kit()->diagmatrix(gdim_);
    RefDiagSCMatrix Aij;
    //RefSCMatrix Uij = U.kit()->matrix(gdim_, gdim_);
    RefSCMatrix Uij;
    const unsigned int nij = oodim_.n();
    for (unsigned int ij=0; ij<nij; ++ij) {
      i=ij/nocc_act_;
      j=ij%nocc_act_;
      if(i!=j) {
        Bij = B.kit()->symmmatrix(two_gdim);
        Aij = A.kit()->diagmatrix(two_gdim);
        Uij = U.kit()->matrix(two_gdim, two_gdim);
      }
      else {
        Bij = B.kit()->symmmatrix(gdim_);
        Aij = A.kit()->diagmatrix(gdim_);
        Uij = U.kit()->matrix(gdim_, gdim_);
      }

      get(ij, A, Aij);
      get(ij, U, Uij);
      Bij.assign(0.0);
      Bij.accumulate_transform(Uij, Aij);
      put(ij, B, Bij);
    }
  }

  void MP2R12EnergyUtil_Diag_DifferentSpin::solve_linear_system(const RefSymmSCMatrix& A,
                                                                 RefSCMatrix& X,
                                                                 const RefSCMatrix& B) const {
    check_dims(A);
    check_dims(X);
    check_dims(B);
    int i,j;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    //RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    RefSymmSCMatrix Aij;
    //RefSCVector Bij = A.kit()->vector(gdim_);
    RefSCVector Bij;
    //RefSCVector Xij = A.kit()->vector(gdim_);
    RefSCVector Xij;
    const int noo = oodim_.n();
    for (int ij=0; ij<noo; ij++) {
      i=ij/nocc_act_;
      j=ij%nocc_act_;
      if(i!=j) {
        Aij = A.kit()->symmmatrix(two_gdim);
        Bij = A.kit()->vector(two_gdim);
        Xij = A.kit()->vector(two_gdim);
      }
      else {
        Aij = A.kit()->symmmatrix(gdim_);
        Bij = A.kit()->vector(gdim_);
        Xij = A.kit()->vector(gdim_);
      }

      get(ij, A, Aij);
      get(ij, B, Bij);
      sc::lapack_linsolv_symmnondef(Aij, Xij, Bij);
      put(ij, X, Xij);
    }
  }

  void MP2R12EnergyUtil_Diag_DifferentSpin::solve_linear_system(unsigned int ij,
                                                                const RefSymmSCMatrix& A,
                                                                RefSCMatrix& X,
                                                                const RefSCMatrix& B) const {
    check_dims(A);
    check_dims(X);
    check_dims(B);
    int i,j;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    //RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    RefSymmSCMatrix Aij;
    //RefSCVector Bij = A.kit()->vector(gdim_);
    RefSCVector Bij;
    //RefSCVector Xij = A.kit()->vector(gdim_);
    RefSCVector Xij;
    const int noo = oodim_.n();
    i=ij/nocc_act_;
    j=ij%nocc_act_;
    if(i!=j) {
      Aij = A.kit()->symmmatrix(two_gdim);
      Bij = A.kit()->vector(two_gdim);
      Xij = A.kit()->vector(two_gdim);
    }
    else {
      Aij = A.kit()->symmmatrix(gdim_);
      Bij = A.kit()->vector(gdim_);
      Xij = A.kit()->vector(gdim_);
    }

    get(ij, A, Aij);
    get(ij, B, Bij);
    sc::lapack_linsolv_symmnondef(Aij, Xij, Bij);
    put(ij, X, Xij);
  }

  void MP2R12EnergyUtil_Diag_DifferentSpin::times(const RefSymmSCMatrix& A,
                                                   const RefSCMatrix& x,
                                                   RefSCMatrix& y) const {
    check_dims(A);
    check_dims(x);
    check_dims(y);
    int i,j;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    //RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    RefSymmSCMatrix Aij;
    //RefSCVector xij = A.kit()->vector(gdim_);
    RefSCVector xij;
    //RefSCVector yij = A.kit()->vector(gdim_);
    RefSCVector yij;
    const int noo = oodim_.n();
    for (int ij=0; ij<noo; ij++) {
      i=ij/nocc_act_;
      j=ij%nocc_act_;
      if(i!=j) {
        Aij = A.kit()->symmmatrix(two_gdim);
        xij = A.kit()->vector(two_gdim);
        yij = A.kit()->vector(two_gdim);
      }
      else {
        Aij = A.kit()->symmmatrix(gdim_);
        xij = A.kit()->vector(gdim_);
        yij = A.kit()->vector(gdim_);
      }

      get(ij, A, Aij);
      get(ij, x, xij);
      yij = Aij * xij;
      put(ij, y, yij);
    }
  }

  void MP2R12EnergyUtil_Diag_DifferentSpin::times(unsigned int ij,
                                                  const RefSymmSCMatrix& A,
                                                  const RefSCMatrix& x,
                                                  RefSCMatrix& y) const {
    check_dims(A);
    check_dims(x);
    check_dims(y);
    int i,j;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    //RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    RefSymmSCMatrix Aij;
    //RefSCVector xij = A.kit()->vector(gdim_);
    RefSCVector xij;
    //RefSCVector yij = A.kit()->vector(gdim_);
    RefSCVector yij;
    const int noo = oodim_.n();
    i=ij/nocc_act_;
    j=ij%nocc_act_;
    if(i!=j) {
      Aij = A.kit()->symmmatrix(two_gdim);
      xij = A.kit()->vector(two_gdim);
      yij = A.kit()->vector(two_gdim);
    }
    else {
      Aij = A.kit()->symmmatrix(gdim_);
      xij = A.kit()->vector(gdim_);
      yij = A.kit()->vector(gdim_);
    }

    get(ij, A, Aij);
    get(ij, x, xij);
    yij = Aij * xij;
    put(ij, y, yij);
  }

  RefSCVector MP2R12EnergyUtil_Diag_DifferentSpin::dot_product(const RefSCMatrix& A,
                                                                const RefSCMatrix& B) const {
    check_dims(A);
    check_dims(B);
    if (A.coldim().n() != B.coldim().n() ||A.rowdim().n() != B.rowdim().n() ||A.coldim().n() != oodim_.n() ||A.rowdim().n() != f12dim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::dot_product -- dimentions do not match",__FILE__,__LINE__);
    const int noo = oodim_.n();
    RefSCVector result = A.kit()->vector(oodim_);
    int i,j;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    //RefSCVector Aij = A.kit()->vector(gdim_);
    RefSCVector Aij;
    //RefSCVector Bij = A.kit()->vector(gdim_);
    RefSCVector Bij;
    for (int ij=0; ij<noo; ij++) {
      i=ij/nocc_act_;
      j=ij%nocc_act_;
      if(i!=j){
        Aij = A.kit()->vector(two_gdim);
        Bij = A.kit()->vector(two_gdim);
      }
      else {
        Aij = A.kit()->vector(gdim_);
        Bij = A.kit()->vector(gdim_);
      }

      get(ij, A, Aij);
      get(ij, B, Bij);
      result(ij) = Aij.dot(Bij);
    }

    return result;
  }

  void MP2R12EnergyUtil_Diag_DifferentSpin::print(const char* label, const RefSCMatrix& A, std::ostream& os) const {
    os << indent << label << ":"<< endl;
    int i,j;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    //RefSCVector Aij = A.kit()->vector(gdim_);
    RefSCVector Aij;
    const unsigned int noo = oodim_.n();
    for (unsigned int ij=0; ij<noo; ++ij) {
      i=ij/nocc_act_;
      j=ij%nocc_act_;
      if(i!=j){
        Aij = A.kit()->vector(two_gdim);
      }
      else {
        Aij = A.kit()->vector(gdim_);
      }

      get(ij, A, Aij);
      ostringstream oss;
      oss << "Block "<< ij;
      Aij.print(oss.str().c_str(), os);
    }
  }

  void MP2R12EnergyUtil_Diag_DifferentSpin::print(const char* label, const RefSymmSCMatrix& A, std::ostream& os) const {
    os << indent << label << ":"<< endl;
    int i,j;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    //RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    RefSymmSCMatrix Aij;
    const unsigned int noo = oodim_.n();
    for (unsigned int ij=0; ij<noo; ++ij) {
      i=ij/nocc_act_;
      j=ij%nocc_act_;
      if(i!=j) {
        Aij = A.kit()->symmmatrix(two_gdim);
      }
      else {
        Aij = A.kit()->symmmatrix(gdim_);
      }

      get(ij, A, Aij);
      ostringstream oss;
      oss << "Block "<< ij;
      Aij.print(oss.str().c_str(), os);
    }
  }

  void MP2R12EnergyUtil_Diag_DifferentSpin::print(const char* label, const RefDiagSCMatrix& A, std::ostream& os) const {
    os << indent << label << ":"<< endl;
    int i,j;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    //RefDiagSCMatrix Aij = A.kit()->diagmatrix(gdim_);
    RefDiagSCMatrix Aij;
    const unsigned int noo = oodim_.n();
    for (unsigned int ij=0; ij<noo; ++ij) {
      i=ij/nocc_act_;
      j=ij%nocc_act_;
      if(i!=j) {
        Aij = A.kit()->diagmatrix(two_gdim);
      }
      else {
        Aij = A.kit()->diagmatrix(gdim_);
      }

      get(ij, A, Aij);
      ostringstream oss;
      oss << "Block "<< ij;
      Aij.print(oss.str().c_str(), os);
    }
  }

  /****************************************
   * class MP2R12EnergyUtil_Diag_SameSpin *
   ****************************************/

  MP2R12EnergyUtil_Diag_SameSpin::MP2R12EnergyUtil_Diag_SameSpin(const RefSCDimension& oodim,
                                                                   const RefSCDimension& xydim,
                                                                   const RefSCDimension& f12dim,
                                                                   const unsigned int nocc_act) :
    MP2R12EnergyUtil_Diag(oodim,xydim,f12dim,nocc_act) {
    gdim_ = new SCDimension(nf12_);
    if (f12dim_.n()%xydim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::MP2R12EnergyUtil -- rank of f12dim must be divisible by rank of xydim",__FILE__,__LINE__);
    if (oodim_.n() != xydim_.n())
      throw ProgrammingError("MP2R12EnergyUtil_Diag_SameSpin::MP2R12EnergyUtil_Diag_SameSpin -- number of generating pairs must be as nij if diagonal ansatz is chosen",__FILE__,__LINE__);
  }

  // number of blocks should only be needed for diagonal ansatze
  unsigned int MP2R12EnergyUtil_Diag_SameSpin::nrowblks(const RefSCMatrix& A) const {
    check_dims(A);
    return A.rowdim().n()/oodim_.n();
  }

  // number of blocks should only be needed for diagonal ansatze
  unsigned int MP2R12EnergyUtil_Diag_SameSpin::ncolblks(const RefSCMatrix& A) const {
    check_dims(A);
    return A.coldim().n()/oodim_.n();
  }

  // number of blocks should only be needed for diagonal ansatze
  unsigned int MP2R12EnergyUtil_Diag_SameSpin::nblks(const RefSymmSCMatrix& A) const {
    check_dims(A);
    return A.dim().n()/oodim_.n();
  }

  void MP2R12EnergyUtil_Diag_SameSpin::get(unsigned int ij, const RefSCMatrix& A, const RefSCVector& Aij) const {
    const unsigned int nij = oodim_.n();
    const int nrow = A.rowdim().n();
    const int ncol = A.coldim().n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::get -- ij >= nij",__FILE__,__LINE__);
    if (nrow != f12dim_.n() || ncol != oodim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::get -- incorrect dimensions of A",__FILE__,__LINE__);
    for (unsigned int g=0; g<nf12_; g++) {
      Aij.set_element(g, A.get_element(g*nij+ij, ij));
    }
  }

  void MP2R12EnergyUtil_Diag_SameSpin::get(unsigned int ij, const RefSCMatrix& A, const RefSCMatrix& Aij) const {
    const unsigned int nij = oodim_.n();
    const int nrow = A.rowdim().n();
    const int ncol = A.coldim().n();
    const int nrb = nrowblks(A);
    const int ncb = ncolblks(A);
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::get -- ij >= nij",__FILE__,__LINE__);
    if (nrow != f12dim_.n() || ncol != oodim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::get -- incorrect dimensions of A",__FILE__,__LINE__);
    if(nrb!=ncb){
      throw ProgrammingError("MP2R12EnergyUtil_Diag::get -- not implemented for case of different number of column and row blocks",__FILE__,__LINE__);
    }
    unsigned int gij = 0;
    for (unsigned int g=0; g<nrb; ++g, gij+=nij) {
      unsigned int fij = 0;
      for (unsigned int f=0; f<ncb; ++f, fij+=nij) {
        //Aij.set_element(g, f, A.get_element(gij, fij));
        Aij.set_element(g, f, A.get_element(gij+ij, fij+ij));
      }
    }
  }

  void MP2R12EnergyUtil_Diag_SameSpin::get(unsigned int ij, const RefSymmSCMatrix& A, const RefSymmSCMatrix& Aij) const {
    const unsigned int nij = oodim_.n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::get -- ij >= nij",__FILE__,__LINE__);
    if (A.dim().n() != f12dim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::get -- dimension of A does not match",__FILE__,__LINE__);
    int i,j;
    for (unsigned int g=0; g<nf12_; g++) {
      for (unsigned int f=0; f<=g; f++) {
        Aij.set_element(g, f, A.get_element(g*nij+ij, f*nij+ij));
      }
    }
  }

  void MP2R12EnergyUtil_Diag_SameSpin::get(unsigned int ij, const RefDiagSCMatrix& A, const RefDiagSCMatrix& Aij) const {
    const unsigned int nij = oodim_.n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::get -- ij >= nij",__FILE__,__LINE__);
    if (A.dim().n() != f12dim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::get -- dimension of A does not match",__FILE__,__LINE__);
    unsigned int gij;
    int i,j;
    gij=0;
    for (unsigned int g=0; g<nf12_; g++, gij+=nij) {
      Aij.set_element(g, A.get_element(gij+ij));
    }
  }

  void MP2R12EnergyUtil_Diag_SameSpin::put(unsigned int ij, const RefSCMatrix& A, const RefSCVector& Aij) const {
    const unsigned int nij = oodim_.n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::put -- ij >= nij",__FILE__,__LINE__);
    const int nrow = A.rowdim().n();
    const int ncol = A.coldim().n();
    if (nrow != f12dim_.n() || ncol != oodim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::put -- incorrect dimensions of A",__FILE__,__LINE__);
    int i,j;
    for (unsigned int g=0; g<nf12_; g++) {
      A.set_element(g*nij+ij, ij, Aij.get_element(g));
    }
  }

  void MP2R12EnergyUtil_Diag_SameSpin::put(unsigned int ij, const RefSCMatrix& A, const RefSCMatrix& Aij) const {
    const unsigned int nij = oodim_.n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::put -- ij >= nij",__FILE__,__LINE__);
    const unsigned int nrb = nrowblks(A);
    const unsigned int ncb = ncolblks(A);
    if(nrb!=ncb){
      throw ProgrammingError("MP2R12EnergyUtil_Diag::put -- only implemented for the case of same number of row and column blocks.",__FILE__,__LINE__);
    }
    int i,j;
    unsigned int gij = 0;
    for (unsigned int g=0; g<nrb; ++g, gij+=nij) {
      unsigned int fij = 0;
      for (unsigned int f=0; f<ncb; ++f, fij+=nij) {
        A.set_element(gij+ij, fij+ij, Aij.get_element(g, f));
      }
    }
  }

  void MP2R12EnergyUtil_Diag_SameSpin::put(unsigned int ij, const RefSymmSCMatrix& A, const RefSymmSCMatrix& Aij) const {
    const unsigned int nij = oodim_.n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::put -- ij >= nij",__FILE__,__LINE__);
    if (A.dim().n() != f12dim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::put -- dimension of A does not match",__FILE__,__LINE__);
    int i,j;
    for (unsigned int g=0; g<nf12_; g++) {
      for (unsigned int f=0; f<=g; f++) {
        A.set_element(g*nij+ij, f*nij+ij, Aij.get_element(g, f));
      }
    }
  }

  void MP2R12EnergyUtil_Diag_SameSpin::put(unsigned int ij, const RefDiagSCMatrix& A, const RefDiagSCMatrix& Aij) const {
    const unsigned int nij = oodim_.n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::get -- ij >= nij",__FILE__,__LINE__);
    if (A.dim().n() != f12dim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::get -- dimension of A does not match",__FILE__,__LINE__);
    int i,j;
    unsigned int gij = 0;
    for (unsigned int g=0; g<nf12_; g++, gij+=nij) {
      A.set_element(gij+ij, Aij.get_element(g));
    }
  }

  void MP2R12EnergyUtil_Diag_SameSpin::invert(RefSymmSCMatrix& A) const {
    check_dims(A);
    int i,j;
    RefSymmSCMatrix Aij;
    //RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    const int noo = oodim_.n();
    for (int ij=0; ij<noo; ij++) {
      Aij = A.kit()->symmmatrix(gdim_);

      get(ij, A, Aij);
      Aij->gen_invert_this();
      put(ij, A, Aij);
    }
  }

  RefDiagSCMatrix MP2R12EnergyUtil_Diag_SameSpin::eigenvalues(const RefSymmSCMatrix& A) const {
    check_dims(A);
    RefDiagSCMatrix evals = A.kit()->diagmatrix(f12dim_);
    //RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    int i,j;
    RefSymmSCMatrix Aij;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    const int noo = oodim_.n();
    for (int ij=0; ij<noo; ij++) {
      Aij = A.kit()->symmmatrix(gdim_);

      get(ij, A, Aij);
      RefDiagSCMatrix evals_ij = Aij.eigvals();
      put(ij, evals, evals_ij);
    }

    return evals;
  }

  void MP2R12EnergyUtil_Diag_SameSpin::diagonalize(const RefSymmSCMatrix& A, RefDiagSCMatrix& evals, RefSCMatrix& evecs) const {
    check_dims(A);
    evals = A.kit()->diagmatrix(A.dim());
    evecs = A.kit()->matrix(A.dim(), A.dim());

    //RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    RefSymmSCMatrix Aij;
    //RefDiagSCMatrix evals_ij = A.kit()->diagmatrix(gdim_);
    RefDiagSCMatrix evals_ij;
    //RefSCMatrix evecs_ij = A.kit()->matrix(gdim_, gdim_);
    RefSCMatrix evecs_ij;
    int i,j;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    const unsigned int nij = oodim_.n();
    for (unsigned int ij=0; ij<nij; ++ij) {
      Aij = A.kit()->symmmatrix(gdim_);
      evals_ij = A.kit()->diagmatrix(gdim_);
      evecs_ij = A.kit()->matrix(gdim_, gdim_);

      get(ij, A, Aij);
      Aij.diagonalize(evals_ij, evecs_ij);
      put(ij, evals, evals_ij);
      put(ij, evecs, evecs_ij);
    }
  }

  void MP2R12EnergyUtil_Diag_SameSpin::transform(const RefSymmSCMatrix& B, const RefDiagSCMatrix& A, const RefSCMatrix& U) const {
    check_dims(B);
    check_dims(U);
    B.assign(0.0);
    int i,j;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    //RefSymmSCMatrix Bij = B.kit()->symmmatrix(gdim_);
    RefSymmSCMatrix Bij;
    //RefDiagSCMatrix Aij = A.kit()->diagmatrix(gdim_);
    RefDiagSCMatrix Aij;
    //RefSCMatrix Uij = U.kit()->matrix(gdim_, gdim_);
    RefSCMatrix Uij;
    const unsigned int nij = oodim_.n();
    for (unsigned int ij=0; ij<nij; ++ij) {
      Bij = B.kit()->symmmatrix(gdim_);
      Aij = A.kit()->diagmatrix(gdim_);
      Uij = U.kit()->matrix(gdim_, gdim_);

      get(ij, A, Aij);
      get(ij, U, Uij);
      Bij.assign(0.0);
      Bij.accumulate_transform(Uij, Aij);
      put(ij, B, Bij);
    }
  }

  void MP2R12EnergyUtil_Diag_SameSpin::solve_linear_system(const RefSymmSCMatrix& A,
                                                            RefSCMatrix& X,
                                                            const RefSCMatrix& B) const {
    check_dims(A);
    check_dims(X);
    check_dims(B);
    int i,j;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    //RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    RefSymmSCMatrix Aij;
    //RefSCVector Bij = A.kit()->vector(gdim_);
    RefSCVector Bij;
    //RefSCVector Xij = A.kit()->vector(gdim_);
    RefSCVector Xij;
    const int noo = oodim_.n();
    for (int ij=0; ij<noo; ij++) {
      Aij = A.kit()->symmmatrix(gdim_);
      Bij = A.kit()->vector(gdim_);
      Xij = A.kit()->vector(gdim_);

      get(ij, A, Aij);
      get(ij, B, Bij);
      sc::lapack_linsolv_symmnondef(Aij, Xij, Bij);
      put(ij, X, Xij);
    }
  }

  void MP2R12EnergyUtil_Diag_SameSpin::solve_linear_system(unsigned int ij,
                                                           const RefSymmSCMatrix& A,
                                                           RefSCMatrix& X,
                                                           const RefSCMatrix& B) const {
    check_dims(A);
    check_dims(X);
    check_dims(B);
    int i,j;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    //RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    RefSymmSCMatrix Aij;
    //RefSCVector Bij = A.kit()->vector(gdim_);
    RefSCVector Bij;
    //RefSCVector Xij = A.kit()->vector(gdim_);
    RefSCVector Xij;
    const int noo = oodim_.n();
    Aij = A.kit()->symmmatrix(gdim_);
    Bij = A.kit()->vector(gdim_);
    Xij = A.kit()->vector(gdim_);

    get(ij, A, Aij);
    get(ij, B, Bij);
    sc::lapack_linsolv_symmnondef(Aij, Xij, Bij);
    put(ij, X, Xij);
  }

  void MP2R12EnergyUtil_Diag_SameSpin::times(const RefSymmSCMatrix& A,
                                              const RefSCMatrix& x,
                                              RefSCMatrix& y) const {
    check_dims(A);
    check_dims(x);
    check_dims(y);
    int i,j;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    //RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    RefSymmSCMatrix Aij;
    //RefSCVector xij = A.kit()->vector(gdim_);
    RefSCVector xij;
    //RefSCVector yij = A.kit()->vector(gdim_);
    RefSCVector yij;
    const int noo = oodim_.n();
    for (int ij=0; ij<noo; ij++) {
      Aij = A.kit()->symmmatrix(gdim_);
      xij = A.kit()->vector(gdim_);
      yij = A.kit()->vector(gdim_);

      get(ij, A, Aij);
      get(ij, x, xij);
      yij = Aij * xij;
      put(ij, y, yij);
    }
  }

  void MP2R12EnergyUtil_Diag_SameSpin::times(unsigned int ij,
                                             const RefSymmSCMatrix& A,
                                             const RefSCMatrix& x,
                                             RefSCMatrix& y) const {
    check_dims(A);
    check_dims(x);
    check_dims(y);
    int i,j;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    //RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    RefSymmSCMatrix Aij;
    //RefSCVector xij = A.kit()->vector(gdim_);
    RefSCVector xij;
    //RefSCVector yij = A.kit()->vector(gdim_);
    RefSCVector yij;
    const int noo = oodim_.n();
    Aij = A.kit()->symmmatrix(gdim_);
    xij = A.kit()->vector(gdim_);
    yij = A.kit()->vector(gdim_);

    get(ij, A, Aij);
    get(ij, x, xij);
    yij = Aij * xij;
    put(ij, y, yij);
  }

  RefSCVector MP2R12EnergyUtil_Diag_SameSpin::dot_product(const RefSCMatrix& A,
                                                           const RefSCMatrix& B) const {
    check_dims(A);
    check_dims(B);
    if (A.coldim().n() != B.coldim().n() ||A.rowdim().n() != B.rowdim().n() ||A.coldim().n() != oodim_.n() ||A.rowdim().n() != f12dim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::dot_product -- dimentions do not match",__FILE__,__LINE__);
    const int noo = oodim_.n();
    RefSCVector result = A.kit()->vector(oodim_);
    int i,j;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    //RefSCVector Aij = A.kit()->vector(gdim_);
    RefSCVector Aij;
    //RefSCVector Bij = A.kit()->vector(gdim_);
    RefSCVector Bij;
    for (int ij=0; ij<noo; ij++) {
      Aij = A.kit()->vector(gdim_);
      Bij = A.kit()->vector(gdim_);

      get(ij, A, Aij);
      get(ij, B, Bij);
      result(ij) = Aij.dot(Bij);
    }

    return result;
  }

  void MP2R12EnergyUtil_Diag_SameSpin::print(const char* label, const RefSCMatrix& A, std::ostream& os) const {
    os << indent << label << ":"<< endl;
    int i,j;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    //RefSCVector Aij = A.kit()->vector(gdim_);
    RefSCVector Aij;
    const unsigned int noo = oodim_.n();
    for (unsigned int ij=0; ij<noo; ++ij) {
      Aij = A.kit()->vector(gdim_);

      get(ij, A, Aij);
      ostringstream oss;
      oss << "Block "<< ij;
      Aij.print(oss.str().c_str(), os);
    }
  }

  void MP2R12EnergyUtil_Diag_SameSpin::print(const char* label, const RefSymmSCMatrix& A, std::ostream& os) const {
    os << indent << label << ":"<< endl;
    int i,j;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    //RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    RefSymmSCMatrix Aij;
    const unsigned int noo = oodim_.n();
    for (unsigned int ij=0; ij<noo; ++ij) {
      Aij = A.kit()->symmmatrix(gdim_);

      get(ij, A, Aij);
      ostringstream oss;
      oss << "Block "<< ij;
      Aij.print(oss.str().c_str(), os);
    }
  }

  void MP2R12EnergyUtil_Diag_SameSpin::print(const char* label, const RefDiagSCMatrix& A, std::ostream& os) const {
    os << indent << label << ":"<< endl;
    int i,j;
    RefSCDimension two_gdim=new SCDimension(2*gdim_.n());
    //RefDiagSCMatrix Aij = A.kit()->diagmatrix(gdim_);
    RefDiagSCMatrix Aij;
    const unsigned int noo = oodim_.n();
    for (unsigned int ij=0; ij<noo; ++ij) {
      Aij = A.kit()->diagmatrix(gdim_);

      get(ij, A, Aij);
      ostringstream oss;
      oss << "Block "<< ij;
      Aij.print(oss.str().c_str(), os);
    }
  }


  Ref<MP2R12EnergyUtil_Diag> generate_MP2R12EnergyUtil_Diag(SpinCase2 spincase2,
                                                            const RefSCDimension& oodim,
                                                            const RefSCDimension& xydim,
                                                            const RefSCDimension& f12dim,
                                                            const unsigned int nocc_act){
    Ref<MP2R12EnergyUtil_Diag> mp2r12energyutil_diag;
    ExEnv::out0() << "SpinCase2 = " << spincase2 << endl;
    if(spincase2==AlphaBeta) {
      ExEnv::out0() << "generate_MP2R12EnergyUtil_Diag -- generate object of type MP2R12EnergyUtil_Diag_DifferentSpin."
                    << endl;
      mp2r12energyutil_diag = new MP2R12EnergyUtil_Diag_DifferentSpin(oodim,xydim,f12dim,nocc_act);
    }
    else {
      ExEnv::out0() << "generate_MP2R12EnergyUtil_Diag -- generate object of type MP2R12EnergyUtil_Diag_SameSpin."
                    << endl;
      mp2r12energyutil_diag = new MP2R12EnergyUtil_Diag_SameSpin(oodim,xydim,f12dim,nocc_act);
    }

    return(mp2r12energyutil_diag);
  }

  /// fills C with
  void firstorder_cusp_coefficients(const SpinCase2 &spincase2,
                                    RefSCMatrix& C,
                                    const Ref<OrbitalSpace>& i1,
                                    const Ref<OrbitalSpace>& i2,
                                    const Ref<R12Technology::CorrelationFactor>& corrfactor) {
    RefSCDimension xydim = C.coldim();
    RefSCDimension f12dim = C.rowdim();
    const unsigned int nxy = xydim.n();
    const unsigned int nf12 = f12dim.n();
    const unsigned int ngem = nf12/nxy;  // number of geminals
    MPQC_ASSERT(ngem == 1);

    RefSCDimension gdim=new SCDimension(ngem);
    RefSCDimension two_gdim=new SCDimension(2*gdim.n());

    SpinMOPairIter ij_iter(i1->rank(), i2->rank(), spincase2);
    Ref<MP2R12EnergyUtil_Diag> util = generate_MP2R12EnergyUtil_Diag(spincase2,xydim,xydim,f12dim,i2->rank());

    const double Cp_ij_ij = 1.0/2.0;
    const double Cm_ij_ij = 1.0/4.0;

    C.assign(0.0);

    for(ij_iter.start(); int(ij_iter); ij_iter.next()) {
      int ij = ij_iter.ij();
      int i = ij_iter.i();
      int j = ij_iter.j();

      if ((spincase2==AlphaBeta) && (i!=j)){
        RefSCVector C_pair = C.kit()->vector(two_gdim);
        for(int g=0; g<ngem; g++){ // mixed singlet/triplet
          const int goffset=2*g;
          C_pair.set_element(goffset  ,0.5*(Cp_ij_ij+Cm_ij_ij));
          C_pair.set_element(goffset+1,0.5*(Cp_ij_ij-Cm_ij_ij));
        }
        util->put(ij,C,C_pair);
      }
      else {  // (spincase2!=AlphaBeta) || (i==j)
        RefSCVector  C_pair = C.kit()->vector(gdim);
        if(spincase2==AlphaBeta){
          for(int g=0; g<ngem; g++) { // pure singlet
            C_pair.set_element(g,Cp_ij_ij);
          }
          util->put(ij,C,C_pair);
        }
        else {  // spincase2!=AlphaBeta
          for(int g=0; g<ngem; g++) { // pure triplet
            C_pair.set_element(g,Cm_ij_ij);
          }
          util->put(ij,C,C_pair);
        }
      }
    }

  } // end of firstorder_cusp_coefficients()

}

