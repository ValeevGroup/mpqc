//
// mp2r12_energy_util.cc
//
// Copyright (C) 2006 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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

#include <chemistry/qc/mbptr12/mp2r12_energy_util.h>

using namespace sc;

MP2R12EnergyUtil_base::MP2R12EnergyUtil_base() {
}

MP2R12EnergyUtil_base::~MP2R12EnergyUtil_base() {
}

////

namespace sc {
  
  MP2R12EnergyUtil_Diag::MP2R12EnergyUtil_Diag(const RefSCDimension& oodim,
                       const RefSCDimension& xydim,
                       const RefSCDimension& f12dim) :
  oodim_(oodim), xydim_(xydim), f12dim_(f12dim), nf12_(f12dim.n()/xydim.n())
  {
    gdim_ = new SCDimension(nf12_);
    if (f12dim_.n()%xydim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::MP2R12EnergyUtil -- rank of f12dim must be divisible by rank of xydim",__FILE__,__LINE__);
    if (oodim_.n() != xydim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::MP2R12EnergyUtil -- number of generating pairs must be as nij if diagonal ansatz is chosen",__FILE__,__LINE__);
  }

  MP2R12EnergyUtil_Nondiag::MP2R12EnergyUtil_Nondiag(const RefSCDimension& oodim,
                                                     const RefSCDimension& xydim,
                                                     const RefSCDimension& f12dim) :
  oodim_(oodim), xydim_(xydim), f12dim_(f12dim), nf12_(f12dim.n()/xydim.n())
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

  void MP2R12EnergyUtil_Nondiag::check_dims(const RefSCMatrix& A) const
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

  void MP2R12EnergyUtil_Nondiag::check_dims(const RefSymmSCMatrix& A) const
  {
    const int n = A.dim().n();
    if (n != f12dim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::check_dims -- dimension does not match",__FILE__,__LINE__);
  }
  
  // number of blocks should only be needed for diagonal ansatze
  unsigned int MP2R12EnergyUtil_Diag::nrowblks(const RefSCMatrix& A) const {
    check_dims(A);
    return A.rowdim().n()/oodim_.n();
  }
  
  unsigned int MP2R12EnergyUtil_Diag::ncolblks(const RefSCMatrix& A) const {
    check_dims(A);
    return A.coldim().n()/oodim_.n();
  }
  
  unsigned int MP2R12EnergyUtil_Diag::nblks(const RefSymmSCMatrix& A) const {
    check_dims(A);
    return A.dim().n()/oodim_.n();
  }
  
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
  
  // put/get can only be implemented when Diag=true
  void MP2R12EnergyUtil_Diag::get(unsigned int ij,
                                              const RefSCMatrix& A,
                                              const RefSCVector& Aij) const {
    const unsigned int nij = oodim_.n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::get -- ij >= nij",__FILE__,__LINE__);
    const int nrow = A.rowdim().n();
    const int ncol = A.coldim().n();
    if (nrow != f12dim_.n() || ncol != oodim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::get -- incorrect dimensions of A",__FILE__,__LINE__);
    for (unsigned int g=0; g<nf12_; g++) {
      Aij.set_element(g, A.get_element(g*nij+ij, ij));
    }
  }
  
  void MP2R12EnergyUtil_Diag::get(unsigned int ij,
                                              const RefSCMatrix& A,
                                              const RefSCMatrix& Aij) const {
    const unsigned int nij = oodim_.n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::get -- ij >= nij",__FILE__,__LINE__);
    const int nrb = nrowblks(A);
    const int ncb = ncolblks(A);
    unsigned int gij = 0;
    for (unsigned int g=0; g<nrb; ++g, gij+=nij) {
      unsigned int fij = 0;
      for (unsigned int f=0; f<ncb; ++f, fij+=nij) {
        Aij.set_element(g, f, A.get_element(gij, fij));
      }
    }
  }
  
  void MP2R12EnergyUtil_Diag::get(unsigned int ij,
                                              const RefSymmSCMatrix& A,
                                              const RefSymmSCMatrix& Aij) const {
    const unsigned int nij = oodim_.n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::get -- ij >= nij",__FILE__,__LINE__);
    if (A.dim().n() != f12dim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::get -- dimension of A does not match",__FILE__,__LINE__);
    for (unsigned int g=0; g<nf12_; g++)
      for (unsigned int f=0; f<=g; f++) {
        Aij.set_element(g, f, A.get_element(g*nij+ij, f*nij+ij));
      }
  }
  
  void MP2R12EnergyUtil_Diag::get(unsigned int ij,
                                              const RefDiagSCMatrix& A,
                                              const RefDiagSCMatrix& Aij) const {
    const unsigned int nij = oodim_.n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::get -- ij >= nij",__FILE__,__LINE__);
    if (A.dim().n() != f12dim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::get -- dimension of A does not match",__FILE__,__LINE__);
    unsigned int gij = ij;
    for (unsigned int g=0; g<nf12_; g++, gij+=nij)
      Aij.set_element(g, A.get_element(gij));
  }
  
  void MP2R12EnergyUtil_Diag::put(unsigned int ij,
                                              const RefSCMatrix& A,
                                              const RefSCVector& Aij) const {
    const unsigned int nij = oodim_.n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::put -- ij >= nij",__FILE__,__LINE__);
    const int nrow = A.rowdim().n();
    const int ncol = A.coldim().n();
    if (nrow != f12dim_.n() || ncol != oodim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::put -- incorrect dimensions of A",__FILE__,__LINE__);
    for (unsigned int g=0; g<nf12_; g++) {
      A.set_element(g*nij+ij, ij, Aij.get_element(g));
    }
  }
  
  void MP2R12EnergyUtil_Diag::put(unsigned int ij,
                                              const RefSCMatrix& A,
                                              const RefSCMatrix& Aij) const {
    const unsigned int nij = oodim_.n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::put -- ij >= nij",__FILE__,__LINE__);
    const unsigned int nrb = nrowblks(A);
    const unsigned int ncb = ncolblks(A);
    unsigned int gij = 0;
    for (unsigned int g=0; g<nrb; ++g, gij+=nij) {
      unsigned int fij = 0;
      for (unsigned int f=0; f<ncb; ++f, fij+=nij) {
        A.set_element(gij, fij, Aij.get_element(g, f));
      }
    }
  }
  
  void MP2R12EnergyUtil_Diag::put(unsigned int ij,
                                              const RefSymmSCMatrix& A,
                                              const RefSymmSCMatrix& Aij) const {
    const unsigned int nij = oodim_.n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::put -- ij >= nij",__FILE__,__LINE__);
    if (A.dim().n() != f12dim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::put -- dimension of A does not match",__FILE__,__LINE__);
    for (unsigned int g=0; g<nf12_; g++)
      for (unsigned int f=0; f<=g; f++) {
        A.set_element(g*nij+ij, f*nij+ij, Aij.get_element(g, f));
      }
  }
  
  void MP2R12EnergyUtil_Diag::put(unsigned int ij,
                                              const RefDiagSCMatrix& A,
                                              const RefDiagSCMatrix& Aij) const {
    const unsigned int nij = oodim_.n();
    if (ij >= nij)
      throw ProgrammingError("MP2R12EnergyUtil::get -- ij >= nij",__FILE__,__LINE__);
    if (A.dim().n() != f12dim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::get -- dimension of A does not match",__FILE__,__LINE__);
    unsigned int gij = ij;
    for (unsigned int g=0; g<nf12_; g++, gij+=nij)
      A.set_element(gij, Aij.get_element(g));
  }
  
  void MP2R12EnergyUtil_Nondiag::invert(RefSymmSCMatrix& A) const {
    A->gen_invert_this();
  }
  
  void MP2R12EnergyUtil_Diag::invert(RefSymmSCMatrix& A) const {
    check_dims(A);
    RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    const int noo = oodim_.n();
    for (int ij=0; ij<noo; ij++) {
      get(ij, A, Aij);
      Aij->gen_invert_this();
      put(ij, A, Aij);
    }
  }
  
  RefDiagSCMatrix MP2R12EnergyUtil_Nondiag::eigenvalues(
                                                                  const RefSymmSCMatrix& A) const {
    return A.eigvals();
  }
  
  RefDiagSCMatrix MP2R12EnergyUtil_Diag::eigenvalues(
                                                                 const RefSymmSCMatrix& A) const {
    check_dims(A);
    RefDiagSCMatrix evals = A.kit()->diagmatrix(f12dim_);
    RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    const int noo = oodim_.n();
    for (int ij=0; ij<noo; ij++) {
      get(ij, A, Aij);
      RefDiagSCMatrix evals_ij = Aij.eigvals();
      put(ij, evals, evals_ij);
    }
    
    return evals;
  }
  
  void MP2R12EnergyUtil_Nondiag::diagonalize(
                                                       const RefSymmSCMatrix& A,
                                                       RefDiagSCMatrix& evals,
                                                       RefSCMatrix& evecs) const {
    evals = A.kit()->diagmatrix(A.dim());
    evecs = A.kit()->matrix(A.dim(), A.dim());
    A.diagonalize(evals, evecs);
  }
  
  void MP2R12EnergyUtil_Diag::diagonalize(const RefSymmSCMatrix& A,
                                                      RefDiagSCMatrix& evals,
                                                      RefSCMatrix& evecs) const {
    check_dims(A);
    evals = A.kit()->diagmatrix(A.dim());
    evecs = A.kit()->matrix(A.dim(), A.dim());
    
    RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    RefDiagSCMatrix evals_ij = A.kit()->diagmatrix(gdim_);
    RefSCMatrix evecs_ij = A.kit()->matrix(gdim_, gdim_);
    const unsigned int nij = oodim_.n();
    for (unsigned int ij=0; ij<nij; ++ij) {
      get(ij, A, Aij);
      Aij.diagonalize(evals_ij, evecs_ij);
      put(ij, evals, evals_ij);
      put(ij, evecs, evecs_ij);
    }
  }
  
  void MP2R12EnergyUtil_Nondiag::transform(const RefSymmSCMatrix& B,
                                                     const RefDiagSCMatrix& A,
                                                     const RefSCMatrix& U) const {
    B.assign(0.0);
    B.accumulate_transform(U, A);
  }
  
  void MP2R12EnergyUtil_Diag::transform(const RefSymmSCMatrix& B,
                                                    const RefDiagSCMatrix& A,
                                                    const RefSCMatrix& U) const {
    check_dims(B);
    check_dims(U);
    B.assign(0.0);
    RefSymmSCMatrix Bij = B.kit()->symmmatrix(gdim_);
    RefDiagSCMatrix Aij = A.kit()->diagmatrix(gdim_);
    RefSCMatrix Uij = U.kit()->matrix(gdim_, gdim_);
    const unsigned int nij = oodim_.n();
    for (unsigned int ij=0; ij<nij; ++ij) {
      get(ij, A, Aij);
      get(ij, U, Uij);
      Bij.assign(0.0);
      Bij.accumulate_transform(Uij, Aij);
      put(ij, B, Bij);
    }
  }
  
  // Solves A*X = B
  void MP2R12EnergyUtil_Nondiag::solve_linear_system(
                                                               const RefSymmSCMatrix& A,
                                                               RefSCMatrix& X,
                                                               const RefSCMatrix& B) const {
    sc::exp::lapack_linsolv_symmnondef(A, X, B);
  }
  
  void MP2R12EnergyUtil_Diag::solve_linear_system(
                                                              const RefSymmSCMatrix& A,
                                                              RefSCMatrix& X,
                                                              const RefSCMatrix& B) const {
    check_dims(A);
    check_dims(X);
    check_dims(B);
    RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    RefSCVector Bij = A.kit()->vector(gdim_);
    RefSCVector Xij = A.kit()->vector(gdim_);
    const int noo = oodim_.n();
    for (int ij=0; ij<noo; ij++) {
      get(ij, A, Aij);
      get(ij, B, Bij);
      sc::exp::lapack_linsolv_symmnondef(Aij, Xij, Bij);
      put(ij, X, Xij);
    }
  }
  
  // computes y = A x
  void MP2R12EnergyUtil_Nondiag::times(const RefSymmSCMatrix& A,
                                                 const RefSCMatrix& x,
                                                 RefSCMatrix& y) const {
    y = A*x;
  }
  
  void MP2R12EnergyUtil_Diag::times(const RefSymmSCMatrix& A,
                                                const RefSCMatrix& x,
                                                RefSCMatrix& y) const {
    check_dims(A);
    check_dims(x);
    check_dims(y);
    RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    RefSCVector xij = A.kit()->vector(gdim_);
    RefSCVector yij = A.kit()->vector(gdim_);
    const int noo = oodim_.n();
    for (int ij=0; ij<noo; ij++) {
      get(ij, A, Aij);
      get(ij, x, xij);
      yij = Aij * xij;
      put(ij, y, yij);
    }
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
  
  RefSCVector MP2R12EnergyUtil_Diag::dot_product(
                                                             const RefSCMatrix& A,
                                                             const RefSCMatrix& B) const {
    check_dims(A);
    check_dims(B);
    if (A.coldim().n() != B.coldim().n() ||A.rowdim().n() != B.rowdim().n() ||A.coldim().n() != oodim_.n() ||A.rowdim().n() != f12dim_.n())
      throw ProgrammingError("MP2R12EnergyUtil::dot_product -- dimentions do not match",__FILE__,__LINE__);
    const int noo = oodim_.n();
    RefSCVector result = A.kit()->vector(oodim_);
    RefSCVector Aij = A.kit()->vector(gdim_);
    RefSCVector Bij = A.kit()->vector(gdim_);
    for (int ij=0; ij<noo; ij++) {
      get(ij, A, Aij);
      get(ij, B, Bij);
      result(ij) = Aij.dot(Bij);
    }
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
  
  void MP2R12EnergyUtil_Diag::print(const char* label,
                                                const RefSCMatrix& A,
                                                std::ostream& os) const {
    os << indent << label << ":"<< endl;
    RefSCVector Aij = A.kit()->vector(gdim_);
    const unsigned int noo = oodim_.n();
    for (unsigned int ij=0; ij<noo; ++ij) {
      get(ij, A, Aij);
      ostringstream oss;
      oss << "Block "<< ij;
      Aij.print(oss.str().c_str(), os);
    }
  }
  
  void MP2R12EnergyUtil_Diag::print(const char* label,
                                                const RefSymmSCMatrix& A,
                                                std::ostream& os) const {
    os << indent << label << ":"<< endl;
    RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
    const unsigned int noo = oodim_.n();
    for (unsigned int ij=0; ij<noo; ++ij) {
      get(ij, A, Aij);
      ostringstream oss;
      oss << "Block "<< ij;
      Aij.print(oss.str().c_str(), os);
    }
  }
  
  void MP2R12EnergyUtil_Diag::print(const char* label,
                                                const RefDiagSCMatrix& A,
                                                std::ostream& os) const {
    os << indent << label << ":"<< endl;
    RefDiagSCMatrix Aij = A.kit()->diagmatrix(gdim_);
    const unsigned int noo = oodim_.n();
    for (unsigned int ij=0; ij<noo; ++ij) {
      get(ij, A, Aij);
      ostringstream oss;
      oss << "Block "<< ij;
      Aij.print(oss.str().c_str(), os);
    }
  }

}

