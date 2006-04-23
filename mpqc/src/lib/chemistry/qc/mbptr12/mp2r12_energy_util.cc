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

MP2R12EnergyUtil_base::MP2R12EnergyUtil_base()
{}

MP2R12EnergyUtil_base::~MP2R12EnergyUtil_base()
{}

////

  // put/get can only be implemented when Diag=true
  template<>
    void MP2R12EnergyUtil<true>::get(unsigned int ij, const RefSCMatrix& A, RefSCVector& Aij) const
    {
      const unsigned int nij = oodim_.n();
      if (ij >= nij)
        throw ProgrammingError("MP2R12EnergyUtil::get -- ij >= nij",__FILE__,__LINE__);
      const int nrow = A.rowdim().n();
      const int ncol = A.coldim().n();
      if (nrow != f12dim_.n() || ncol != oodim_.n())
        throw ProgrammingError("MP2R12EnergyUtil::get -- incorrect dimensions of A",__FILE__,__LINE__);
      for(unsigned int g=0; g<nf12_; g++) {
        Aij.set_element(g,A.get_element(g*nij+ij,ij));
      }
    }
  template<>
    void MP2R12EnergyUtil<true>::get(unsigned int ij, const RefSymmSCMatrix& A, RefSymmSCMatrix& Aij) const
    {
      const unsigned int nij = oodim_.n();
      if (ij >= nij)
        throw ProgrammingError("MP2R12EnergyUtil::get -- ij >= nij",__FILE__,__LINE__);
      if (A.dim().n() != f12dim_.n())
        throw ProgrammingError("MP2R12EnergyUtil::get -- dimension of A does not match",__FILE__,__LINE__);
      for(unsigned int g=0; g<nf12_; g++)
        for(unsigned int f=0; f<=g; f++) {
          Aij.set_element(g,f,A.get_element(g*nij+ij,f*nij+ij));
        }
    }
  template<>
    void MP2R12EnergyUtil<true>::put(unsigned int ij, RefSCMatrix& A, const RefSCVector& Aij) const
    {
      const unsigned int nij = oodim_.n();
      if (ij >= nij)
        throw ProgrammingError("MP2R12EnergyUtil::get -- ij >= nij",__FILE__,__LINE__);
      const int nrow = A.rowdim().n();
      const int ncol = A.coldim().n();
      if (nrow != f12dim_.n() || ncol != oodim_.n())
        throw ProgrammingError("MP2R12EnergyUtil::get -- incorrect dimensions of A",__FILE__,__LINE__);
      for(unsigned int g=0; g<nf12_; g++) {
        A.set_element(g*nij+ij,ij,Aij.get_element(g));
      }
    }
  template<>
    void MP2R12EnergyUtil<true>::put(unsigned int ij, RefSymmSCMatrix& A, const RefSymmSCMatrix& Aij) const
    {
      const unsigned int nij = oodim_.n();
      if (ij >= nij)
        throw ProgrammingError("MP2R12EnergyUtil::get -- ij >= nij",__FILE__,__LINE__);
      if (A.dim().n() != f12dim_.n())
        throw ProgrammingError("MP2R12EnergyUtil::get -- dimension of A does not match",__FILE__,__LINE__);
      for(unsigned int g=0; g<nf12_; g++)
        for(unsigned int f=0; f<=g; f++) {
          A.set_element(g*nij+ij,f*nij+ij,Aij.get_element(g,f));
        }
    }

  template <>
    void MP2R12EnergyUtil<false>::invert(RefSymmSCMatrix& A) const {
      A->gen_invert_this();
    }
  template <>
    void MP2R12EnergyUtil<true>::invert(RefSymmSCMatrix& A) const {
      check_dims(A);
      RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
      const int noo = oodim_.n();
      for(int ij=0; ij<noo; ij++) {
        get(ij,A,Aij);
        Aij->gen_invert_this();
        put(ij,A,Aij);
      }
    }

  // Solves A*X = B
  template <>
    void MP2R12EnergyUtil<false>::solve_linear_system(const RefSymmSCMatrix& A, RefSCMatrix& X, const RefSCMatrix& B) const
    {
      sc::exp::lapack_linsolv_symmnondef(A, X, B);
    }
  template <>
    void MP2R12EnergyUtil<true>::solve_linear_system(const RefSymmSCMatrix& A, RefSCMatrix& X, const RefSCMatrix& B) const
    {
      check_dims(A);
      check_dims(X);
      check_dims(B);
      RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
      RefSCVector Bij = A.kit()->vector(gdim_);
      RefSCVector Xij = A.kit()->vector(gdim_);
      const int noo = oodim_.n();
      for(int ij=0; ij<noo; ij++) {
        get(ij,A,Aij);
        get(ij,B,Bij);
        sc::exp::lapack_linsolv_symmnondef(Aij, Xij, Bij);
        put(ij,X,Xij);
      }
    }
  
  // computes y = A x
  template <>
    void MP2R12EnergyUtil<false>::times(const RefSymmSCMatrix& A, const RefSCMatrix& x, RefSCMatrix& y) const
    {
      y = A*x;
    }
  template <>
    void MP2R12EnergyUtil<true>::times(const RefSymmSCMatrix& A, const RefSCMatrix& x, RefSCMatrix& y) const
    {
      check_dims(A);
      check_dims(x);
      check_dims(y);
      RefSymmSCMatrix Aij = A.kit()->symmmatrix(gdim_);
      RefSCVector xij = A.kit()->vector(gdim_);
      RefSCVector yij = A.kit()->vector(gdim_);
      const int noo = oodim_.n();
      for(int ij=0; ij<noo; ij++) {
        get(ij,A,Aij);
        get(ij,x,xij);
        yij = Aij * xij;
        put(ij,y,yij);
      }
    }
  
  template <>
    RefSCVector MP2R12EnergyUtil<false>::dot_product(const RefSCMatrix& A, const RefSCMatrix& B) const
    {
      check_dims(A);
      check_dims(B);
      if (A.coldim().n() != B.coldim().n() ||
          A.rowdim().n() != B.rowdim().n() ||
          A.coldim().n() != oodim_.n() ||
          A.rowdim().n() != f12dim_.n()
         )
         throw ProgrammingError("MP2R12EnergyUtil::dot_product -- dimentions do not match",__FILE__,__LINE__);
      RefSCMatrix AB = A.t() * B;
      const int noo = oodim_.n();
      RefSCVector result = AB.kit()->vector(oodim_);
      for(int ij=0; ij<noo; ij++)
        result(ij) = AB.get_element(ij,ij);
      return result;
    }
  template <>
    RefSCVector MP2R12EnergyUtil<true>::dot_product(const RefSCMatrix& A, const RefSCMatrix& B) const
    {
      check_dims(A);
      check_dims(B);
      if (A.coldim().n() != B.coldim().n() ||
          A.rowdim().n() != B.rowdim().n() ||
          A.coldim().n() != oodim_.n() ||
          A.rowdim().n() != f12dim_.n()
         )
        throw ProgrammingError("MP2R12EnergyUtil::dot_product -- dimentions do not match",__FILE__,__LINE__);
      const int noo = oodim_.n();
      RefSCVector result = A.kit()->vector(oodim_);
      RefSCVector Aij = A.kit()->vector(gdim_);
      RefSCVector Bij = A.kit()->vector(gdim_);
      for(int ij=0; ij<noo; ij++) {
        get(ij,A,Aij);
        get(ij,B,Bij);
        result(ij) = Aij.dot(Bij);
      }
      return result;
    }
