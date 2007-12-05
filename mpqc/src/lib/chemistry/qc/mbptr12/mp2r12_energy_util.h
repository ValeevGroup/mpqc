//
// mp2r12_energy_util.h
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
#pragma interface
#endif

#ifndef _chemistry_qc_mbptr12_mp2r12energyutil_h
#define _chemistry_qc_mbptr12_mp2r12energyutil_h

#include <util/ref/ref.h>
#include <math/scmat/matrix.h>
#include <util/class/scexception.h>
#include <chemistry/qc/mbptr12/spin.h>
#include <chemistry/qc/mbptr12/svd.h>

namespace sc {
  
  /** Class MP2R12EnergyUtil_base is the abstract interface to utility functions used by MP2R12Energy derivatives.
  */
  class MP2R12EnergyUtil_base : virtual public RefCount {
    protected:
      /// number of active occupied orbitals
      unsigned int nocc_act_;
      /// number of ij pairs
      RefSCDimension oodim_;
      /// number of xy pairs
      RefSCDimension xydim_;
      /// number of geminals per pair times number of xy pairs
      RefSCDimension f12dim_;
      /// number of geminals
      RefSCDimension gdim_;
      /// number of geminals
      unsigned int nf12_;

    public:
    MP2R12EnergyUtil_base();
    MP2R12EnergyUtil_base(const RefSCDimension& oodim,
                          const RefSCDimension& xydim,
                          const RefSCDimension& f12dim,
                          const unsigned int nocc_act);
    virtual ~MP2R12EnergyUtil_base();
    
    /// Checks if matrix A has proper dimensions. Throw, if not.
    void check_dims(const RefSCMatrix& A) const;
    /// Checks if matrix A has proper dimensions. Throw, if not.
    void check_dims(const RefSymmSCMatrix& A) const;
    
    /// Prints A
    virtual void print(const char* label, const RefSCMatrix& A, std::ostream& os = ExEnv::out0()) const = 0;
    /// Prints A
    virtual void print(const char* label, const RefSymmSCMatrix& A, std::ostream& os = ExEnv::out0()) const = 0;
    /// Prints A
    virtual void print(const char* label, const RefDiagSCMatrix& A, std::ostream& os = ExEnv::out0()) const = 0;
    /// Inverts A in-place
    virtual void invert(RefSymmSCMatrix& A) const = 0;
    /// Computes eigenvalues of A
    virtual RefDiagSCMatrix eigenvalues(const RefSymmSCMatrix& A) const = 0;
    /// Computes eigenvalues and eigenvectors of A. evals and evecs don't have to be allocated
    virtual void diagonalize(const RefSymmSCMatrix& A, RefDiagSCMatrix& evals, RefSCMatrix& evecs) const = 0;
    /// B = U * A * U.t()
    virtual void transform(const RefSymmSCMatrix& B, const RefDiagSCMatrix& A, const RefSCMatrix& U) const = 0;
    /// Solves A*X = B
    virtual void solve_linear_system(const RefSymmSCMatrix& A,
                                     RefSCMatrix& X,
                                     const RefSCMatrix& B) const = 0;
    virtual void solve_linear_system(unsigned int ij,
                                     const RefSymmSCMatrix& A,
                                     RefSCMatrix& X,
                                     const RefSCMatrix& B) const = 0;
    /// computes y = A x
    virtual void times(const RefSymmSCMatrix& A,
                       const RefSCMatrix& x,
                       RefSCMatrix& y) const = 0;
    /// computes y = A x
    virtual void times(unsigned int ij,
                       const RefSymmSCMatrix& A,
                       const RefSCMatrix& x,
                       RefSCMatrix& y) const = 0;
    /// Computes "dot" product of A and B: tr[i] = sum_j A[j][i] B[j][i]
    virtual RefSCVector dot_product(const RefSCMatrix& A,
                                    const RefSCMatrix& B) const = 0;
  };
  
  class MP2R12EnergyUtil_Diag : virtual public RefCount {
    protected:
      /// number of active occupied orbitals
      unsigned int nocc_act_;
      /// number of ij pairs
      RefSCDimension oodim_;
      /// number of xy pairs
      RefSCDimension xydim_;
      /// number of geminals per pair times number of xy pairs
      RefSCDimension f12dim_;
      /// number of geminals
      RefSCDimension gdim_;
      /// number of geminals
      unsigned int nf12_;

    public:
    MP2R12EnergyUtil_Diag();
    MP2R12EnergyUtil_Diag(const RefSCDimension& oodim,
                          const RefSCDimension& xydim,
                          const RefSCDimension& f12dim,
                          const unsigned int nocc_act);
    virtual ~MP2R12EnergyUtil_Diag();
    
    /// Checks if matrix A has proper dimensions. Throw, if not.
    void check_dims(const RefSCMatrix& A) const;
    /// Checks if matrix A has proper dimensions. Throw, if not.
    void check_dims(const RefSymmSCMatrix& A) const;
    
    /// Prints A
    virtual void print(const char* label, const RefSCMatrix& A, std::ostream& os = ExEnv::out0()) const = 0;
    /// Prints A
    virtual void print(const char* label, const RefSymmSCMatrix& A, std::ostream& os = ExEnv::out0()) const = 0;
    /// Prints A
    virtual void print(const char* label, const RefDiagSCMatrix& A, std::ostream& os = ExEnv::out0()) const = 0;
    
    /// gets ij block of A
    virtual void get(unsigned int ij, const RefSCMatrix& A, const RefSCVector& Aij) const = 0;
    /// gets ij block of A
    virtual void get(unsigned int ij, const RefSCMatrix& A, const RefSCMatrix& Aij) const =0;
    /// gets ij block of A
    virtual void get(unsigned int ij, const RefSymmSCMatrix& A, const RefSymmSCMatrix& Aij) const = 0;
    /// gets ij block of A
    virtual void get(unsigned int ij, const RefDiagSCMatrix& A, const RefDiagSCMatrix& Aij) const =0;
    
    /// puts ij block into A
    virtual void put(unsigned int ij, const RefSCMatrix& A, const RefSCVector& Aij) const =0;
    /// puts ij block into A
    virtual void put(unsigned int ij, const RefSCMatrix& A, const RefSCMatrix& Aij) const =0;
    /// puts ij block into A
    virtual void put(unsigned int ij, const RefSymmSCMatrix& A, const RefSymmSCMatrix& Aij) const = 0;
    /// puts ij block into A
    virtual void put(unsigned int ij, const RefDiagSCMatrix& A, const RefDiagSCMatrix& Aij) const = 0;
    
    /// Inverts A in-place
    virtual void invert(RefSymmSCMatrix& A) const = 0;
    /// Computes eigenvalues of A
    virtual RefDiagSCMatrix eigenvalues(const RefSymmSCMatrix& A) const = 0;
    /// Computes eigenvalues and eigenvectors of A. evals and evecs don't have to be allocated
    virtual void diagonalize(const RefSymmSCMatrix& A, RefDiagSCMatrix& evals, RefSCMatrix& evecs) const = 0;
    /// B = U * A * U.t()
    virtual void transform(const RefSymmSCMatrix& B, const RefDiagSCMatrix& A, const RefSCMatrix& U) const = 0;
    /// Solves A*X = B
    virtual void solve_linear_system(const RefSymmSCMatrix& A,
                                     RefSCMatrix& X,
                                     const RefSCMatrix& B) const = 0;
    virtual void solve_linear_system(unsigned int ij,
                                     const RefSymmSCMatrix& A,
                                     RefSCMatrix& X,
                                     const RefSCMatrix& B) const = 0;
    /// computes y = A x
    virtual void times(const RefSymmSCMatrix& A,
                       const RefSCMatrix& x,
                       RefSCMatrix& y) const = 0;
    /// computes y = A x
    virtual void times(unsigned int ij,
                       const RefSymmSCMatrix& A,
                       const RefSCMatrix& x,
                       RefSCMatrix& y) const = 0;
    /// Computes "dot" product of A and B: tr[i] = sum_j A[j][i] B[j][i]
    virtual RefSCVector dot_product(const RefSCMatrix& A,
                                    const RefSCMatrix& B) const = 0;
  };
  
  /** Class MP2R12EnergyUtil provides some misc functions to operate on (blocked)
      ijxy and xyxy matrices. If Diag is true, then each xy block is the same size as ij, and ijxy and xyxy blocks are "diagonal".
      f12-f12 matrix then has nf12 by nf12 blocks, each with diagonal structure.
  */
  
  class MP2R12EnergyUtil_Diag_DifferentSpin : public MP2R12EnergyUtil_base {
    private:
      SpinCase2 spincase2_;
    public:
      /// oodim is a dimension of ij (active occupied pairs)
      /// xydim is a dimension of xy (orbital products used to generate geminal space)
      /// f12dim is has rank nf12*nxy.
      /// nocc_act is the number of active occupied orbitals.
      MP2R12EnergyUtil_Diag_DifferentSpin(const RefSCDimension& oodim,
                                           const RefSCDimension& xydim,
                                           const RefSCDimension& f12dim,
                                           const unsigned int nocc_act);
      ~MP2R12EnergyUtil_Diag_DifferentSpin() {}
      
      /// Number of oo blocks in row dimension of A
      unsigned int nrowblks(const RefSCMatrix& A) const;
      /// Number of oo blocks in column dimension of A
      unsigned int ncolblks(const RefSCMatrix& A) const;
      /// Number of oo blocks in dimension of A
      unsigned int nblks(const RefSymmSCMatrix& A) const;
      
      /// gets ij block of A
      void get(unsigned int ij, const RefSCMatrix& A, const RefSCVector& Aij) const;
      /// gets ij block of A
      void get(unsigned int ij, const RefSCMatrix& A, const RefSCMatrix& Aij) const;
      /// gets ij block of A
      void get(unsigned int ij, const RefSymmSCMatrix& A, const RefSymmSCMatrix& Aij) const;
      /// gets ij block of A
      void get(unsigned int ij, const RefDiagSCMatrix& A, const RefDiagSCMatrix& Aij) const;
      
      /// puts ij block into A
      void put(unsigned int ij, const RefSCMatrix& A, const RefSCVector& Aij) const;
      /// puts ij block into A
      void put(unsigned int ij, const RefSCMatrix& A, const RefSCMatrix& Aij) const;
      /// puts ij block into A
      void put(unsigned int ij, const RefSymmSCMatrix& A, const RefSymmSCMatrix& Aij) const;
      /// puts ij block into A
      void put(unsigned int ij, const RefDiagSCMatrix& A, const RefDiagSCMatrix& Aij) const;
      
      /// Inverts A in-place
      void invert(RefSymmSCMatrix& A) const;
      /// Computes eigenvalues of A
      RefDiagSCMatrix eigenvalues(const RefSymmSCMatrix& A) const;
      /// Computes eigenvalues and eigenvectors of A. evals and evecs don't have to be allocated
      void diagonalize(const RefSymmSCMatrix& A, RefDiagSCMatrix& evals, RefSCMatrix& evecs) const;
      /// B = U * A * U.t()
      void transform(const RefSymmSCMatrix& B, const RefDiagSCMatrix& A, const RefSCMatrix& U) const;
      /// Solves A*X = B
      void solve_linear_system(const RefSymmSCMatrix& A,
                               RefSCMatrix& X,
                               const RefSCMatrix& B) const;
      void solve_linear_system(unsigned int ij,
                               const RefSymmSCMatrix& A,
                               RefSCMatrix& X,
                               const RefSCMatrix& B) const;
      /// computes y = A x
      void times(const RefSymmSCMatrix& A,
                 const RefSCMatrix& x,
                 RefSCMatrix& y) const;
      /// computes y = A x
      void times(unsigned int ij,
                 const RefSymmSCMatrix& A,
                 const RefSCMatrix& x,
                 RefSCMatrix& y) const;
      /// Computes "dot" product of A and B: tr[i] = sum_j A[j][i] B[j][i]
      RefSCVector dot_product(const RefSCMatrix& A,
                              const RefSCMatrix& B) const;
      
      /// Prints A
      void print(const char* label, const RefSCMatrix& A, std::ostream& os = ExEnv::out0()) const;
      /// Prints A
      void print(const char* label, const RefSymmSCMatrix& A, std::ostream& os = ExEnv::out0()) const;
      /// Prints A
      void print(const char* label, const RefDiagSCMatrix& A, std::ostream& os = ExEnv::out0()) const;
  };
  
  class MP2R12EnergyUtil_Diag_SameSpin : public MP2R12EnergyUtil_base {
    private:
      SpinCase2 spincase2_;
    public:
      /// oodim is a dimension of ij (active occupied pairs)
      /// xydim is a dimension of xy (orbital products used to generate geminal space)
      /// f12dim is has rank nf12*nxy.
      /// nocc_act is the number of active occupied orbitals.
      MP2R12EnergyUtil_Diag_SameSpin(const RefSCDimension& oodim,
                                      const RefSCDimension& xydim,
                                      const RefSCDimension& f12dim,
                                      const unsigned int nocc_act);
      ~MP2R12EnergyUtil_Diag_SameSpin() {}
      
      /// Number of oo blocks in row dimension of A
      unsigned int nrowblks(const RefSCMatrix& A) const;
      /// Number of oo blocks in column dimension of A
      unsigned int ncolblks(const RefSCMatrix& A) const;
      /// Number of oo blocks in dimension of A
      unsigned int nblks(const RefSymmSCMatrix& A) const;
      
      /// gets ij block of A
      void get(unsigned int ij, const RefSCMatrix& A, const RefSCVector& Aij) const;
      /// gets ij block of A
      void get(unsigned int ij, const RefSCMatrix& A, const RefSCMatrix& Aij) const;
      /// gets ij block of A
      void get(unsigned int ij, const RefSymmSCMatrix& A, const RefSymmSCMatrix& Aij) const;
      /// gets ij block of A
      void get(unsigned int ij, const RefDiagSCMatrix& A, const RefDiagSCMatrix& Aij) const;
      
      /// puts ij block into A
      void put(unsigned int ij, const RefSCMatrix& A, const RefSCVector& Aij) const;
      /// puts ij block into A
      void put(unsigned int ij, const RefSCMatrix& A, const RefSCMatrix& Aij) const;
      /// puts ij block into A
      void put(unsigned int ij, const RefSymmSCMatrix& A, const RefSymmSCMatrix& Aij) const;
      /// puts ij block into A
      void put(unsigned int ij, const RefDiagSCMatrix& A, const RefDiagSCMatrix& Aij) const;
      
      /// Inverts A in-place
      void invert(RefSymmSCMatrix& A) const;
      /// Computes eigenvalues of A
      RefDiagSCMatrix eigenvalues(const RefSymmSCMatrix& A) const;
      /// Computes eigenvalues and eigenvectors of A. evals and evecs don't have to be allocated
      void diagonalize(const RefSymmSCMatrix& A, RefDiagSCMatrix& evals, RefSCMatrix& evecs) const;
      /// B = U * A * U.t()
      void transform(const RefSymmSCMatrix& B, const RefDiagSCMatrix& A, const RefSCMatrix& U) const;
      /// Solves A*X = B
      void solve_linear_system(const RefSymmSCMatrix& A,
                               RefSCMatrix& X,
                               const RefSCMatrix& B) const;
      void solve_linear_system(unsigned int ij,
                               const RefSymmSCMatrix& A,
                               RefSCMatrix& X,
                               const RefSCMatrix& B) const;
      /// computes y = A x
      void times(const RefSymmSCMatrix& A,
                 const RefSCMatrix& x,
                 RefSCMatrix& y) const;
      void times(unsigned int ij,
                 const RefSymmSCMatrix& A,
                 const RefSCMatrix& x,
                 RefSCMatrix& y) const;
      /// Computes "dot" product of A and B: tr[i] = sum_j A[j][i] B[j][i]
      RefSCVector dot_product(const RefSCMatrix& A,
                              const RefSCMatrix& B) const;
      
      /// Prints A
      void print(const char* label, const RefSCMatrix& A, std::ostream& os = ExEnv::out0()) const;
      /// Prints A
      void print(const char* label, const RefSymmSCMatrix& A, std::ostream& os = ExEnv::out0()) const;
      /// Prints A
      void print(const char* label, const RefDiagSCMatrix& A, std::ostream& os = ExEnv::out0()) const;
  };
  
  class MP2R12EnergyUtil_Diag_DifferentSpin_indep : public MP2R12EnergyUtil_Diag {
    private:
      SpinCase2 spincase2_;
    public:
      /// oodim is a dimension of ij (active occupied pairs)
      /// xydim is a dimension of xy (orbital products used to generate geminal space)
      /// f12dim is has rank nf12*nxy.
      /// nocc_act is the number of active occupied orbitals.
      MP2R12EnergyUtil_Diag_DifferentSpin_indep(const RefSCDimension& oodim,
                                           const RefSCDimension& xydim,
                                           const RefSCDimension& f12dim,
                                           const unsigned int nocc_act);
      ~MP2R12EnergyUtil_Diag_DifferentSpin_indep() {}
      
      /// Number of oo blocks in row dimension of A
      unsigned int nrowblks(const RefSCMatrix& A) const;
      /// Number of oo blocks in column dimension of A
      unsigned int ncolblks(const RefSCMatrix& A) const;
      /// Number of oo blocks in dimension of A
      unsigned int nblks(const RefSymmSCMatrix& A) const;
      
      /// gets ij block of A
      void get(unsigned int ij, const RefSCMatrix& A, const RefSCVector& Aij) const;
      /// gets ij block of A
      void get(unsigned int ij, const RefSCMatrix& A, const RefSCMatrix& Aij) const;
      /// gets ij block of A
      void get(unsigned int ij, const RefSymmSCMatrix& A, const RefSymmSCMatrix& Aij) const;
      /// gets ij block of A
      void get(unsigned int ij, const RefDiagSCMatrix& A, const RefDiagSCMatrix& Aij) const;
      
      /// puts ij block into A
      void put(unsigned int ij, const RefSCMatrix& A, const RefSCVector& Aij) const;
      /// puts ij block into A
      void put(unsigned int ij, const RefSCMatrix& A, const RefSCMatrix& Aij) const;
      /// puts ij block into A
      void put(unsigned int ij, const RefSymmSCMatrix& A, const RefSymmSCMatrix& Aij) const;
      /// puts ij block into A
      void put(unsigned int ij, const RefDiagSCMatrix& A, const RefDiagSCMatrix& Aij) const;
      
      /// Inverts A in-place
      void invert(RefSymmSCMatrix& A) const;
      /// Computes eigenvalues of A
      RefDiagSCMatrix eigenvalues(const RefSymmSCMatrix& A) const;
      /// Computes eigenvalues and eigenvectors of A. evals and evecs don't have to be allocated
      void diagonalize(const RefSymmSCMatrix& A, RefDiagSCMatrix& evals, RefSCMatrix& evecs) const;
      /// B = U * A * U.t()
      void transform(const RefSymmSCMatrix& B, const RefDiagSCMatrix& A, const RefSCMatrix& U) const;
      /// Solves A*X = B
      void solve_linear_system(const RefSymmSCMatrix& A,
                               RefSCMatrix& X,
                               const RefSCMatrix& B) const;
      void solve_linear_system(unsigned int ij,
                               const RefSymmSCMatrix& A,
                               RefSCMatrix& X,
                               const RefSCMatrix& B) const;
      /// computes y = A x
      void times(const RefSymmSCMatrix& A,
                 const RefSCMatrix& x,
                 RefSCMatrix& y) const;
      /// computes y = A x
      void times(unsigned int ij,
                 const RefSymmSCMatrix& A,
                 const RefSCMatrix& x,
                 RefSCMatrix& y) const;
      /// Computes "dot" product of A and B: tr[i] = sum_j A[j][i] B[j][i]
      RefSCVector dot_product(const RefSCMatrix& A,
                              const RefSCMatrix& B) const;
      
      /// Prints A
      void print(const char* label, const RefSCMatrix& A, std::ostream& os = ExEnv::out0()) const;
      /// Prints A
      void print(const char* label, const RefSymmSCMatrix& A, std::ostream& os = ExEnv::out0()) const;
      /// Prints A
      void print(const char* label, const RefDiagSCMatrix& A, std::ostream& os = ExEnv::out0()) const;
  };
  
  class MP2R12EnergyUtil_Diag_SameSpin_indep : public MP2R12EnergyUtil_Diag {
    private:
      SpinCase2 spincase2_;
    public:
      /// oodim is a dimension of ij (active occupied pairs)
      /// xydim is a dimension of xy (orbital products used to generate geminal space)
      /// f12dim is has rank nf12*nxy.
      /// nocc_act is the number of active occupied orbitals.
      MP2R12EnergyUtil_Diag_SameSpin_indep(const RefSCDimension& oodim,
                                      const RefSCDimension& xydim,
                                      const RefSCDimension& f12dim,
                                      const unsigned int nocc_act);
      ~MP2R12EnergyUtil_Diag_SameSpin_indep() {}
      
      /// Number of oo blocks in row dimension of A
      unsigned int nrowblks(const RefSCMatrix& A) const;
      /// Number of oo blocks in column dimension of A
      unsigned int ncolblks(const RefSCMatrix& A) const;
      /// Number of oo blocks in dimension of A
      unsigned int nblks(const RefSymmSCMatrix& A) const;
      
      /// gets ij block of A
      void get(unsigned int ij, const RefSCMatrix& A, const RefSCVector& Aij) const;
      /// gets ij block of A
      void get(unsigned int ij, const RefSCMatrix& A, const RefSCMatrix& Aij) const;
      /// gets ij block of A
      void get(unsigned int ij, const RefSymmSCMatrix& A, const RefSymmSCMatrix& Aij) const;
      /// gets ij block of A
      void get(unsigned int ij, const RefDiagSCMatrix& A, const RefDiagSCMatrix& Aij) const;
      
      /// puts ij block into A
      void put(unsigned int ij, const RefSCMatrix& A, const RefSCVector& Aij) const;
      /// puts ij block into A
      void put(unsigned int ij, const RefSCMatrix& A, const RefSCMatrix& Aij) const;
      /// puts ij block into A
      void put(unsigned int ij, const RefSymmSCMatrix& A, const RefSymmSCMatrix& Aij) const;
      /// puts ij block into A
      void put(unsigned int ij, const RefDiagSCMatrix& A, const RefDiagSCMatrix& Aij) const;
      
      /// Inverts A in-place
      void invert(RefSymmSCMatrix& A) const;
      /// Computes eigenvalues of A
      RefDiagSCMatrix eigenvalues(const RefSymmSCMatrix& A) const;
      /// Computes eigenvalues and eigenvectors of A. evals and evecs don't have to be allocated
      void diagonalize(const RefSymmSCMatrix& A, RefDiagSCMatrix& evals, RefSCMatrix& evecs) const;
      /// B = U * A * U.t()
      void transform(const RefSymmSCMatrix& B, const RefDiagSCMatrix& A, const RefSCMatrix& U) const;
      /// Solves A*X = B
      void solve_linear_system(const RefSymmSCMatrix& A,
                               RefSCMatrix& X,
                               const RefSCMatrix& B) const;
      void solve_linear_system(unsigned int ij,
                               const RefSymmSCMatrix& A,
                               RefSCMatrix& X,
                               const RefSCMatrix& B) const;
      /// computes y = A x
      void times(const RefSymmSCMatrix& A,
                 const RefSCMatrix& x,
                 RefSCMatrix& y) const;
      void times(unsigned int ij,
                 const RefSymmSCMatrix& A,
                 const RefSCMatrix& x,
                 RefSCMatrix& y) const;
      /// Computes "dot" product of A and B: tr[i] = sum_j A[j][i] B[j][i]
      RefSCVector dot_product(const RefSCMatrix& A,
                              const RefSCMatrix& B) const;
      
      /// Prints A
      void print(const char* label, const RefSCMatrix& A, std::ostream& os = ExEnv::out0()) const;
      /// Prints A
      void print(const char* label, const RefSymmSCMatrix& A, std::ostream& os = ExEnv::out0()) const;
      /// Prints A
      void print(const char* label, const RefDiagSCMatrix& A, std::ostream& os = ExEnv::out0()) const;
  };
  
    class MP2R12EnergyUtil_Nondiag : public MP2R12EnergyUtil_base {
    public:
    /// oodim is a dimension of ij (active occupied pairs)
    /// xydim is a dimension of xy (orbital products used to generate geminal space)
    /// f12dim is has rank nf12*nxy.
    MP2R12EnergyUtil_Nondiag(const RefSCDimension& oodim,
                             const RefSCDimension& xydim,
                             const RefSCDimension& f12dim,
                             const unsigned int nocc_act);
    ~MP2R12EnergyUtil_Nondiag() {}
    
    /// Prints A
    void print(const char* label, const RefSCMatrix& A, std::ostream& os = ExEnv::out0()) const;
    /// Prints A
    void print(const char* label, const RefSymmSCMatrix& A, std::ostream& os = ExEnv::out0()) const;
    /// Prints A
    void print(const char* label, const RefDiagSCMatrix& A, std::ostream& os = ExEnv::out0()) const;
    /// Inverts A in-place
    void invert(RefSymmSCMatrix& A) const;
    /// Computes eigenvalues of A
    RefDiagSCMatrix eigenvalues(const RefSymmSCMatrix& A) const;
    /// Computes eigenvalues and eigenvectors of A. evals and evecs don't have to be allocated
    void diagonalize(const RefSymmSCMatrix& A, RefDiagSCMatrix& evals, RefSCMatrix& evecs) const;
    /// B = U * A * U.t()
    void transform(const RefSymmSCMatrix& B, const RefDiagSCMatrix& A, const RefSCMatrix& U) const;
    /// Solves A*X = B
    void solve_linear_system(const RefSymmSCMatrix& A,
                             RefSCMatrix& X,
                             const RefSCMatrix& B) const;
    void solve_linear_system(unsigned int ij,
                             const RefSymmSCMatrix& A,
                             RefSCMatrix& X,
                             const RefSCMatrix& B) const;
    /// computes y = A x
    void times(const RefSymmSCMatrix& A,
               const RefSCMatrix& x,
               RefSCMatrix& y) const;
    void times(unsigned int ij,
               const RefSymmSCMatrix& A,
               const RefSCMatrix& x,
               RefSCMatrix& y) const;
    /// Computes "dot" product of A and B: tr[i] = sum_j A[j][i] B[j][i]
    RefSCVector dot_product(const RefSCMatrix& A,
                            const RefSCMatrix& B) const;

    /// Number of oo blocks in row dimension of A
    unsigned int nrowblks(const RefSCMatrix& A) const;
    /// Number of oo blocks in column dimension of A
    unsigned int ncolblks(const RefSCMatrix& A) const;
    /// Number of oo blocks in dimension of A
    unsigned int nblks(const RefSymmSCMatrix& A) const;
    
  };

  Ref<MP2R12EnergyUtil_Diag> generate_MP2R12EnergyUtil_Diag(SpinCase2 spincase2,
                                                            const RefSCDimension& oodim,
                                                            const RefSCDimension& xydim,
                                                            const RefSCDimension& f12dim,
                                                            const unsigned int nocc_act);
  
}

#endif // include guard
