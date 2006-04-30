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
#include <chemistry/qc/mbptr12/svd.h>

namespace sc {
  
  /** Class MP2R12EnergyUtil_base is the abstract interface to some utility functions.
  */
  class MP2R12EnergyUtil_base : virtual public RefCount {
    public:
    MP2R12EnergyUtil_base();
    virtual ~MP2R12EnergyUtil_base();
    
    /// Prints A
    virtual void print(const char* label, const RefSCMatrix& A) const = 0;
    /// Prints A
    virtual void print(const char* label, const RefSymmSCMatrix& A) const = 0;
    /// Inverts A in-place
    virtual void invert(RefSymmSCMatrix& A) const = 0;
    /// Solves A*X = B
    virtual void solve_linear_system(const RefSymmSCMatrix& A,
                                     RefSCMatrix& X,
                                     const RefSCMatrix& B) const = 0;
    /// computes y = A x
    virtual void times(const RefSymmSCMatrix& A,
                       const RefSCMatrix& x,
                       RefSCMatrix& y) const = 0;
    /// Computes "dot" product of A and B: tr[i] = sum_j A[j][i] B[j][i]
    virtual RefSCVector dot_product(const RefSCMatrix& A,
                                    const RefSCMatrix& B) const = 0;
  };
  
  /** Class MP2R12EnergyUtil provides some misc functions to operate on (blocked)
      occ-occ matrices. If Diag is true, then each occ-occ block has nonzero diagonal.
      f12-f12 matrix then has nf12 by nf12 blocks, each with diagonal structure.
  */
  template <bool Diag>
    class MP2R12EnergyUtil : public MP2R12EnergyUtil_base {
    public:
    /// oodim is a dimension of active occupied pairs. f12dim is has rank nf12*noo.
    MP2R12EnergyUtil(const RefSCDimension& oodim,
                     const RefSCDimension& f12dim);
    ~MP2R12EnergyUtil() {}
    
    /// Prints A
    void print(const char* label, const RefSCMatrix& A) const;
    /// Prints A
    void print(const char* label, const RefSymmSCMatrix& A) const;
    /// Inverts A in-place
    void invert(RefSymmSCMatrix& A) const;
    /// Solves A*X = B
    void solve_linear_system(const RefSymmSCMatrix& A,
                             RefSCMatrix& X,
                             const RefSCMatrix& B) const;
    /// computes y = A x
    void times(const RefSymmSCMatrix& A,
               const RefSCMatrix& x,
               RefSCMatrix& y) const;
    /// Computes "dot" product of A and B: tr[i] = sum_j A[j][i] B[j][i]
    RefSCVector dot_product(const RefSCMatrix& A,
                            const RefSCMatrix& B) const;
    
    private:
    // Number of pairs
    RefSCDimension oodim_;
    // number of geminals per pair times number of pairs
    RefSCDimension f12dim_;
    // number of geminals
    RefSCDimension gdim_;
    // number of geminals
    unsigned int nf12_;
    
    /// Checks if matrix A has proper dimensions. Throw, if not.
    void check_dims(const RefSCMatrix& A) const;
    /// Checks if matrix A has proper dimensions. Throw, if not.
    void check_dims(const RefSymmSCMatrix& A) const;
    /// Number of oo blocks in row dimension of A
    unsigned int nrowblks(const RefSCMatrix& A) const { check_dims(A); return A.rowdim().n()/oodim_.n(); }
    /// Number of oo blocks in column dimension of A
    unsigned int ncolblks(const RefSCMatrix& A) const { check_dims(A); return A.coldim().n()/oodim_.n(); }
    /// Number of oo blocks in dimension of A
    unsigned int nblks(const RefSymmSCMatrix& A) const { check_dims(A); return A.dim().n()/oodim_.n(); }
    
    /// gets ij block of A
    void get(unsigned int ij, const RefSCMatrix& A, RefSCVector& Aij) const;
    /// gets ij block of A
    void get(unsigned int ij, const RefSymmSCMatrix& A, RefSymmSCMatrix& Aij) const;
    /// puts ij block into A
    void put(unsigned int ij, RefSCMatrix& A, const RefSCVector& Aij) const;
    /// puts ij block into A
    void put(unsigned int ij, RefSymmSCMatrix& A, const RefSymmSCMatrix& Aij) const;
  };
  
  template<bool Diag>
    MP2R12EnergyUtil<Diag>::MP2R12EnergyUtil(const RefSCDimension& oodim,
                                       const RefSCDimension& f12dim) :
       oodim_(oodim), f12dim_(f12dim), nf12_(f12dim.n()/oodim.n())
    {
      gdim_ = new SCDimension(nf12_);
      if (f12dim_.n()%oodim_.n())
        throw ProgrammingError("MP2R12EnergyUtils::MP2R12EnergyUtils -- rank of f12dim must be divisible by rank of oodim",__FILE__,__LINE__);
    }
  
  template <bool Diag>
    void MP2R12EnergyUtil<Diag>::check_dims(const RefSCMatrix& A) const
    {
      const int nrow = A.rowdim().n();
      const int ncol = A.coldim().n();
      if (nrow != f12dim_.n() && nrow != oodim_.n())
        throw ProgrammingError("MP2R12EnergyUtil::check_dims -- row dimension does not match",__FILE__,__LINE__);
      if (ncol != f12dim_.n() && ncol != oodim_.n())
        throw ProgrammingError("MP2R12EnergyUtil::check_dims -- column dimension does not match",__FILE__,__LINE__);
    }
  template <bool Diag>
    void MP2R12EnergyUtil<Diag>::check_dims(const RefSymmSCMatrix& A) const
    {
      const int n = A.dim().n();
      if (n != f12dim_.n() && n != oodim_.n())
        throw ProgrammingError("MP2R12EnergyUtil::check_dims -- dimension does not match",__FILE__,__LINE__);
    }

  // Prints out matrix A
  template <bool Diag>
    void MP2R12EnergyUtil<Diag>::print(const char* label, const RefSCMatrix& A) const {
      A.print(label);
    }
  template <bool Diag>
    void MP2R12EnergyUtil<Diag>::print(const char* label, const RefSymmSCMatrix& A) const {
      A.print(label);
    }
  
}

#endif // include guard
