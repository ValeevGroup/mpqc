//
// int2e.h
//
// Copyright (C) 2001 Edward Valeev
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

#ifndef _chemistry_qc_libint2_int2e_h
#define _chemistry_qc_libint2_int2e_h

#include <limits.h>

#include <util/ref/ref.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/libint2/bounds.h>

namespace sc {

class Integral;
class Log2Bounds;

/** Int2eLibint2 is an interface to various specializations of two-electron
    integral evaluators implemented in Libint2. It is used by TwoBodyIntLibint2 and TwoBodyDerivIntLibint2 to
    implement IntegralLibint2. */
class Int2eLibint2: public RefCount {
  protected:
    Integral *integral_;

    Ref<GaussianBasisSet> bs1_;
    Ref<GaussianBasisSet> bs2_;
    Ref<GaussianBasisSet> bs3_;
    Ref<GaussianBasisSet> bs4_;

    Ref<Log2Bounds> bounds_;

    GaussianShell *int_shell1_;
    GaussianShell *int_shell2_;
    GaussianShell *int_shell3_;
    GaussianShell *int_shell4_;

    // Whether shells can be permuted
    int permute_;
    // Whether redundant integrals are needed or only unique ones
    int redundant_;

    /*--- Storage related stuff ---*/
    // Available storage
    size_t storage_;
    // Storage currently used
    size_t storage_used_;
    // Checks if too much storage is used
    void check_storage_() const;
    // Reports minimum "significant" storage needed to initialize the Int2e base
    static size_t storage_required_(const Ref<GaussianBasisSet>& b1,
				    const Ref<GaussianBasisSet>& b2 = 0,
				    const Ref<GaussianBasisSet>& b3 = 0,
				    const Ref<GaussianBasisSet>& b4 = 0);

    /*--- Scratch ---*/
    std::vector<double> tformbuf_;    // stores one partially transformed contraction quartet

    /*--- helper functions ---*/
    // cart.->sph.harm. transform functions
    void transform_contrquartets_(double *,double *);
    // normalize cartesian quartets, if needed
    void norm_contrcart1_(double* data);
    template <unsigned int ntypes> void norm_contrcart_(double* data);
    // sort from by-contraction-quartet to shell-quartet order
    void sort_contrquartets_to_shellquartet_(double *,double *);
    // permute perm_ints_ into target_int_buf_
    void permute_target_(double *, double *, int, int, int);
    void permute_1234_to_1243_(double *, double *);
    void permute_1234_to_2134_(double *, double *);
    void permute_1234_to_2143_(double *, double *);
    void permute_1234_to_3412_(double *, double *);
    void permute_1234_to_3421_(double *, double *);
    void permute_1234_to_4312_(double *, double *);
    void permute_1234_to_4321_(double *, double *);
    // retrieve nonredundant integrals
    void get_nonredundant_ints_(double *, double *, int, int, int);
    
  public:
    Int2eLibint2(Integral *,
	       const Ref<GaussianBasisSet>&,
	       const Ref<GaussianBasisSet>&,
	       const Ref<GaussianBasisSet>&,
	       const Ref<GaussianBasisSet>&,
	       size_t storage);
    ~Int2eLibint2();

    /// "clones" this engine, all precomputed data is shallow-copied

    /// the default implementation throws
    virtual Ref<Int2eLibint2> clone();

    /// provides the bounds evaluator for log2_bound
    void bounds(const Ref<Log2Bounds>&);
    /// returns the bounds evaluator for log2_bound
    const Ref<Log2Bounds>& bounds() const { return bounds_; }

    ///  Reports how much storage is actually used at a given time
    size_t storage_used() const { return storage_used_; }

    /// Whether redundant integrals are returned
    int redundant() const { return redundant_; }
    /// Set redundant flag
    void set_redundant(int flag) { redundant_ = flag; }

    /// Whether shells can be permuted
    int permute() const { return permute_; }
    /// Set shell permutation flag
    void set_permute(int flag) { permute_ = flag; }

    /// Evaluate the target quartet of integrals
    virtual void compute_quartet(int *, int*, int*, int*) =0;
    /// Returns the location of the buffer with target integrals
    virtual double *buffer(unsigned int) const =0;
    /// Computes log2 bound. By default Int2eLibint2 implementations are not required to compute bounds
    /// hence usually this will return a large value.
    virtual int log2_bound(int s1, int s2, int s3, int s4);


    Ref<GaussianBasisSet> basis()
    {
      if (bs1_==bs2_ && bs1_ == bs3_ && bs1_ == bs4_) return bs1_;
      return 0;
    }
    Ref<GaussianBasisSet> basis1() { return bs1_; }
    Ref<GaussianBasisSet> basis2() { return bs2_; }
    Ref<GaussianBasisSet> basis3() { return bs3_; }
    Ref<GaussianBasisSet> basis4() { return bs4_; }

  protected:
    /// shallow-copies \c other
    Int2eLibint2(const Int2eLibint2& other);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
