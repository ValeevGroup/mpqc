//
// grt.h
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

#ifndef _chemistry_qc_libint2_grt_h
#define _chemistry_qc_libint2_grt_h

#include <limits.h>

#include <util/ref/ref.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/libint2/shellpairs.h>
#include <chemistry/qc/basis/fjt.h>
#include <chemistry/qc/libint2/int2e.h>
#include <libint2.h>

namespace sc {

class Integral;

/** GRTLibint2 is a specialization of Int2eLibint2 that computes two-electron integrals specific to linear R12 methods */
class GRTLibint2: public Int2eLibint2 {
  private:
    // Number of integral types produced
//    static const int num_te_types_ = 4;
#define num_te_types_ 4

    // Storage for target integrals
    double *target_ints_buffer_[num_te_types_];

    /*--- Intermediate scratch arrays (may be used in new[] and delete[]) ---*/
    double *cart_ints_[num_te_types_];       // cartesian integrals, in by-contraction-quartet order
    double *sphharm_ints_;    // transformed integrals, in by-contraction-quartet order
    double *perm_ints_;       // redundant target integrals in shell quartet order, shells permuted

    /*--- Pointers to scratch arrays (never used in new[] and delete[]) ---*/
    double *prim_ints_[num_te_types_];       // this points to the appropriate location for raw integrals
    double *contr_quartets_[num_te_types_];
    double *shell_quartet_[num_te_types_];

    /*--- Precomputed data ---*/
    Ref<ShellPairsLibint2> shell_pairs12_;
    Ref<ShellPairsLibint2> shell_pairs34_;

    /*--- Internally used "interfaces" ---*/
    struct {
      int p12, p34, p13p24;           // flags indicating if functions were permuted
      ShellPairLibint2 *shell_pair12, *shell_pair34;   // Shell pairs corresponding to the original
                                                     // (before permutation) order of shell
      int *op1, *op2, *op3, *op4;     // pointers to the primitive indices in the original order
      /////////// The rest of data has been permuted according to p12, p34, p13p24
      double A[3], B[3], C[3], D[3];
      double AB2, CD2;
      int gc1, gc2, gc3, gc4;
      int p1, p2, p3, p4;
      int am;
    } quartet_info_;
    typedef Libint_t prim_data;
    void grt_quartet_data_(prim_data *Data, double scale);
    /*--- Compute engines ---*/
    Libr12_t Libr12_;
    Ref<FJT> Fm_Eval_;

  public:
    GRTLibint2(Integral *,
	     const Ref<GaussianBasisSet>&,
	     const Ref<GaussianBasisSet>&,
	     const Ref<GaussianBasisSet>&,
	     const Ref<GaussianBasisSet>&,
	     size_t storage);
    ~GRTLibint2();

    double *buffer(TwoBodyOper::type te_type) const {
      if (te_type == TwoBodyOper::eri ||
	  te_type == TwoBodyOper::r12 ||
	  te_type == TwoBodyOper::r12t1 ||
	  te_type == TwoBodyOper::r12t2)
	return target_ints_buffer_[te_type];
      else
	return 0;
    }

    static size_t storage_required(const Ref<GaussianBasisSet>& b1,
				   const Ref<GaussianBasisSet>& b2 = 0,
				   const Ref<GaussianBasisSet>& b3 = 0,
				   const Ref<GaussianBasisSet>& b4 = 0);

    // evaluate integrals
    void compute_quartet(int*, int*, int*, int*);
};

/* Libr12StaticInterface is an initializer class for the static part
   of libr12's interface (one per executable) */
class Libr12StaticInterface {
    bool ready;

    public:
    Libr12StaticInterface() { init_libr12_base(); ready = true; }
    ~Libr12StaticInterface() { ready = false; }
};

}

#include <chemistry/qc/libint2/grt_quartet_data.h>

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
