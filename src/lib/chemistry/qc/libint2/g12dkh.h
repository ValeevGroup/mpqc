//
// g12dkh.h
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

#include <limits.h>
#include <stdexcept>

#include <util/ref/ref.h>
#include <util/misc/scexception.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/intparams.h>
#include <chemistry/qc/basis/intdescr.h>
#include <chemistry/qc/libint2/shellpairs.h>
#include <chemistry/qc/basis/fjt.h>
#include <chemistry/qc/libint2/int2e.h>
#include <libint2.h>

#if LIBINT2_SUPPORT_G12DKH
#ifndef _chemistry_qc_libint2_g12dkh_h
#define _chemistry_qc_libint2_g12dkh_h

namespace sc {

class Integral;

/** G12DKHLibint2 is a specialization of Int2eLibint2 that computes two-electron integrals specific
    to relativistic explicitly correlated methods which use Gaussian geminals.

    G12DKHLibint2 can compute integrals with the same geminal factor in bra and ket.
 */
class G12DKHLibint2: public Int2eLibint2 {
  private:
  /** Number of integral types produced. Produces g12p4g12_m_g12t1g12t1 integrals */
    static const int num_te_types_ = TwoBodyIntDescrG12DKH::num_intsets;

    typedef IntParamsG12::PrimitiveGeminal PrimitiveGeminal;
    typedef IntParamsG12::ContractedGeminal ContractedGeminal;
    // the geminal in the bra
    ContractedGeminal geminal_bra_;
    // the geminal in the ket -- must be equal to bra at the moment
    ContractedGeminal geminal_ket_;
    
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
    void g12dkh_quartet_data_(prim_data *Data, double scale, double gamma_bra, double gamma_ket);
    /*--- Compute engines ---*/
    Libint_t Libint_;
    Ref<Fjt> Fm_Eval_;

    class ExpensiveMath {
    public:
      ExpensiveMath();
      double fac[4*LIBINT2_MAX_AM_G12DKH+1];
      double bc[4*LIBINT2_MAX_AM_G12DKH+1][4*LIBINT2_MAX_AM_G12DKH+1];
    };
    ExpensiveMath ExpMath_;
  
  public:
    ///
    G12DKHLibint2(Integral *,
		 const Ref<GaussianBasisSet>&,
		 const Ref<GaussianBasisSet>&,
		 const Ref<GaussianBasisSet>&,
		 const Ref<GaussianBasisSet>&,
		 size_t storage,
		 const ContractedGeminal& gbra
	);
    ~G12DKHLibint2();

    double *buffer(unsigned int t) const {
	return target_ints_buffer_[t];
    }

    static size_t storage_required(const Ref<GaussianBasisSet>& b1,
				   const Ref<GaussianBasisSet>& b2 = 0,
				   const Ref<GaussianBasisSet>& b3 = 0,
				   const Ref<GaussianBasisSet>& b4 = 0);
    
    // evaluate integrals
    void compute_quartet(int*, int*, int*, int*);
};

}

#include <chemistry/qc/libint2/g12dkh_quartet_data.h>

#endif // header guard
#endif // if LIBINT2_SUPPORT_G12DKH

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
