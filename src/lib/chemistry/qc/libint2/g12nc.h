//
// g12nc.h
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
#include <libint2/boys.h>
#include <chemistry/qc/libint2/core_ints_engine.h>

#ifndef _chemistry_qc_libint2_g12nc_h
#define _chemistry_qc_libint2_g12nc_h

namespace sc {

class Integral;

/** G12NCLibint2 is a specialization of Int2eLibint2 that computes two-electron integrals specific
    to explicitly correlated methods which use Gaussian geminals (formulation without commutators).

    G12NCLibint2 can compute integrals with 1 or 2 geminals. All 2-geminal integrals can be represented as 1-geminals integrals
    of a product of the original 2 geminals. For example, overlap of 2 geminals (g12*g12') is directly reduced to an integral over 1 geminal (G12)
    whose exponent is a sum of exponents of g12 and g12'. The following integrals over 2 geminals are needed:
    <list type="bullet">
      <item> g12*g12' = G12. Returned as r12_0_g12. </item>
      <item> [g12,[t1,g12']] = -2 exp(g12) * exp(g12') r12^2 G12. Returned as g12t1g12. </item>
      <item> Since g12[ti,g12'] - g12'[ti,g12] = [ti,G12] * (exp(g12') - exp(g12))/(exp(g12')+exp(g12)), need to compute
             (exp(g12') - exp(g12))/(exp(g12')+exp(g12)) * G12. Returned as anti_g12g12.</item>
    </list>
 */
class G12NCLibint2: public Int2eLibint2 {
  private:
  /** Number of integral types produced. Produces eri, r12_m1_g12, r12_0_g12, g12t1g12, anti_g12g12 integrals */
    static const int num_te_types_ = TwoBodyIntDescrG12NC::num_intsets;

    typedef IntParamsG12::PrimitiveGeminal PrimitiveGeminal;
    typedef IntParamsG12::ContractedGeminal ContractedGeminal;
    // the geminal in the bra
    ContractedGeminal geminal_bra_;
    // the geminal in the ket (can be null)
    ContractedGeminal geminal_ket_;

    // types of integrals that can be computed
    enum OperType {coulomb, f12_coulomb, f12, f12_2, f12_T1_f12};
    
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

    /// initializes Data for a given otype
    /// if otype -- coulomb gbra and gket can be 0, etc.
    void g12nc_quartet_data_(prim_data *Data, double scale, OperType otype,
                             const ContractedGeminal* gbra,
                             const ContractedGeminal* gket);
    /*--- Compute engines ---*/
    std::vector<Libint_t> Libint_;
    typedef ::libint2::FmEval_Chebyshev7<double> _FmEvalType;
    typedef CoreIntsEngine<_FmEvalType>::Engine FmEvalType;
    Ref<FmEvalType> Fm_Eval_;
    double* Fm_table_;

    class ExpensiveMath {
    public:
      ExpensiveMath();
      double fac[4*LIBINT2_MAX_AM_eri+1];
      double bc[4*LIBINT2_MAX_AM_eri+1][4*LIBINT2_MAX_AM_eri+1];
    };
    ExpensiveMath ExpMath_;
  
  public:
    /// When integrals with 1 geminal are needed, gket should be IntParamsG12::null_geminal
    G12NCLibint2(Integral *,
		 const Ref<GaussianBasisSet>&,
		 const Ref<GaussianBasisSet>&,
		 const Ref<GaussianBasisSet>&,
		 const Ref<GaussianBasisSet>&,
		 size_t storage,
		 const ContractedGeminal& gbra,
		 const ContractedGeminal& gket
	);
    ~G12NCLibint2();

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

#include <chemistry/qc/libint2/g12nc_quartet_data.h>

#endif // header guard

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
