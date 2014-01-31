//
// int1e.h
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

#ifndef _chemistry_qc_libint2_int1e_h
#define _chemistry_qc_libint2_int1e_h

#include <util/ref/ref.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/fjt.h>
#include <libint2/boys.h>

namespace sc {

class Integral;

/** Int1eLibint2 is used by OneBodyIntLibint2 and OneBodyDerivIntLibint2 to
    implement IntegralLibint2. */
class Int1eLibint2: public RefCount {
  protected:
    Integral *integral_;

    Ref<GaussianBasisSet> bs1_;
    Ref<GaussianBasisSet> bs2_;

    Ref<IntParams> operset_params_;
    Ref<EfieldDotVectorData> EdotV_origin_;
    Ref<PointChargeData> Q_origin_;

    bool need_overlap_;
    bool need_coulomb_;

    int ntypes_;                  // Number of integrals held together
                                  // usually 1, 3 for dipole, 6 for quadrupole, etc.

    double *target_ints_buffer_;  // buffer of integrals in a shell doublet
    int max_doublet_size_;

    /*--- Internal scratch ---*/
    double *cart_ints_;           // buffer of cartesian integrals in by-contraction-doublets order
    double *sphharm_ints_;        // buffer of integrals over spherical harmonics in by-contraction-doublets order
    double *tformbuf_;            // temporary scratch used in cart.->sph.harm. transforms
    int max_cart_doublet_size_;
    
    /*--- Pointers to scratch arrays (never used in new[] and delete[]) ---*/
    double *prim_ints_;       // this points to the appropriate location for raw integrals
    double *contr_doublets_;
    double *shell_doublet_;

    /*--- Internally used "interfaces" ---*/
    GaussianShell *int_shell1_;
    GaussianShell *int_shell2_;
    struct {
      double A[3], B[3];
      double AB2;
      int gc1, gc2;
      int p1, p2;
      int am;
    } doublet_info_;
    // Buffers for primitive overlap integrals over 1D Gaussians
    double **OIX_, **OIY_, **OIZ_;
    // Buffer for primitive nuclear attraction integrals
    // Such buffer contains all integrals including intermediates
    // These integrals are nonseparable, hence the first and second indices are composite
    double ***AI0_;
    int indmax_;         // that's the range of the first 2 indices of AI0
    // Compute engines
    void AI_OSrecurs_(double ***AI0, double PA[3], double PB[3],
		      double PC[3], double gamma, int iang, int jang);
    void OI_OSrecurs_(double **OIX, double **OIY, double **OIZ, double PA[3], double PB[3],
		      double gamma, int lmaxi, int lmaxj);
    ::libint2::FmEval_Chebyshev3* Fm_Eval_;
    double* Fm_table_;

    // tasks common to different types of integral evaluation
    void compute_doublet_info_(int, int);
    void zero_buffers_();
    void transform_contrquartets_(double *, double *);
    void sort_contrdoublets_to_shelldoublet_(double *, double *);
    void norm_contrcart1_(double* data);
    template <unsigned int ntypes> void norm_contrcart_(double* data);
    // same tasks, but adapted for the case of several sets of integrals (multipole, etc.)
    void zero_buffers_vec_(const int);
    void transform_contrquartets_vec_(const int, double *, double *);
    void sort_contrdoublets_to_shelldoublet_vec_(const int, double *, double *);
    // specialized versions of compute routines
    void overlap_full_general_();
    void overlap_sameam_general_();
    void kinetic_full_general_();
    void kinetic_sameam_general_();
    void nuclear_full_general_();
    void nuclear_sameam_general_();
    void hcore_full_general_();
    void hcore_sameam_general_();
    void edipole_full_general_();
    void equadrupole_full_general_();
    void p4_full_general_();

    // Utility functions
    double **init_block_(int, int);     // Used to allocate OIXYZ
    void free_block_(double **);
    double ***init_box_(int, int, int);   // Used to allocate AI0;
    void free_box_(double ***);
  
  public:
    Int1eLibint2(Integral *,
            const Ref<GaussianBasisSet>&,
            const Ref<GaussianBasisSet>&,
            int order, bool need_overlap, bool need_coulomb, int ntypes);
    ~Int1eLibint2();

    void set_params(const Ref<IntParams>& p);
    void set_EdotV_origin(const Ref<EfieldDotVectorData>&);
    void set_Q_origin(const Ref<PointChargeData>&);
    Ref<IntParams> params();
    Ref<IntParamsOrigin> multipole_origin();
    Ref<EfieldDotVectorData> EdotV_origin();
    Ref<PointChargeData> Q_origin();

    double *buffer() { return target_ints_buffer_; }
    Ref<GaussianBasisSet> basis() { if (bs1_==bs2_) return bs1_; return 0; }
    Ref<GaussianBasisSet> basis1() { return bs1_; }
    Ref<GaussianBasisSet> basis2() { return bs2_; }

    void kinetic(int ish, int jsh);
    void nuclear(int ish, int jsh);
    void overlap(int ish, int jsh);
    void hcore(int ish, int jsh);
    void edipole(int ish, int jsh);
    void equadrupole(int ish, int jsh);
    void p4(int ish, int jsh);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
