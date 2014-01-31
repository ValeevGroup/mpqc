//
// int1e.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
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

#ifndef _chemistry_qc_int1e_h
#define _chemistry_qc_int1e_h

#include <util/ref/ref.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/fjt.h>
#include <chemistry/qc/intv3/array.h>

namespace sc {

class Integral;

/** Int1eV3 is a class wrapper for the one body part of the C language
    IntV3 library.  It is used by OneBodyIntV3 and OneBodyDerivIntV3 to
    implement IntegralV3. */
class Int1eV3: public RefCount {
  protected:
    Integral *integral_;

    Ref<GaussianBasisSet> bs1_;
    Ref<GaussianBasisSet> bs2_;
    double *fjttable_;
    Ref<FJT> fjt_;
    int bs1_shell_offset_;
    int bs2_shell_offset_;
    int bs1_func_offset_;
    int bs2_func_offset_;
    int bs1_prim_offset_;
    int bs2_prim_offset_;

    // statics from comp_1e.c:
  protected:
    double oo2zeta_a;
    double oo2zeta_b;
    double sMus[3];
    double sTs;
    double xi;
    double A[3];
    double B[3];
    double C[3];
    double ss;
    double PmA[3];
    double PmB[3];
    double PmC[3];
    double zeta;
    double oo2zeta;
    GaussianShell *gshell1, *gshell2;
    int exponent_weighted;
    int scale_shell_result;
    double result_scale_factor;
    int three_center;
    Ref<GaussianBasisSet> third_centers;
    int third_centernum;
    int init_order;
    double *buff;
    double *cartesianbuffer;
    double *cartesianbuffer_scratch;
    int mu;
    IntV3Arraydoublep3 inter;
    IntV3Arraydoublep3 efield_inter;

  protected:
    void accum_shell_1der(
        double *buff, int ish, int jsh,
        Ref<GaussianBasisSet> dercs, int centernum,
        double (Int1eV3::*)(int,int,int,int,int,int,int,int)
        );
    void accum_shell_block_1der(
        double *buff, int ish, int jsh,
        Ref<GaussianBasisSet> dercs, int centernum,
        void (Int1eV3::*shell_block_function)
                                  (int gc1, int a, int gc2, int b,
                                   int gcsize2, int gcoff1, int gcoff2,
                                   double coef, double *buffer)
        );
    double comp_shell_overlap(int gc1, int i1, int j1, int k1,
                              int gc2, int i2, int j2, int k2);
    double comp_prim_overlap(int i1, int j1, int k1,
                             int i2, int j2, int k2);
    double comp_shell_kinetic(int gc1, int i1, int j1, int k1,
                              int gc2, int i2, int j2, int k2);
    double comp_prim_kinetic(int i1, int j1, int k1,
                             int i2, int j2, int k2);
    double comp_shell_nuclear(int gc1, int i1, int j1, int k1,
                              int gc2, int i2, int j2, int k2);
    void accum_shell_efield(double *buff, int ish, int jsh);
    void accum_shell_block_efield(double *buff, int ish, int jsh);
    double comp_prim_nuclear(int i1, int j1, int k1,
                             int i2, int j2, int k2, int m);
    void comp_shell_efield(double *efield,
                           int gc1, int i1, int j1, int k1,
                           int gc2, int i2, int j2, int k2);
    void comp_shell_block_efield(int gc1, int a, int gc2, int b,
                                 int gcsize2, int gcoff1, int gcoff2,
                                 double coef, double *buffer);
    double comp_prim_efield(int xyz, int i1, int j1, int k1,
                            int i2, int j2, int k2, int m);
    void comp_shell_dipole(double* dipole,
                           int gc1, int i1, int j1, int k1,
                           int gc2, int i2, int j2, int k2);
    double comp_prim_dipole(int axis,
                            int i1, int j1, int k1,
                            int i2, int j2, int k2);
    void comp_shell_block_nuclear(int gc1, int a, int gc2, int b,
                                  int gcsize2, int gcoff1, int gcoff2,
                                  double coef, double *buffer);
    void comp_shell_block_p_dot_nuclear_p(int gc1, int a, int gc2, int b,
                                  int gcsize2, int gcoff1, int gcoff2,
                                  double coef, double *buffer);
    void comp_prim_block_nuclear(int a, int b);
    void comp_prim_block_nuclear_build_a(int a, int b, int m);
    void comp_prim_block_nuclear_build_b(int b, int m);
    void comp_prim_block_efield(int a, int b);
    void comp_prim_block_efield_build_a(int a, int b, int m);
    void comp_prim_block_efield_build_b(int b, int m);
    // routines from comp_1e:
  protected:
    void int_accum_shell_overlap_1der(int ish, int jsh,
                                      Ref<GaussianBasisSet> dercs,
                                      int centernum);
    void int_done_1e();
    void int_initialize_1e(int flags, int order);
#if 0
    double int_prim_overlap(shell_t *pshell1, shell_t *pshell2,
                            double *pA, double *pB,
                            int prim1, int prim2,
                            int i1, int j1, int k1,
                            int i2, int j2, int k2);
#endif
    void int_accum_shell_kinetic(int ish, int jsh);
    void int_accum_shell_kinetic_1der(int ish, int jsh,
                                      Ref<GaussianBasisSet> dercs,
                                      int centernum);
    void int_accum_shell_nuclear_1der(int ish, int jsh,
                                      Ref<GaussianBasisSet> dercs,
                                      int centernum);
    void int_accum_shell_nuclear_hfc_1der(int ish, int jsh,
                                          Ref<GaussianBasisSet> dercs,
                                          int centernum);
    void int_accum_shell_nuclear_hf_1der(int ish, int jsh,
                                         Ref<GaussianBasisSet> dercs,
                                         int centernum);
    void int_accum_shell_nuclear_nonhf_1der(int ish, int jsh,
                                            Ref<GaussianBasisSet> dercs,
                                            int centernum);
    void int_accum_shell_efield(int ish, int jsh,
                                double *position);
    void int_accum_shell_point_charge(int ish, int jsh,
                                      int ncharge, const double* charge,
                                      const double*const* position);
    void int_shell_nuclear_hf_1der(int ish, int jsh,
                                   Ref<GaussianBasisSet> dercs,
                                   int centernum);
    void int_shell_nuclear_nonhf_1der(int ish, int jsh,
                                      Ref<GaussianBasisSet> dercs,
                                      int centernum);
    void int_accum_shell_dipole(int ish, int jsh,
                                double *com);

    // from offsets.cc
  protected:
    void int_initialize_offsets1();
    void int_done_offsets1();

    // from tformv3.cc
  protected:
    double *source;
    int nsourcemax;
    // transform implementation functions:
    void transform_init();
    void transform_done();
    void source_space(int nsource);
    void copy_to_source(double *integrals, int nsource);
    void do_transform_1e(Integral *integ,
                         double *integrals,
                         GaussianShell *sh1, GaussianShell *sh2,
                         int chunk);
    void transform_1e(Integral *integ,
                      double *integrals, double *target,
                      GaussianShell *sh1, GaussianShell *sh2, int chunk);
    void accum_transform_1e(Integral *integ,
                            double *integrals, double *target,
                            GaussianShell *sh1, GaussianShell *sh2, int chunk);

    // functions for general use outside of tformv3.cc:
    void transform_1e(Integral*integ,
                      double *integrals, double *target,
                      GaussianShell *sh1, GaussianShell *sh2);
    void accum_transform_1e(Integral*integ,
                            double *integrals, double *target,
                            GaussianShell *sh1, GaussianShell *sh2);
    void transform_1e_xyz(Integral*integ,
                          double *integrals, double *target,
                          GaussianShell *sh1, GaussianShell *sh2);
    void accum_transform_1e_xyz(Integral*integ,
                                double *integrals, double *target,
                                GaussianShell *sh1, GaussianShell *sh2);

  public:
    Int1eV3(Integral *,
            const Ref<GaussianBasisSet>&,
            const Ref<GaussianBasisSet>&,
            int order);
    ~Int1eV3();

    double *buffer() { return buff; }
    Ref<GaussianBasisSet> basis() { if (bs1_==bs2_) return bs1_; return 0; }
    Ref<GaussianBasisSet> basis1() { return bs1_; }
    Ref<GaussianBasisSet> basis2() { return bs2_; }

    void kinetic(int ish, int jsh);
    void nuclear_slow(int ish, int jsh);
    void nuclear(int ish, int jsh);
    void p_dot_nuclear_p(int ish, int jsh);
    void overlap(int ish, int jsh);
    void hcore(int ish, int jsh);
    void efield(int ish, int jsh, const double position[3]);
    void point_charge(int ish, int jsh,
                      int ncharge, const double* charge,
                      const double*const* position);
    void dipole(int ish, int jsh,
                double *com);

    void hcore_1der(int ish, int jsh,
                    int dercs, int centernum);
    void kinetic_1der(int ish, int jsh,
                      int dercs, int centernum);
    void nuclear_1der(int ish, int jsh,
                      int dercs, int centernum);
    void overlap_1der(int ish, int jsh,
                      int dercs, int centernum);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
