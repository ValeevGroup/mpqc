//
// int2e.h
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

#ifdef __GNUG__
#pragma interface
#endif

#ifndef _chemistry_qc_intv3_int2e_h
#define _chemistry_qc_intv3_int2e_h

#include <limits.h>

#include <util/ref/ref.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/oint3/build.h>
#include <chemistry/qc/intv3/fjt.h>
#include <chemistry/qc/intv3/types.h>
#include <chemistry/qc/intv3/storage.h>
#include <chemistry/qc/intv3/array.h>
#include <chemistry/qc/intv3/macros.h>

namespace sc {

class Integral;

#define CHECK_INTEGRAL_ALGORITHM 0

/** Int2eV3 is a class wrapper for the two body part of the C language
    IntV3 library.  It is used by TwoBodyIntV3 and TwoBodyDerivIntV3 to
    implement IntegralV3. */
class Int2eV3: public RefCount {
  protected:
    Integral *integral_;

    BuildIntV3 build;
    Ref<IntegralStorer> storer;

    Ref<GaussianBasisSet> bs1_;
    Ref<GaussianBasisSet> bs2_;
    Ref<GaussianBasisSet> bs3_;
    Ref<GaussianBasisSet> bs4_;

    // the permuted bases
    GaussianBasisSet *pbs1_;
    GaussianBasisSet *pbs2_;
    GaussianBasisSet *pbs3_;
    GaussianBasisSet *pbs4_;

    Ref<MessageGrp> grp_;

    int bs1_shell_offset_;
    int bs2_shell_offset_;
    int bs3_shell_offset_;
    int bs4_shell_offset_;
    int bs1_func_offset_;
    int bs2_func_offset_;
    int bs3_func_offset_;
    int bs4_func_offset_;
    int bs1_prim_offset_;
    int bs2_prim_offset_;
    int bs3_prim_offset_;
    int bs4_prim_offset_;

    // statics from vrr.cc
  public:
    enum { STORAGE_CHUNK = 4096 };
  protected:
    struct store_list {
        void* data[STORAGE_CHUNK];
        struct store_list* p;
    };
    typedef struct store_list store_list_t;
    int n_store_last;
    store_list_t* store;
    typedef int (BuildIntV3::*intfunc)();
    intfunc build_routine[4][4][4][4][2];
    /* Offset shell numbers. */
    int osh1, osh2, osh3, osh4;
    /* Offset primitive numbers. */
    int opr1, opr2, opr3, opr4;
    /* Saved initialization parameters used to free data. */
    int saved_am12,saved_am34,saved_ncon;
    /* Stores the length of the inner loop for integral contraction. */
    IntV3Arrayint3 contract_length;

    // statics from hrr.cc
  protected:
    /* The general contraction numbers. */
    int g1,g2,g3,g4;
    /* A[] - B[] */
    double AmB[3];
    /* C[] - D[] */
    double CmD[3];
    int eAB;
    double *buf34;
    double *buf12;
    double *bufshared;

    int redundant_;
    int permute_;

  protected:
    Ref<FJT> fjt_;

    int *int_shell_to_prim;
    IntV3Arraydouble2 int_shell_r;
    IntV3Arraydouble2 int_prim_zeta;
    IntV3Arraydouble2 int_prim_k;
    IntV3Arraydouble2 int_prim_oo2zeta;
    IntV3Arraydouble3 int_prim_p;

    double *int_buffer;
    double *int_derint_buffer;

    Ref<GaussianBasisSet> int_cs1;
    Ref<GaussianBasisSet> int_cs2;
    Ref<GaussianBasisSet> int_cs3;
    Ref<GaussianBasisSet> int_cs4;

    GaussianShell *int_shell1;
    GaussianShell *int_shell2;
    GaussianShell *int_shell3;
    GaussianShell *int_shell4;

    IntV3Arraydoublep2 ****e0f0_con_ints_array;  /* The contr. int. inter. */

    int int_expweight1; // For exponent weighted contractions.
    int int_expweight2; // For exponent weighted contractions.
    int int_expweight3; // For exponent weighted contractions.
    int int_expweight4; // For exponent weighted contractions.

    // These are used to compute two and three center electron repulsion
    // integrals.  int_unit2 is 1 if shell 2 is to have value one everywhere
    // and int_unit4 is 1 if shell4 is to be a unit function.  Otherwise,
    // they should be zero.
    //

    int int_unit2;
    int int_unit4;
    GaussianShell* int_unit_shell;

    int int_integral_storage;
    int int_store1;
    int int_store2;
    int int_derivative_bounds;

    // locals from vrr.cc
  protected:
    void add_store(void *p);
    void free_store();
    void _free_store(store_list_t* s, int n);
    void build_not_using_gcs(int nc1, int nc2, int nc3, int nc4,
                             int minam1, int minam3, int maxam12, int maxam34,
                             int dam1, int dam2, int dam3, int dam4, int eAB);
    void build_using_gcs(int nc1, int nc2, int nc3, int nc4,
                         int minam1, int minam3, int maxam12, int maxam34,
                         int dam1, int dam2, int dam3, int dam4, int eAB);
    void gen_prim_intermediates(int pr1, int pr2, int pr3, int pr4, int am);
    void gen_prim_intermediates_with_norm(int pr1, int pr2, int pr3, int pr4,
                                 int am, double norm);
    void gen_shell_intermediates(int sh1, int sh2, int sh3, int sh4);
    void blockbuildprim(int minam1, int maxam12, int minam3, int maxam34);
    void blockbuildprim_1(int am12min, int am12max, int am34, int m);
    void blockbuildprim_3(int am34min, int am34max, int m);

    // globals from vrr.cc
  protected:
    void int_init_buildgc(int order,
                          int am1, int am2, int am3, int am4,
                          int nc1, int nc2, int nc3, int nc4);
    void int_done_buildgc();
    void int_buildgcam(int minam1, int minam2, int minam3, int minam4,
                       int maxam1, int maxam2, int maxam3, int maxam4,
                       int dam1, int dam2, int dam3, int dam4,
                       int sh1, int sh2, int sh3, int sh4,
                       int eAB);

    // globals from print2e.cc
  protected:
    void int_offset_print(std::ostream &,
                          double *buffer,
                          Ref<GaussianBasisSet> c1, int s1,
                          Ref<GaussianBasisSet> c2, int s2,
                          Ref<GaussianBasisSet> c3, int s3,
                          Ref<GaussianBasisSet> c4, int s4);
    void int_offset_print_n(std::ostream &, double *buffer,
                            int n1, int n2, int n3, int n4,
                            int o1, int o2, int o3, int o4,
                            int e12, int e13e24, int e34);
    void int_print(std::ostream &, double *buffer,
                   Ref<GaussianBasisSet> c1, int s1,
                   Ref<GaussianBasisSet> c2, int s2,
                   Ref<GaussianBasisSet> c3, int s3,
                   Ref<GaussianBasisSet> c4, int s4);
    void int_print_n(std::ostream &, double *buffer,
                     int n1, int n2, int n3, int n4,
                     int e12, int e13e24, int e34);
    void int_print_intermediates(std::ostream &);

    // locals from hrr.cc
  protected:
    void shiftam_12(double *I0100, double *I1000, double *I0000,
                    int am1, int am2, int am3, int am4);
    void shiftam_12eAB(double *I0100, double *I1000, double *I0000,
                       int am1, int am2, int am3, int am4);
    void shiftam_34(double *I0001, double *I0010, double *I0000,
                    int am1, int am2, int am3, int am4);
        
    // globals from hrr.cc
  protected:
    void int_init_shiftgc(int order, int am1, int am2, int am3, int am4);
    void int_done_shiftgc();
    double *int_shiftgcam(int gc1, int gc2, int gc3, int gc4,
                          int tam1, int tam2, int tam3, int tam4, int peAB);

    // locals from init2e.cc
  protected:
    void alloc_inter(int nprim,int nshell);
    void compute_shell_1(Ref<GaussianBasisSet> cs, int, int);
    void compute_prim_2(Ref<GaussianBasisSet> cs1,int,int,
                        Ref<GaussianBasisSet> cs2,int,int);


    // globals from init2e.cc
  protected:
    double *int_initialize_erep(size_t storage, int order,
                                const Ref<GaussianBasisSet> &cs1,
                                const Ref<GaussianBasisSet> &cs2,
                                const Ref<GaussianBasisSet> &cs3,
                                const Ref<GaussianBasisSet> &cs4);
    void int_done_erep();

    // from tformv3.cc
  protected:
    double *source;
    double *target;
    double *scratch;
    int nsourcemax;
    // transform implementation functions:
    void transform_init();
    void transform_done();
    void source_space(int nsource);
    void copy_to_source(double *integrals, int nsource);
    void do_gencon_sparse_transform_2e(Integral*integ,
                                       double *integrals, double *target,
                                       int index,
                                       GaussianShell *sh1, GaussianShell *sh2,
                                       GaussianShell *sh3, GaussianShell *sh4);
    // functions for general use outside of tformv3.cc:
    // integrals and target may overlap
    void transform_2e_slow(Integral *,
                      double *integrals, double *target,
                      GaussianShell *sh1, GaussianShell *sh2,
                      GaussianShell *sh3, GaussianShell *sh4);
    void transform_2e(Integral *,
                      double *integrals, double *target,
                      GaussianShell *sh1, GaussianShell *sh2,
                      GaussianShell *sh3, GaussianShell *sh4);

    // locals from comp2e.cc
  protected:
    void compute_erep(int flags, int *psh1, int *psh2, int *psh3, int *psh4,
                      int dam1, int dam2, int dam3, int dam4);
    void compute_erep_1der(int flags, double *buffer,
                           int *psh1, int *psh2, int *psh3, int *psh4,
                           int dercenter);
    void nonredundant_erep(double *buffer, int e12, int e34, int e13e24,
                           int n1, int n2, int n3, int n4,
                           int *red_off, int *nonred_off);
    void compute_erep_bound1der(int flags, double *buffer,
                                int *psh1, int *psh2, int *psh3, int *psh4);

    // globals from comp2e.cc
  protected:
    void int_erep_bound1der(int flags, int bsh1, int bsh2, int *size);


    // global vars from bounds.h
  protected:
    typedef signed char int_bound_t;
    enum { int_bound_min = SCHAR_MIN, int_bound_max = SCHAR_MAX };
    int_bound_t int_Q;    
    int_bound_t int_R;    
    int_bound_t *int_Qvec;
    int_bound_t *int_Rvec;

    // global routines from bounds.cc
  protected:
    void int_init_bounds_nocomp();
    void int_init_bounds_1der_nocomp();
    void int_bounds_comp(int s1, int s2);
    void int_bounds_1der_comp(int s1, int s2);
    int int_erep_2bound(int s1, int s2);
    int int_erep_0bound_1der();
    int int_erep_2bound_1der(int s1, int s2);

    // local routines from bounds.cc
  protected:
    void compute_bounds(int_bound_t *overall, int_bound_t *vec, int flag);
    void compute_bounds_shell(int_bound_t *overall, int_bound_t *vec,
                              int flag, int sh1, int sh2);

    // global routines from storage.cc
  protected:
    int int_have_stored_integral(int sh1,int sh2,int sh3,int sh4,
                                 int p12,int p34,int p13p24);
    void int_store_integral(int sh1,int sh2,int sh3,int sh4,
                            int p12,int p34,int p13p24,
                            int size);

    // from offsets.cc
  protected:
    void int_initialize_offsets2();
    void int_done_offsets2();

    // from comp2e3c.cc
  protected:
    void make_int_unit_shell();
    void delete_int_unit_shell();

  protected:
    // for intermediate storage:
    int used_storage_;
    int used_storage_build_;
    int used_storage_shift_;

  public:
    // bs4 must be null for 3 center integrals
    // bs2 must be null for 2 center integrals
    // bs1 and bs3 must be nonnull.
    Int2eV3(Integral *,
            const Ref<GaussianBasisSet>&bs1,
            const Ref<GaussianBasisSet>&bs2,
            const Ref<GaussianBasisSet>&bs3,
            const Ref<GaussianBasisSet>&bs4,
            int order, size_t storage);
    ~Int2eV3();

    // storage.cc: for the storage of integrals
    void init_storage(int size);
    void done_storage();

    // for intermediate storage
    int storage_used() { return used_storage_; }

    // bounds.cc
    void init_bounds();
    void init_bounds_1der();
    void done_bounds();
    void done_bounds_1der();
    // Covert a bound to/from the log of the bound (returns 2^bound)
    // replace:
    //double int_bound_to_double(int bound);
    //double int_bound_double(int value);
    //int int_bound_log(double value);
    static double logbound_to_bound(int);
    static int bound_to_logbound(double value);

    // If redundant is false the redundant integrals that arise
    // when a shell index is repeated are stored.
    // The default is true.
    int redundant() { return redundant_; }
    void set_redundant(int i) { redundant_ = i; }

    // If permute is true the routines are allowed to permute indices.
    // The default is false.
    int permute() { return permute_; }
    void set_permute(int i) { permute_ = i; }

    int used_storage() const { return used_storage_; }

    // from comp2e.cc
    void erep(int &psh1, int &psh2, int &psh3, int &psh4);
    void erep(int *shells, int  *sizes);
    void erep_all1der(int &psh1, int &psh2, int &psh3, int &psh4,
                      der_centersv3_t *der_centers);
    void erep_all1der(int *shells, int  *sizes,
                      der_centersv3_t *dercenters);

    // from comp2e3c.cc
    void erep_2center(int &psh1, int &psh2);
    void erep_2center(int *shells, int  *sizes);
    void erep_3center(int &psh1, int &psh2, int &psh3);
    void erep_3center(int *shells, int  *sizes);

    // from bounds.cc
    int erep_4bound(int s1, int s2, int s3, int s4);
    int erep_4bound_1der(int s1, int s2, int s3, int s4);

    double *buffer() { return int_buffer; }

    Ref<GaussianBasisSet> basis()
    {
      if (bs1_==bs2_ && bs1_ == bs3_ && bs1_ == bs4_) return bs1_;
      return 0;
    }
    Ref<GaussianBasisSet> basis1() { return bs1_; }
    Ref<GaussianBasisSet> basis2() { return bs2_; }
    Ref<GaussianBasisSet> basis3() { return bs3_; }
    Ref<GaussianBasisSet> basis4() { return bs4_; }

    Ref<GaussianBasisSet> cs1() const { return int_cs1; }
    Ref<GaussianBasisSet> cs2() const { return int_cs2; }
    Ref<GaussianBasisSet> cs3() const { return int_cs3; }
    Ref<GaussianBasisSet> cs4() const { return int_cs4; }

    GaussianBasisSet * pcs1() const { return int_cs1.pointer(); }
    GaussianBasisSet * pcs2() const { return int_cs2.pointer(); }
    GaussianBasisSet * pcs3() const { return int_cs3.pointer(); }
    GaussianBasisSet * pcs4() const { return int_cs4.pointer(); }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
