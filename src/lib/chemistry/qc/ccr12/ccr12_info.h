//
// ccr12_info.h -- common utilities for all CC/CC-R12 methods
//
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: TS & EFV
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

#ifndef _chemistry_qc_ccr12_ccr12_info_h
#define _chemistry_qc_ccr12_ccr12_info_h

#include <string>
#include <vector>
#include <util/misc/compute.h>
#include <util/group/memory.h>
#include <util/group/message.h>
#include <util/group/thread.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/ccr12/tensor.h>
#include <chemistry/qc/ccr12/tensorextrap.h>

namespace sc {

/** CCR12_Info is the compilation of members that are used in CC and CC-R12 methods.
    It should be strictly owned by these classes. Usually initialized by CCR12::compute() */
class CCR12_Info : virtual public RefCount {
  protected:
    const Ref<SCF> ref_;
    const Ref<MemoryGrp>& mem_;

    const std::string theory_;
    const std::string perturbative_;

    // common tensors
    Ref<Tensor> d_f1;
    Ref<Tensor> d_v2;  bool need_w2_;  bool need_w1_;
    Ref<Tensor> d_t1;  bool need_t1_;
    Ref<Tensor> d_t2;  bool need_t2_;
    Ref<Tensor> d_gt2; bool need_gt2_;
    Ref<Tensor> d_t3;  bool need_t3_;
    Ref<Tensor> d_t4;  bool need_t4_;
    Ref<Tensor> d_vr2; bool need_VpA_;// V^gg_ii ; sometimes V with CABS indices is not needed.
    Ref<Tensor> d_vd2; // V^ii_gg
    Ref<Tensor> d_bs2; // B^ii_ii
    Ref<Tensor> d_xs2; // X^ii_ii
    Ref<Tensor> d_ps2; // P^ii_ii
    Ref<Tensor> d_fr2; bool need_FAA_; bool need_F_; // F12^AA_ii --- sometimes F with two CABS indices is not needed.
    Ref<Tensor> d_fd2; // F12^ii_AA

    Ref<Tensor> d_vd2_gen;

    Ref<Tensor> d_qy; // F12 * gt2
    Ref<Tensor> d_qx; // F12 * gt2 with two CABS (for CCSD-R12)
    Ref<Tensor> d_ly; // gl2 * F12^dagger
    Ref<Tensor> d_lx; // gl2 * F12^dagger with two CABS (for CCSD-R12)

    Ref<Tensor> d_lambda1;
    Ref<Tensor> d_lambda2;
    Ref<Tensor> d_glambda2;
    Ref<Tensor> d_lambda3;
    Ref<Tensor> d_lambda4;

    // computes source integrals
    void compute_source_integrals_rhf();
    void compute_source_integrals_uhf();
    void compute_source_integrals_rhf_r12();

    // this is the OrbitalSpace representing the correlated MO space assumed by SMITH
    Ref<OrbitalSpace> corr_space_;  // full space
    Ref<OrbitalSpace> aobs_space_; // alpha spin only, no CABS
    Ref<OrbitalSpace> bobs_space_; // beta spin only, no CABS
    void compute_corr_space();

    // Fgg
    RefSCMatrix F_[NSpinCases1];
    // FgA and FAg
    RefSCMatrix FpA_[NSpinCases1], FAp_[NSpinCases1];
    // source <pp|ERI|pp> integrals
    Ref<DistArray4> pppp_acc_[NSpinCases2];
    // this computes R12 intermediates
    Ref<R12IntEval> r12int_eval_;
    // V_iigg
    RefSCMatrix Vgg_[NSpinCases2];
    // source <pp|ERI|pA> integrals
    Ref<DistArray4> pppA_acc_[NSpinCases2];
    // source <ii|F12|aA> integrals
    Ref<DistArray4> iiaA_acc_[NSpinCases2];


    long irrep_f_;
    long irrep_v_;
    long irrep_t_;
    long irrep_e_;
    long irrep_y_;

    /// offset routines for const tensors
    void offset_f1();
    void offset_v2();
    void offset_vr2();
    void offset_vd2();
    void offset_fr2();
    void offset_fd2();
    void offset_vd2_gen(bool need_cabs, bool need_xx);

    // tiling data
    int maxtilesize_;
    std::vector<long> sym_;
    std::vector<long> range_;
    std::vector<long> spin_;
    std::vector<long> alpha_;
    std::vector<long> offset_;
    long noab_, nvab_, ncab_;

    // global data
    bool restricted_;
    int  nirrep_;
    bool fixed_;

    // orbital data
    int nfzc_, nfzv_;
    int naoa_, naob_;
    int nava_, navb_;
    int nria_, nrib_;
    std::vector<unsigned int> mosyma_,mosymb_;
    std::vector<unsigned int> cabssyma_,cabssymb_;
    const Ref<R12WavefunctionWorld>& r12world_;
    size_t workmemsize_;

    // orbital data after sorted
    std::vector<double> orbital_evl_sorted_;
    std::vector<long> orbital_spin_sorted_;
    std::vector<long> orbital_sym_sorted_;
    std::vector<long> momap_;

    // functions for initialization
    void set_naocc(int i, int j){naoa_=i; naob_=j;};
    void set_navir(int i, int j){nava_=i; navb_=j;};
    void set_ncabs(int i, int j){nria_=i; nrib_=j;};
    void set_mosym(std::vector<unsigned int> i, std::vector<unsigned int> j){mosyma_=i; mosymb_=j;};
    void set_cabssym(std::vector<unsigned int> i, std::vector<unsigned int> j){cabssyma_=i; cabssymb_=j;};
    void set_fixed(bool fixed){fixed_=fixed;};
    void determine_tilesizes();
    int  determine_tilesize_each(int, int, int, std::vector<unsigned int>&, bool, int&, int);
    void determine_maxtilesize(double);
    void print_tile_info();
    void needs();
    void orbital_energies();

    /// utilities for efficient sort_indices(n) to be implemented
    void transpose(  double* dest,const double* source,const int dim1,const int dim2,const double factor);
    void transpose_1(double* dest,const double* source,const int dim1,const int dim2);

    /// local functions for jacobi_t2_and_gt2_ (i.e. a non-iterative MP2-R12 solver)
    void form_ad(Ref<Tensor>& out);
    void form_ca(const Ref<Tensor>&, const Ref<Tensor>&, Ref<Tensor>&);
    void form_adt(const Ref<Tensor>&, const Ref<Tensor>&, Ref<Tensor>&);

    /// B and X intermediate in RefSymmSCMatrix; to be used in a certain class of methods
    RefSymmSCMatrix B_;
    RefSymmSCMatrix X_;
    void retrieve_B_and_X_ii();
    RefSymmSCMatrix B_ip_;
    RefSymmSCMatrix X_ip_;
    void retrieve_B_and_X_ip();

    RefDiagSCMatrix bdiag_;
    RefSCMatrix lmatrix_;

  public:
    CCR12_Info(const Ref<R12WavefunctionWorld>&,const Ref<MemoryGrp>&,size_t,
               const Ref<SCF>,int,int,int,long,long,int,int,
               std::string,std::string,int);
    ~CCR12_Info();

    void print(std::ostream&);

    double magnetic_moment() const { return ref_->magnetic_moment(); }

    /// constants used in initialization
    const Ref<R12IntEval>& r12eval() const { return r12int_eval_; }
    const Ref<R12WavefunctionWorld>& r12world() { return r12world_; }
    const Ref<SCF> ref(){return ref_;};
    bool fixed() const {return fixed_;};
    long nirrep()const {return (long)nirrep_;};
    int nirrep_int() const {return nirrep_;};
    int nfzc() const {return nfzc_;};
    int nfzv() const {return nfzv_;};
    int naoa() const {return naoa_;};
    int nava() const {return nava_;};
    int naob() const {return naob_;};
    int navb() const {return navb_;};
    int nria() const {return nria_;};
    int nrib() const {return nrib_;};
    int mosyma(int i) const {return static_cast<int>(mosyma_[i]);};
    int mosymb(int i) const {return static_cast<int>(mosymb_[i]);};
    int cabssyma(int i) const {return static_cast<int>(cabssyma_[i]);};
    int cabssymb(int i) const {return static_cast<int>(cabssymb_[i]);};


    /// Constants used in amplitude evaluators
    long noab() const {return noab_;};
    long nvab() const {return nvab_;};
    long ncab() const {return ncab_;};
    long nab()  const {return noab_+nvab_+ncab_;};
    bool restricted() const {return restricted_;};

    long maxtilesize() const { return (long)maxtilesize_; };
    long get_sym(long tile)   const {return sym_[tile];};
    long get_range(long tile) const {return range_[tile];};
    long get_spin(long tile)  const {return spin_[tile];};
    long get_offset(long tile)const {return offset_[tile];};
    long get_alpha(long tile) const {return alpha_[tile];};
    double get_orb_energy(long orb) const {return orbital_evl_sorted_[orb];};

    Ref<Tensor> f1() const {return d_f1;};
    Ref<Tensor> v2() const {return d_v2;};
    Ref<Tensor> t1() const {return d_t1;};
    Ref<Tensor> t2() const {return d_t2;};
    Ref<Tensor> gt2()const {return d_gt2;};
    Ref<Tensor> t3() const {return d_t3;};
    Ref<Tensor> t4() const {return d_t4;};

    Ref<Tensor> fr2() const {return d_fr2;};
    Ref<Tensor> fd2() const {return d_fd2;};
    Ref<Tensor> vr2() const {return d_vr2;};
    Ref<Tensor> vd2() const {return d_vd2;};
    Ref<Tensor> xs2() const {return d_xs2;};
    Ref<Tensor> bs2() const {return d_bs2;};
    Ref<Tensor> ps2() const {return d_ps2;};

    Ref<Tensor> vd2_gen() const {return d_vd2_gen; };

    Ref<Tensor> qy() const {return d_qy;};
    Ref<Tensor> qx() const {return d_qx;};
    Ref<Tensor> ly() const {return d_ly;};
    Ref<Tensor> lx() const {return d_lx;};

    Ref<Tensor> lambda1() const {return d_lambda1;};
    Ref<Tensor> lambda2() const {return d_lambda2;};
    Ref<Tensor> glambda2() const {return d_glambda2;};
    Ref<Tensor> lambda3() const {return d_lambda3;};

    const Ref<MemoryGrp>& mem() const {return mem_;};
    long irrep_f() const {return irrep_f_;};
    long irrep_v() const {return irrep_v_;};
    long irrep_t() const {return irrep_t_;};
    long irrep_e() const {return irrep_e_;};
    long irrep_y() const {return irrep_y_;};


    /// Functions used in amplitude evaluators
    void restricted_2(const long, const long, long&, long& ) const;
    void restricted_4(const long, const long, const long, const long, long&, long&, long&, long&) const;
    void restricted_6(const long, const long, const long, const long, const long, const long,
                      long&, long&, long&, long&, long&, long&) const;
    void restricted_8(const long, const long, const long, const long, const long, const long, const long, const long,
                      long&, long&, long&, long&, long&, long&, long&, long&) const;
    void sort_indices0(const double* a,double* b,const double facter) const {*b=facter*(*a);};
    void sort_indices2(const double*, double*, const long, const long,
                       const int, const int, const double) const;
    void sort_indices4(const double*, double*, const long, const long, const long, const long,
                       const int, const int, const int, const int, const double) const;
    void sort_indices6(const double*, double*, const long, const long, const long, const long, const long, const long,
                       const int, const int, const int, const int, const int, const int, const double) const;
    void sort_indices8(const double*, double*, const long, const long, const long, const long, const long, const long, const long, const long,
                       const int, const int, const int, const int, const int, const int, const int, const int, const double) const;
    void sort_indices_acc6(const double*, double*, const long, const long, const long, const long, const long, const long,
                           const int, const int, const int, const int, const int, const int, const double) const;
    void sort_indices_acc8(const double*, double*, const long, const long, const long, const long, const long, const long, const long, const long,
                           const int, const int, const int, const int, const int, const int, const int, const int, const double) const;
    void smith_dgemm(const long,const long,const long,const double,const double*,const long,
                     const double*,const long,const double,double*,const long) const;


    /// Constatnts used in specific (i.e. derived) CC-R12 object
    bool need_t1() const {return need_t1_;};
    bool need_t2() const {return need_t2_;};
    bool need_t3() const {return need_t3_;};
    bool need_t4() const {return need_t4_;};
    bool need_w1() const {return need_w1_;};
    bool need_w2() const {return need_w2_;};
    bool need_gt2()const {return need_gt2_;};
    bool need_FAA()const  {return need_FAA_;};
    bool need_VpA()const  {return need_VpA_;};


    /// Functions used in specific (i.e. derived) CC-R12 object
    double get_e(const Ref<Tensor>&);
    void jacobi_t1(const Ref<Tensor>& r1){jacobi_t1_(r1,d_t1);};
    void jacobi_t1_(const Ref<Tensor>&,Ref<Tensor>&);
    void jacobi_t2(const Ref<Tensor>& r2){jacobi_t2_(r2,d_t2);};
    void jacobi_t2_(const Ref<Tensor>&,Ref<Tensor>&);
    void jacobi_t3(const Ref<Tensor>& r3){jacobi_t3_(r3,d_t3);};
    void jacobi_t3_(const Ref<Tensor>&,Ref<Tensor>&);
    void jacobi_t4(const Ref<Tensor>& r4){jacobi_t4_(r4,d_t4);};
    void jacobi_t4_(const Ref<Tensor>&,Ref<Tensor>&);
    void jacobi_t2_and_gt2(const Ref<Tensor>& r2,const Ref<Tensor>& gr2) {jacobi_t2_and_gt2_(r2,d_t2,gr2,d_gt2);};
    void jacobi_t2_and_gt2_(const Ref<Tensor>&,Ref<Tensor>&,const Ref<Tensor>&,Ref<Tensor>&);

    void jacobi_lambda1(const Ref<Tensor>& lr1) {jacobi_l1_(lr1,d_lambda1);};
    void jacobi_l1_(const Ref<Tensor>&,Ref<Tensor>&);
    void jacobi_lambda2(const Ref<Tensor>& lr2) {jacobi_l2_(lr2,d_lambda2);};
    void jacobi_l2_(const Ref<Tensor>&,Ref<Tensor>&);
/*  void jacobi_gl2(Ref<Tensor>&,Ref<Tensor>&);
    void jacobi_l3(Ref<Tensor>&,Ref<Tensor>&);
    void jacobi_l4(Ref<Tensor>&,Ref<Tensor>&); */

    /// Functions for ddot of the tensors in an accurate way (used for energy evaluation).
    double energy_lagrangian_r2(const Ref<Tensor>&) const;
    double energy_lagrangian_r3(const Ref<Tensor>&) const;


    //// DIIS utilities // could not locate in tensor.h..
    Ref<SCExtrapData>  edata(Ref<Tensor> t){return new TensorExtrapData(t);};
    Ref<SCExtrapError> eerr(Ref<Tensor> t) {return new TensorExtrapError(t);};

    /// offset routines
    void offset_t1(Ref<Tensor>&, bool);
    void offset_t2(Ref<Tensor>&, bool);
    void offset_gt2(Ref<Tensor>&, bool);
    void offset_t3(Ref<Tensor>&, bool);
    void offset_t4(Ref<Tensor>&, bool);
    void offset_e(Ref<Tensor>&);

    void offset_x_gen(Ref<Tensor>& t, const bool need_xx, const bool lprint = false);

    void offset_qy();
    void update_qy();
    void offset_qx();
    void update_qx();
    void offset_ly();
    void update_ly();
    void offset_lx();
    void update_lx();

    void offset_l1(Ref<Tensor>&);
    void offset_l2(Ref<Tensor>&);

    /// guess routine
    void guess_t2(Ref<Tensor>& d_t2_);
    void guess_t2_r12(Ref<Tensor>& d_t2_, Ref<Tensor>& d_gt2_);
    void guess_lambda1(Ref<Tensor>& d_lambda1_);
    void guess_lambda1() {guess_lambda1(d_lambda1);};
    void guess_lambda2(Ref<Tensor>& d_lambda2_);
    void guess_lambda2() {guess_lambda2(d_lambda2);};
    void guess_glambda2(Ref<Tensor>& d_glambda2_);
    void guess_glambda2() {guess_glambda2(d_glambda2);};

    /// fill-in routines
    void fill_in_f1();
    void fill_in_v2();
    void fill_in_iiii();
    void fill_in_vr_and_vd();
    void fill_in_fr_and_fd();
    void fill_in_vd2_gen(bool need_cabs, bool need_xx);


    /// utilities for fill-in routines
    long momap(long i) const {return momap_[i];};

    /// utilities for Lambda contribution in fixed-amp approaches
    void prod_iiii(const Ref<Tensor>&, const Ref<Tensor>&, Ref<Tensor>&, const bool transpose = false);

    /// returns B and X intermediate for perturbative methods etc.
    RefSymmSCMatrix B() { return B_; };
    RefSymmSCMatrix X() { return X_; };
    RefSymmSCMatrix B_ip() { return B_ip_; };
    RefSymmSCMatrix X_ip() { return X_ip_; };

    // returns shared pointers of OrbitalSpace objects
    Ref<OrbitalSpace> corr_space() { return corr_space_; };  // full space

    // used in MP2-R12 updates etc.
    void denom_contraction(const Ref<Tensor>&, Ref<Tensor>&);
    void prediagon(RefDiagSCMatrix& eigvals, RefSCMatrix& eigvecs);
    RefSCMatrix lmatrix() { return lmatrix_; };
    RefDiagSCMatrix bdiag() { return bdiag_; };
};

}

#endif
