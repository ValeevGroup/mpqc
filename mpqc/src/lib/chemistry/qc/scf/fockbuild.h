//
// fockbuild.h --- a generic Fock matrix builder
//
// Based on code gbuild.h:
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
// Maintainer: SNL
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

#ifndef _chemistry_qc_scf_fockbuild_h
#define _chemistry_qc_scf_fockbuild_h

#ifdef __GNUC__
#pragma interface
#endif

#include <scconfig.h>
#include <util/group/thread.h>
#include <util/group/message.h>
#include <chemistry/qc/basis/integral.h>

namespace sc {

class FockBuildMatrix {
    bool owns_data_; /** The owner of the data must outlive other
                      *  FockBuildMatrix's that still use data: Ownership
                      *  is not transferred. */
    double **blockpointers_; /** Points to the blocks within the data
                              * array. */
    Ref<SCMatrix> rectmat_;
    Ref<SymmSCMatrix> symmmat_;
    int nI_, nJ_;
    int ndata_;
    Ref<GaussianBasisSet> bs1_;
    Ref<GaussianBasisSet> bs2_;

    inline int offset(int I, int J) const {
      if (!symmetric()) {
          return I * nJ_ + J;
        }
      else {
          return (I*(I+1))/2 + J;
        }
    }
    inline int nblock() const {
      if (!symmetric()) {
          return nI_ * nJ_;
        }
      else {
          return (nI_*(nI_+1))/2;
        }
    }
    void data_to_symmat() const;
    void data_to_rectmat() const;
  public:
    FockBuildMatrix();
    FockBuildMatrix(const FockBuildMatrix &f) { make_reference(f); }
    void operator = (const FockBuildMatrix &f) { make_reference(f); }
    ~FockBuildMatrix();
    bool symmetric() const { return symmmat_.nonnull(); }
    // This will average the off diagonal elements in each
    // diagonal block, iff the matrix is symmetrix.
    void fix_diagonal_blocks() const;
    void clear();
    // These allocate data to hold the SCMatrix.  The data is only
    // copied if copy is true, otherwise it is zeroed.
    void scmat_to_data(const Ref<SymmSCMatrix> &m,
                       const Ref<GaussianBasisSet> &b, bool copy);
    void scmat_to_data(const Ref<SCMatrix> &m,
                       const Ref<GaussianBasisSet> &b1,
                       const Ref<GaussianBasisSet> &b2, bool copy);
    void make_reference(const FockBuildMatrix &);
    // This copies the held data back into the SCMatrix.
    void data_to_scmat() const;
    inline double *block(int I, int J) const {
      return blockpointers_[offset(I,J)];
    }
    /// If the data is not owned, then copy it so it will be.
    void copy_data();
    /// Accumulate fbm into this.
    void accum(const FockBuildMatrix &fbm);
};

class FockContribution: public RefCount {
    double nint_;
  public:
    FockContribution();
    virtual ~FockContribution();
    /** This routine does not permute any indices.  The matrix elements are
        contracted with the integrals are F_IJ and P_KL.  If F is
        symmetric, then J>I will be ignored.
        K>=L.
    */
    virtual void contrib_e_J(double factor,
                             int I, int J, int K, int L,
                             int nI, int nJ, int nK, int nL,
                             const double * restrictxx buf) = 0;
    /** This routine does not permute any indices.  The matrix elements
        are contracted with the integrals are F_IK and P_JL.  If F is
        symmetric, then K>I will be ignored.
    */
    virtual void contrib_e_K(double factor,
                             int I, int J, int K, int L,
                             int nI, int nJ, int nK, int nL,
                             const double * restrictxx buf) = 0;
    virtual void contrib_p12_p13p24_J(double factor,
                                      int I, int J, int K, int L,
                                      int nI, int nJ, int nK, int nL,
                                      const double * restrictxx buf) = 0;
    virtual void contrib_p12_p13p24_K(double factor,
                                      int I, int J, int K, int L,
                                      int nI, int nJ, int nK, int nL,
                                      const double * restrictxx buf) = 0;
    virtual void contrib_p34_p13p24_J(double factor,
                                      int I, int J, int K, int L,
                                      int nI, int nJ, int nK, int nL,
                                      const double * restrictxx buf) = 0;
    virtual void contrib_p34_p13p24_K(double factor,
                                      int I, int J, int K, int L,
                                      int nI, int nJ, int nK, int nL,
                                      const double * restrictxx buf) = 0;
    virtual void contrib_p12_p34_J(double factor,
                                   int I, int J, int K, int L,
                                   int nI, int nJ, int nK, int nL,
                                   const double * restrictxx buf) = 0;
    virtual void contrib_p12_p34_K(double factor,
                                   int I, int J, int K, int L,
                                   int nI, int nJ, int nK, int nL,
                                   const double * restrictxx buf) = 0;
    virtual void contrib_p12_J(double factor,
                               int I, int J, int K, int L,
                               int nI, int nJ, int nK, int nL,
                               const double * restrictxx buf) = 0;
    virtual void contrib_p12_K(double factor,
                               int I, int J, int K, int L,
                               int nI, int nJ, int nK, int nL,
                               const double * restrictxx buf) = 0;
    virtual void contrib_p34_J(double factor,
                               int I, int J, int K, int L,
                               int nI, int nJ, int nK, int nL,
                               const double * restrictxx buf) = 0;
    virtual void contrib_p34_K(double factor,
                               int I, int J, int K, int L,
                               int nI, int nJ, int nK, int nL,
                               const double * restrictxx buf) = 0;
    virtual void contrib_p13p24_J(double factor,
                                  int I, int J, int K, int L,
                                  int nI, int nJ, int nK, int nL,
                                  const double * restrictxx buf) = 0;
    virtual void contrib_p13p24_K(double factor,
                                  int I, int J, int K, int L,
                                  int nI, int nJ, int nK, int nL,
                                  const double * restrictxx buf) = 0;
    virtual void contrib_all_J(double factor,
                               int I, int J, int K, int L,
                               int nI, int nJ, int nK, int nL,
                               const double * restrictxx buf) = 0;
    virtual void contrib_all_K(double factor,
                               int I, int J, int K, int L,
                               int nI, int nJ, int nK, int nL,
                               const double * restrictxx buf) = 0;
    virtual Ref<FockContribution> clone() = 0;

    virtual void set_fmat(int i, const Ref<SCMatrix> &) = 0;
    virtual void set_fmat(int i, const Ref<SymmSCMatrix> &) = 0;

    virtual void set_jmat(int i, const Ref<SCMatrix> &) = 0;
    virtual void set_jmat(int i, const Ref<SymmSCMatrix> &) = 0;

    virtual void set_kmat(int i, const Ref<SCMatrix> &) = 0;
    virtual void set_kmat(int i, const Ref<SymmSCMatrix> &) = 0;

    virtual void set_pmat(int i, const Ref<SCMatrix> &) = 0;
    virtual void set_pmat(int i, const Ref<SymmSCMatrix> &) = 0;

    /** Compute the maximum of the density in each block.  The pmax vector
        holds only the unique elements if symmetric is true. */
    virtual signed char *compute_pmax(bool symmetric = false) const = 0;

    /// Replicate Fock matrix data so multiple threads can accumulate.
    virtual void copy() = 0;
    /** Sum the Fock matrix contributions from different threads.  The
     *  passed specialization type must be the same as the specialization
     *  of this.
     */
    virtual void accum(const Ref<FockContribution> &) = 0;
    /// Push the internal Fock matrix data back into the original object.
    virtual void update() = 0;

    double nint() const { return nint_; }
    double &nint() { return nint_; }
};

/** The FockBuilder class requires a class template, Contribution.  The
    GenericContribution class provides much of the infrastructure needed by
    Contribution class writers.

    In addition to the members that GenericContribution provides,
    Contribution classes must have members that actually do the
    work and sum in the contributions.  Due to two electron integral
    permutation symmetry, one integral may make several contributions.
    A contrib member exists for each subgroup of the full permutation
    group.

    <dl>

    <dt>void f()<dd>The f member.

    </dl>

    A lightweight copy CTOR should also exist for Contribution
    specializations.
*/
class GenericFockContribution: public FockContribution {
  protected:
    int nfmat_;     /// the number of Fock matrices
    std::vector<FockBuildMatrix> jmats_;
    std::vector<FockBuildMatrix> kmats_;
    std::vector<bool> k_is_j_;
    int npmat_;     /// the number of density matrices
    std::vector<FockBuildMatrix> pmats_;
    Ref<GaussianBasisSet> f_b1_, f_b2_, p_b1_, p_b2_;
    bool f_b1_equiv_f_b2, p_b1_equiv_p_b2;
    double nint_;

    GenericFockContribution(int nfmat, int npmat,
                            Ref<GaussianBasisSet> &f_b1,
                            Ref<GaussianBasisSet> &f_b2,
                            Ref<GaussianBasisSet> &p_b1,
                            Ref<GaussianBasisSet> &p_b2);

    void pmax_contrib(const FockBuildMatrix &mat,
                      signed char *pmax,
                      bool symmetric) const;

  public:
    double *jmat_block(int i, int I, int J) {
      return jmats_[i].block(I,J);
    }
    bool jmat_symmetric(int i) const { return jmats_[i].symmetric(); }
    double *kmat_block(int i, int I, int J) {
      return kmats_[i].block(I,J);
    }
    bool kmat_symmetric(int i) const { return kmats_[i].symmetric(); }
    const double *pmat_block(int i, int I, int J) {
      return pmats_[i].block(I,J);
    }
    bool pmat_symmetric(int i) const { return pmats_[i].symmetric(); }

    void set_fmat(int i, const Ref<SCMatrix> &);
    void set_fmat(int i, const Ref<SymmSCMatrix> &);

    void set_jmat(int i, const Ref<SCMatrix> &);
    void set_jmat(int i, const Ref<SymmSCMatrix> &);

    void set_kmat(int i, const Ref<SCMatrix> &);
    void set_kmat(int i, const Ref<SymmSCMatrix> &);

    void set_pmat(int i, const Ref<SCMatrix> &);
    void set_pmat(int i, const Ref<SymmSCMatrix> &);

    void copy();
    void accum(const Ref<FockContribution> &);
    void update();

    signed char* compute_pmax(bool symmetric = false) const;

    ~GenericFockContribution();
};

/** The FockBuildThread class is used to actually build the Fock matrix.
    It is used by the FockBuilder class.
 */
class FockBuildThread : public Thread {
  protected:
    Ref<FockContribution> contrib_;
    Ref<ThreadLock> lock_;
    Ref<Integral> integral_;
    double accuracy_;
    Ref<MessageGrp> msg_;
    int nthread_;
    int threadnum_;
    const signed char *pmax_;

    int can_sym_offset(int i, int j) { return (i*(i+1))/2 + j; }
    int gen_sym_offset(int i, int j) {
      if (i>=j) { return can_sym_offset(i,j); }
      else      { return can_sym_offset(j,i); }
    }
  public:
    /// Each thread must be given a unique contribution, c.
    FockBuildThread(const Ref<MessageGrp> &msg,
                    int nthread,
                    int threadnum,
                    const Ref<FockContribution>&c,
                    const Ref<ThreadLock> &lock,
                    const Ref<Integral> &integral,
                    double acc, const signed char *pmax);
};

/** The FockBuildThread class is used to actually build the Fock matrix.
    It is used by the FockBuilder class.
 */
class FockBuildThread_F11_P11 : public FockBuildThread {
    Ref<GaussianBasisSet> basis_;
    Ref<PetiteList> pl_;
  public:
    /// Each thread must be given a unique contribution, c.
    FockBuildThread_F11_P11(const Ref<MessageGrp> &msg,
                            int nthread,
                            int threadnum,
                            const Ref<FockContribution>&c,
                            const Ref<ThreadLock> &lock,
                            const Ref<Integral> &integral,
                            double acc, const signed char *pmax,
                            const Ref<PetiteList> &pl,
                            const Ref<GaussianBasisSet> &basis1,
                            const Ref<GaussianBasisSet> &basis2/*not used*/,
                            const Ref<GaussianBasisSet> &basis3/*not used*/,
                            const Ref<GaussianBasisSet> &basis4/*not used*/);
    void run();
};

/** This is used to build the Fock matrix when none of the
    basis sets are equivalent.
 */
class FockBuildThread_F12_P34 : public FockBuildThread {
    Ref<GaussianBasisSet> basis1_;
    Ref<GaussianBasisSet> basis2_;
    Ref<GaussianBasisSet> basis3_;
    Ref<GaussianBasisSet> basis4_;
    Ref<PetiteList> pl_;

    void run_J();
    void run_K();

  public:
    /// Each thread must be given a unique contribution, c.
    FockBuildThread_F12_P34(const Ref<MessageGrp> &msg,
                            int nthread,
                            int threadnum,
                            const Ref<FockContribution>&c,
                            const Ref<ThreadLock> &lock,
                            const Ref<Integral> &integral,
                            double acc, const signed char *pmax,
                            const Ref<PetiteList> &pl,
                            const Ref<GaussianBasisSet> &basis1,
                            const Ref<GaussianBasisSet> &basis2,
                            const Ref<GaussianBasisSet> &basis3,
                            const Ref<GaussianBasisSet> &basis4);
    void run();
};

/** This is used to build the Fock matrix.
 */
class FockBuildThread_F11_P22 : public FockBuildThread {
    Ref<GaussianBasisSet> basis1_;
    Ref<GaussianBasisSet> basis2_;
    Ref<GaussianBasisSet> basis3_;
    Ref<GaussianBasisSet> basis4_;
    Ref<PetiteList> pl_;

    void run_J();
    void run_K();

  public:
    /// Each thread must be given a unique contribution, c.
    FockBuildThread_F11_P22(const Ref<MessageGrp> &msg,
                            int nthread,
                            int threadnum,
                            const Ref<FockContribution>&c,
                            const Ref<ThreadLock> &lock,
                            const Ref<Integral> &integral,
                            double acc, const signed char *pmax,
                            const Ref<PetiteList> &pl,
                            const Ref<GaussianBasisSet> &basis1,
                            const Ref<GaussianBasisSet> &basis2,
                            const Ref<GaussianBasisSet> &basis3,
                            const Ref<GaussianBasisSet> &basis4);
    void run();
};

/** The FockBuild class works with the FockBuildThread class to generate
    Fock matrices for both closed shell and open shell methods.  It uses a
    helper class, FockContribution, to do the work of forming
    contributions from the density matrices and integrals and placing these
    into the partial Fock matrices (the G matrices).
*/
class FockBuild: public RefCount {
    Ref<FockContribution> contrib_;
    Ref<GaussianBasisSet> b_f1_;
    Ref<GaussianBasisSet> b_f2_;
    Ref<GaussianBasisSet> b_p1_;
    Ref<GaussianBasisSet> b_p2_;
    Ref<MessageGrp> msg_;
    Ref<ThreadGrp> thr_;
    Ref<Integral> integral_;
    double accuracy_;

    typedef FockBuildThread* (*FBT_CTOR)(const Ref<MessageGrp> &msg,
                                         int nthread,
                                         int threadnum,
                                         const Ref<FockContribution>&c,
                                         const Ref<ThreadLock> &lock,
                                         const Ref<Integral> &integral,
                                         double acc, const signed char *pmax,
                                         const Ref<PetiteList> &pl,
                                         const Ref<GaussianBasisSet> &basis1,
                                         const Ref<GaussianBasisSet> &basis2,
                                         const Ref<GaussianBasisSet> &basis3,
                                         const Ref<GaussianBasisSet> &basis4);


    // Build for the any case.  The thread constructing function is passed in.
    void build_generic(FBT_CTOR);

  public:
    /** Create a FockBuild object using b_f1 as the Fock matrix row
        dimension basis, b_f2 as the Fock matrix column dimension basis,
        b_p1 as the density matrix row dimension, and b_p2 as the density
        matrix column dimension.  If b_f2 is not given, then b_f1 is used.
        If b_p1 is not given, then b_f1 is used.  If b_p2 is not given,
        then b_p1 is used.  If the following parameters are not given, then
        the global defaults are used: The msg parameter specifies the
        MessageGrp, thr gives the ThreadGrp, and integral gives the
        Integral.  */
    FockBuild(const Ref<FockContribution> &contrib,
              double acc,
              const Ref<GaussianBasisSet> &b_f1,
              const Ref<GaussianBasisSet> &b_f2 = 0,
              const Ref<GaussianBasisSet> &b_p1 = 0,
              const Ref<GaussianBasisSet> &b_p2 = 0,
              const Ref<MessageGrp> &msg=MessageGrp::get_default_messagegrp(),
              const Ref<ThreadGrp> &thr=ThreadGrp::get_default_threadgrp(),
              const Ref<Integral> &integral=Integral::get_default_integral());
    virtual ~FockBuild();

    /** Contruct the Fock matrices. */
    void build();

    Ref<FockContribution> contrib() { return contrib_; }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
