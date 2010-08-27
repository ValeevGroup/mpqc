//
// mp2r12_energy.h
//
// Copyright (C) 2003 Edward Valeev
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

#ifdef __GNUG__
#pragma interface
#endif

#include <util/ref/ref.h>
#include <chemistry/qc/mbptr12/r12technology.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/spin.h>
#include <chemistry/qc/mbptr12/twobodygrid.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/mp2r12_energy_util.h>

#ifndef _chemistry_qc_mbptr12_mp2r12energy_h
#define _chemistry_qc_mbptr12_mp2r12energy_h

#define MP2R12ENERGY_CAN_COMPUTE_PAIRFUNCTION 1

namespace sc {
  /**
   * The class R12EnergyIntermediates is the front-end to
   * R12 intermediates.
   */
  class R12EnergyIntermediates : virtual public SavableState {
    private:
      R12Technology::StandardApproximation stdapprox_;
      Ref<R12IntEval> r12eval_;
      bool V_computed_;
      RefSCMatrix V_[NSpinCases2];
      bool X_computed_;
      RefSymmSCMatrix X_[NSpinCases2];
      bool B_computed_;
      RefSymmSCMatrix B_[NSpinCases2];
      bool A_computed_;
      RefSCMatrix A_[NSpinCases2];
    public:
      typedef enum { V=0, X=1, B=2, A=3 } IntermediateType;
      R12EnergyIntermediates(const Ref<R12IntEval>& r12eval,
                             const R12Technology::StandardApproximation stdapp);
      R12EnergyIntermediates(StateIn &si);
      void save_data_state(StateOut &so);
      ~R12EnergyIntermediates(){}
      Ref<R12IntEval> r12eval() const;
      void set_r12eval(Ref<R12IntEval> &r12eval);
      R12Technology::StandardApproximation stdapprox() const;
      bool V_computed() const;
      bool X_computed() const;
      bool B_computed() const;
      bool A_computed() const;
      const RefSCMatrix& get_V(const SpinCase2 &spincase2) const;
      void assign_V(const SpinCase2 &spincase2, const RefSCMatrix& V);
      const RefSymmSCMatrix& get_X(const SpinCase2 &spincase2) const;
      void assign_X(const SpinCase2 &spincase2, const RefSymmSCMatrix& X);
      const RefSymmSCMatrix& get_B(const SpinCase2 &spincase2) const;
      void assign_B(const SpinCase2 &spincase2, const RefSymmSCMatrix& B);
      const RefSCMatrix& get_A(const SpinCase2 &spincase2) const;
      void assign_A(const SpinCase2 &spincase2, const RefSCMatrix& A);
  };

  /** Class MP2R12Energy is the object that computes and maintains MP2-R12 energies */
class MP2R12Energy : virtual public SavableState {
  private:
    /** Computes values of all 2-body products from
        space1 and space2 if electron 1 is at r1 and
        electron 2 is at r2. equiv specifies whether electrons
        are equivalent (same spin) or not */
    RefSCVector compute_2body_values_(bool equiv, const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2,
                                      const SCVector3& r1, const SCVector3& r2) const;

    // computes f12 contributions to the mp2 pair energies
    virtual void compute_ef12() =0;

  protected:
    Ref<R12EnergyIntermediates> r12intermediates_;
    Ref<R12IntEval> r12eval_;
    int debug_;
    bool evaluated_;

    RefSCVector ef12_[NSpinCases2], emp2f12_[NSpinCases2];
    // The coefficients are stored xy by ij, where xy is the geminal-multiplied pair
    RefSCMatrix C_[NSpinCases2];

    // Initialize SCVectors and SCMatrices
    void init();

  public:

    MP2R12Energy(StateIn&);
    MP2R12Energy(const Ref<R12EnergyIntermediates>& r12intermediates,
                 int debug);
    ~MP2R12Energy();

    void save_data_state(StateOut&);
    void obsolete();
    void print(std::ostream&o=ExEnv::out0()) const;

    Ref<R12IntEval> r12eval() const;
    const Ref<R12EnergyIntermediates>& r12intermediates() const;
    R12Technology::StandardApproximation stdapprox() const;
    void set_debug(int debug);
    int get_debug() const;

    /// total correlation energy (including OBS singles for non-Brillouin references)
    double energy();
    const RefSCVector& ef12(SpinCase2 S);
    const RefSCVector& emp2f12(SpinCase2 S);
    double emp2f12tot(SpinCase2 S);
    double ef12tot(SpinCase2 S);
    void print_pair_energies(bool spinadapted, std::ostream&so=ExEnv::out0());

#if MP2R12ENERGY_CAN_COMPUTE_PAIRFUNCTION
    /** Computes values of pair function ij on tbgrid */
    void compute_pair_function(unsigned int i, unsigned int j, SpinCase2 spincase2,
                               const Ref<TwoBodyGrid>& tbgrid);
#endif

    /// Computes the first-order R12 wave function and MP2-R12 energy
    void compute();
    /** Returns the matrix of first-order amplitudes of r12-multiplied occupied orbital pairs.
      */
    RefSCMatrix C(SpinCase2 S);
    /** Returns the matrix of first-order amplitudes of conventional orbital pairs.
     */
    RefSCMatrix T2(SpinCase2 S);
};

/**
 * The class MP2R12Energy_SpinOrbital is the original implementation of MP2R12Energy
 * It supports only the standard orbital-invariant ansatz and the full set of features
 * of R12Technology.
 */
class MP2R12Energy_SpinOrbital : public MP2R12Energy
{
    void compute_ef12();

  public:
    MP2R12Energy_SpinOrbital(StateIn&);
    MP2R12Energy_SpinOrbital(Ref<R12EnergyIntermediates> &r12intermediates, int debug);
    ~MP2R12Energy_SpinOrbital();

    void save_data_state(StateOut&);
};

/**
 * The class MP2R12Energy_SpinOrbital_new is a new version of MP2R12Energy_SpinOrbital
 * and computes diagonal and non diagonal MP2 F12 energies preserving the Spin symmetry
 * of the wavefunction. It also implements fixed-amplitude
 * version of MP2-R12 where the coefficients are determined according to the singlet and
 * triplet electron pair coalescence conditions. For non-diagonal MP2 F12 calculations
 * the old version MP2R12Energy_SpinOrbital should be preferred, since it invokes many
 * security checks still not implemented in this version of MP2R12Energy_SpinOrbital_new.
 */
class MP2R12Energy_SpinOrbital_new : public MP2R12Energy
{
  private:
    void compute_ef12();

    /// \param include_coupling_in_B if true, include A^t . D2^{-1} . A term into B
    ///        (this is only appropriate when optimizing amplitudes)
    RefSymmSCMatrix compute_B_non_pairspecific(const RefSymmSCMatrix &B,
                                               const RefSymmSCMatrix &X,
                                               const RefSCMatrix &V,
                                               const RefSCMatrix &A,
                                               const SpinCase2 &spincase2,
                                               bool include_coupling_in_B = true);
    RefSymmSCMatrix compute_B_pairspecific(const SpinMOPairIter &ij_iter,
                                           const RefSymmSCMatrix &B,
                                           const RefSymmSCMatrix &X,
                                           const RefSCMatrix &V,
                                           const RefSCMatrix &A,
                                           const SpinCase2 &spincase2);
    void determine_C_non_pairspecific(const RefSymmSCMatrix &B_ij,
                                      const RefSCMatrix &V,
                                      const SpinCase2 &spincase2,
                                      const Ref<MP2R12EnergyUtil_Diag> &util);
    void determine_C_pairspecific(const int ij,
                                  const RefSymmSCMatrix &B_ij,
                                  const RefSCMatrix &V,
                                  const SpinCase2 &spincase2);
    void determine_C_fixed_non_pairspecific(const SpinCase2 &spincase2);
    void determine_ef12_hylleraas(const RefSymmSCMatrix &B_ij,
                                  const RefSCMatrix &V,
                                  const SpinCase2 &spincase2,
                                  const Ref<MP2R12EnergyUtil_Diag> &util);
    void compute_MP2R12_nondiag();
    void compute_MP2R12_diag_fullopt();
    void compute_MP2R12_diag_nonfullopt();

    // shortcuts
    bool diag() const;
    bool fixedcoeff() const;

  public:
    MP2R12Energy_SpinOrbital_new(StateIn&);
    MP2R12Energy_SpinOrbital_new(Ref<R12EnergyIntermediates> &r12intermediates,
                                 int debug);
    ~MP2R12Energy_SpinOrbital_new();

    void save_data_state(StateOut&);
};

/**
 * The class MP2R12Energy_Diag is an implementation of MP2R12Energy that supports
 * Ten-no's diagonal orbital-invariant ansatz for closed and open-shells.
 * This class will supercede MP2R12Energy_SpinOrbital_new.
 */
class MP2R12Energy_Diag : public MP2R12Energy
{
    void compute_ef12();
    void activate_ints(const std::string&, const std::string&,
                       const std::string&, const std::string&,
                       const std::string&, Ref<TwoBodyFourCenterMOIntsRuntime>&,
                       Ref<DistArray4>&);
    // Label Y^ij_ij, Y^ij_ji, Y^ji_ij, Y^ji_ji
    enum idx_b1b2k1k2{ij_ij, ij_ji, ji_ij, ji_ji};
    void compute_Y(const int b1b2_k1k2, const double prefactor,
                   const unsigned int oper_idx,
                   Ref<DistArray4>& i1i2i1i2_ints, double* array_i1i2i1i2);
    void compute_YxF(const int b1b2_k1k2, const double prefactor,
                     const unsigned int oper1_idx, const unsigned int oper2_idx,
                     Ref<DistArray4>& i1i2x1x2_ints, Ref<DistArray4>& i2i1x1x2_ints,
                     double* array_ijij);
    void compute_FxT(const int b1b2_k1k2, const unsigned int f12_idx,
                     Ref<DistArray4>& F_ints, const double* Tiiaa,
                     double* V_coupling);
    void compute_VX(const int b1b2_k1k2, std::vector<std::string>& VX_output,
                    const unsigned int oper12_idx, Ref<DistArray4>& Y12_ints,
                    const unsigned int oper1_idx, const unsigned int oper2_idx,
                    std::vector<Ref<DistArray4> >& Y_ints,
                    std::vector<Ref<DistArray4> >& F_ints,
                    double* VX_array);
    void accumulate_P_YxF(std::vector<std::string>& P_output,
                          std::vector<int>& b1b2_k1k2, std::vector<double>& P_prefactor,
                          const unsigned int oper1_idx, const unsigned int oper2_idx,
                          std::vector<Ref<DistArray4> >& Y_ints,
                          std::vector<Ref<DistArray4> >& F_ints,
                          double* P);
  public:
    MP2R12Energy_Diag(StateIn&);
    MP2R12Energy_Diag(Ref<R12EnergyIntermediates> &r12intermediates, int debug);
    ~MP2R12Energy_Diag();

    void save_data_state(StateOut&);
};

Ref<MP2R12Energy> construct_MP2R12Energy(Ref<R12EnergyIntermediates> &r12intermediates,
                                         int debug,
                                         bool use_new_version);

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


