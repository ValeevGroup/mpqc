//
// mp2r12_energy.cc
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

#include <ostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cassert>
#include <util/misc/string.h>
#include <util/misc/formio.h>
#include <util/misc/regtime.h>
#include <util/ref/ref.h>
#include <math/scmat/local.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbptr12/mp2r12_energy.h>
#include <math/mmisc/pairiter.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <chemistry/qc/lcao/transform_factory.h>
#include <math/scmat/svd.h>
#include <math/scmat/util.h>
#include <chemistry/qc/mbptr12/r12_amps.h>
#include <util/misc/print.h>
#include <chemistry/qc/lcao/utils.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}

#define USE_INVERT 0


/********************************
 * class R12EnergyIntermediates *
 ********************************/

static ClassDesc R12EnergyIntermediates_cd(typeid(R12EnergyIntermediates),"R12EnergyIntermediates",
                                           1,"virtual public SavableState",0,0,
                                           create<R12EnergyIntermediates>);

R12EnergyIntermediates::R12EnergyIntermediates(const Ref<R12IntEval>& r12eval,
                                               const R12Technology::StandardApproximation stdapp) {
  stdapprox_=stdapp;
  r12eval_=r12eval;
  V_computed_=false;
  X_computed_=false;
  B_computed_=false;
  A_computed_=false;
  T1_cc_computed_=false;
  T2_cc_computed_=false;
  Onerdm_cc_computed_=false;
}

R12EnergyIntermediates::R12EnergyIntermediates(StateIn &si) {
  int stdapprox; si.get(stdapprox); stdapprox_=(R12Technology::StandardApproximation)stdapprox;
  r12eval_ << SavableState::restore_state(si);
  int V_computed; si.get(V_computed); V_computed_=(bool)V_computed;
  int X_computed; si.get(X_computed); X_computed_=(bool)X_computed;
  int B_computed; si.get(B_computed); B_computed_=(bool)B_computed;
  int A_computed; si.get(A_computed); A_computed_=(bool)A_computed;
  for(int i=0; i<NSpinCases2; i++){
    V_[i].restore(si);
    X_[i].restore(si);
    B_[i].restore(si);
    A_[i].restore(si);
  }
  T1_cc_computed_=false;
  T2_cc_computed_=false;
  Onerdm_cc_computed_=false;
}

void R12EnergyIntermediates::save_data_state(StateOut &so) {
  so.put((int)stdapprox_);
  SavableState::save_state(r12eval_.pointer(),so);
  so.put((int)V_computed_);
  so.put((int)X_computed_);
  so.put((int)B_computed_);
  so.put((int)A_computed_);
  for(int i=0; i<NSpinCases2; i++){
    V_[i].save(so);
    X_[i].save(so);
    B_[i].save(so);
    A_[i].save(so);
  }

  T1_cc_computed_=false;
  T2_cc_computed_=false;
  Onerdm_cc_computed_=false;
}

Ref<R12IntEval> R12EnergyIntermediates::r12eval() const {
  return(r12eval_);
}

void R12EnergyIntermediates::set_r12eval(Ref<R12IntEval> &r12eval) {
  r12eval_=r12eval;
}

R12Technology::StandardApproximation R12EnergyIntermediates::stdapprox() const {
  return(stdapprox_);
}

bool R12EnergyIntermediates::V_computed() const {
  return(V_computed_);
}

bool R12EnergyIntermediates::X_computed() const {
  return(X_computed_);
}

bool R12EnergyIntermediates::B_computed() const {
  return(B_computed_);
}

bool R12EnergyIntermediates::A_computed() const {
  return(A_computed_);
}

bool R12EnergyIntermediates::T1_cc_computed() const {
  return(T1_cc_computed_);
}
bool R12EnergyIntermediates::T2_cc_computed() const {
  return(T2_cc_computed_);
}
bool R12EnergyIntermediates::Onerdm_cc_computed() const {
  return(Onerdm_cc_computed_);
}
bool R12EnergyIntermediates::Onerdm_relax_computed() const {
  return(Onerdm_relax_computed_);
}
const RefSCMatrix& R12EnergyIntermediates::get_V(const SpinCase2 &spincase2) const {
  return(V_[spincase2]);
}

void R12EnergyIntermediates::assign_V(const SpinCase2 &spincase2, const RefSCMatrix& V) {
  V_[spincase2]=V;
  V_computed_ = true;
}

const RefSymmSCMatrix& R12EnergyIntermediates::get_X(const SpinCase2 &spincase2) const {
  return(X_[spincase2]);
}

void R12EnergyIntermediates::assign_X(const SpinCase2 &spincase2, const RefSymmSCMatrix& X) {
  X_[spincase2]=X;
  X_computed_ = true;
}

const RefSymmSCMatrix& R12EnergyIntermediates::get_B(const SpinCase2 &spincase2) const {
  return(B_[spincase2]);
}

void R12EnergyIntermediates::assign_B(const SpinCase2 &spincase2, const RefSymmSCMatrix& B) {
  B_[spincase2]=B;
  B_computed_ = true;
}

const RefSCMatrix& R12EnergyIntermediates::get_A(const SpinCase2 &spincase2) const {
  return(A_[spincase2]);
}

void R12EnergyIntermediates::assign_A(const SpinCase2 &spincase2, const RefSCMatrix& A) {
  A_[spincase2]=A;
  A_computed_ = true;
}

const RefSCMatrix& R12EnergyIntermediates::get_T1_cc(const SpinCase1 &spincase1) const {
  return(T1_cc_[spincase1]);
}

void R12EnergyIntermediates::assign_T1_cc(const SpinCase1 &spincase1, const RefSCMatrix& T1_cc) {
  T1_cc_[spincase1]=T1_cc;
  T1_cc_computed_ = true;
}

const Ref<DistArray4>& R12EnergyIntermediates::get_T2_cc(const SpinCase2 &spincase2) const {
  return(T2_cc_[spincase2]);
}

void R12EnergyIntermediates::assign_T2_cc(const SpinCase2 &spincase2, const Ref<DistArray4>& T2_cc) {
  T2_cc_[spincase2]=T2_cc;
  T2_cc_computed_ = true;
}

const RefSCMatrix& R12EnergyIntermediates::get_1rdm_cc(const SpinCase1 &spincase1) const {
  return(Onerdm_cc_[spincase1]);
}

void R12EnergyIntermediates::assign_1rdm_cc(const SpinCase1 &spincase1, const RefSCMatrix& Onerdm_cc) {
  Onerdm_cc_[spincase1] = Onerdm_cc;
  Onerdm_cc_computed_ = true;
}

const RefSCMatrix& R12EnergyIntermediates::get_1rdm_relax(const SpinCase1 &spincase1) const {
  return(Onerdm_relax_[spincase1]);
}

void R12EnergyIntermediates::assign_1rdm_relax(const SpinCase1 &spincase1, const RefSCMatrix& Onerdm_relax) {
  Onerdm_relax_[spincase1] = Onerdm_relax;
  Onerdm_relax_computed_ = true;
}

/*-------------
  MP2R12Energy
 -------------*/
static ClassDesc MP2R12Energy_cd(
  typeid(MP2R12Energy),"MP2R12Energy",2,"virtual public SavableState",
  0, 0, 0);

MP2R12Energy::MP2R12Energy(const Ref<R12EnergyIntermediates>& r12intermediates,
                           bool include_obs_singles,
                           int debug) :
                             r12intermediates_(r12intermediates),
                             r12eval_(r12intermediates->r12eval()),
                             include_obs_singles_(include_obs_singles),
                             debug_(debug>=0 ? debug : 0),
                             evaluated_(false)
{
  init();
}

MP2R12Energy::MP2R12Energy(StateIn& si) : SavableState(si)
{
  r12eval_ << SavableState::restore_state(si);
  r12intermediates_ << SavableState::restore_state(si);

  si.get(include_obs_singles_);
  si.get(debug_);
  si.get(evaluated_);

  init();

  for(int s=0; s<NSpinCases2; s++) {
    ef12_[s].restore(si);
    emp2f12_[s].restore(si);
    C_[s].restore(si);
  }

}

MP2R12Energy::~MP2R12Energy()
{
  r12eval_ = 0;
  r12intermediates_ = 0;
}

void MP2R12Energy::save_data_state(StateOut& so)
{
  SavableState::save_state(r12eval_.pointer(),so);
  SavableState::save_state(r12intermediates_.pointer(),so);

  so.put(include_obs_singles_);
  so.put(debug_);
  so.put(evaluated_);

  for(int s=0; s<NSpinCases2; s++) {
    ef12_[s].save(so);
    emp2f12_[s].save(so);
    C_[s].save(so);
  }
}

void
MP2R12Energy::init()
{
  const Ref<R12WavefunctionWorld> r12world = r12eval_->r12world();
  if(r12world->r12tech()->ansatz()->orbital_product_gg()==R12Technology::OrbProdgg_pq) {
    throw InputError("MP2R12Energy_SpinOrbital::init_ -- pq Ansatz not allowed for MP2-R12.",__FILE__,__LINE__);
  }
  Ref<SCMatrixKit> kit = new LocalSCMatrixKit;
  for(int s=0; s<NSpinCases2; s++) {
    const bool spin_polarized = r12world->refwfn()->spin_polarized();
    if (spin_polarized || s != BetaBeta) {
      RefSCDimension dim_gg = r12eval()->dim_gg(static_cast<SpinCase2>(s));
      RefSCDimension dim_f12 = r12eval()->dim_f12(static_cast<SpinCase2>(s));
      C_[s] = kit->matrix(dim_f12,dim_gg);  C_[s].assign(0.0);
      ef12_[s] = kit->vector(dim_gg);   ef12_[s].assign(0.0);
      emp2f12_[s] = kit->vector(dim_gg);  emp2f12_[s].assign(0.0);
    }
    else {
      C_[BetaBeta] = C_[AlphaAlpha];
      ef12_[BetaBeta] = ef12_[AlphaAlpha];
      emp2f12_[BetaBeta] = emp2f12_[AlphaAlpha];
    }
  }
}

void MP2R12Energy::obsolete()
{
  evaluated_ = false;
}

Ref<R12IntEval> MP2R12Energy::r12eval() const { return r12eval_; };
const Ref<R12EnergyIntermediates>& MP2R12Energy::r12intermediates() const { return(r12intermediates_); };
R12Technology::StandardApproximation MP2R12Energy::stdapprox() const { return(r12intermediates_->stdapprox()); };
void MP2R12Energy::set_debug(int debug) { debug_ = debug; };
int MP2R12Energy::get_debug() const { return debug_; };

static void print_psi_values(std::ostream& fout, const SCVector3& r1, const SCVector3& r2, double phi_0, double phi_1_mp2, double phi_1_r12)
{
  fout << scprintf("%9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %25.15lf %25.15lf %25.15lf",
                   r1.x(),r1.y(),r1.z(),r2.x(),r2.y(),r2.z(),phi_0,phi_1_mp2,phi_1_r12) << endl;
}

void MP2R12Energy::print(std::ostream& so) const
{
}

namespace {
  // Assigns src to dest safely, i.e. by converting to a double*
  void assign(RefSCVector& dest, const RefSCVector& src) {
    const int n = src.dim().n();
    double* buf = new double[n];
    src->convert(buf);
    dest->assign(buf);
    delete[] buf;
  }
}

void
MP2R12Energy::compute() {
  if (evaluated_)
    return;

  const bool spin_polarized = r12eval_->r12world()->refwfn()->spin_polarized();
  const int num_unique_spincases2 = (spin_polarized ? 3 : 2);

  // compute mp2 energies
  for(int s=0; s<num_unique_spincases2; ++s) {
    const SpinCase2 spincase = static_cast<SpinCase2>(s);
    RefSCVector emp2 = r12eval()->emp2(spincase);
    assign(emp2f12_[s], emp2);
  }

  // f12 energies
  const Ref<R12Technology::GeminalDescriptor> gdescr = r12eval_->r12world()->r12tech()->corrfactor()->geminaldescriptor();
  if (!R12Technology::invalid(gdescr)) { // corr_factor = none? nothing to compute
    compute_ef12();
    for(int s=0; s<num_unique_spincases2; ++s) {
      emp2f12_[s].accumulate(ef12_[s]);
    }
  }

  evaluated_ = true;
}

double MP2R12Energy::emp2f12tot(SpinCase2 s) {
  compute();
  RefSCVector unit = emp2f12_[s].clone();
  unit.assign(1.0);
  return emp2f12_[s].dot(unit);
}

double MP2R12Energy::ef12tot(SpinCase2 s) {
  compute();
  RefSCVector unit = ef12_[s].clone();
  unit.assign(1.0);
  return ef12_[s].dot(unit);
}

void MP2R12Energy::print_pair_energies(bool spinadapted,
                                       double emp2_cabs_singles_energy,
                                       std::ostream& so)
{
  // if CABS singles are requested, OBS singles are included in the CABS singles correction
  if (emp2_cabs_singles_energy != 0.0) MPQC_ASSERT(include_obs_singles_ == false);

  compute();

  std::string SA_str;
  switch (stdapprox()) {
    case R12Technology::StdApprox_Ap:  SA_str = "A'";  break;
    case R12Technology::StdApprox_App: SA_str = "A''";  break;
    case R12Technology::StdApprox_B:   SA_str = "B";   break;
    case R12Technology::StdApprox_C:   SA_str = "C";   break;
    case R12Technology::StdApprox_Cp:  SA_str = "C'";   break;
    default:
      throw InputError("MP2R12Energy::print_pair_energies -- stdapprox_ is not valid",
                       __FILE__,__LINE__);
  }

  const Ref<R12WavefunctionWorld> r12world = r12eval_->r12world();
  const double escf = r12world->refwfn()->energy();
  // WARNING assuming only RHF and ROHF
  const bool spin_polarized = r12world->refwfn()->spin_polarized();
  const int num_unique_spincases2 = (spin_polarized ? 3 : 2);

  // only used if spinadapted == true
  double ef12tot_0;
  double ef12tot_1;
  double emp2f12tot_0;
  double emp2f12tot_1;

  /*---------------------------------------
    Spin-adapt pair energies, if necessary
   ---------------------------------------*/
  if (!spinadapted) {
    for(int s=0; s<num_unique_spincases2; s++) {
      SpinCase2 spincase2 = static_cast<SpinCase2>(s);
      const RefSCVector ef12 = ef12_[s];
      const RefSCVector emp2f12 = emp2f12_[s];
      const Ref<OrbitalSpace> gspace1 = r12eval()->ggspace(case1(spincase2));
      const Ref<OrbitalSpace> gspace2 = r12eval()->ggspace(case2(spincase2));
      SpinMOPairIter ij_iter(gspace1->rank(), gspace2->rank(), spincase2);

      so << endl << indent << prepend_spincase(spincase2,"MBPT2-R12/") << SA_str << " pair energies:" << endl;
      so << indent << scprintf("    i       j        mp2(ij)        f12(ij)      mp2-f12(ij)") << endl;
      so << indent << scprintf("  -----   -----   ------------   ------------   ------------") << endl;
      for(ij_iter.start(); ij_iter; ij_iter.next()) {
        const int i = ij_iter.i();
        const int j = ij_iter.j();
        const int ij = ij_iter.ij();
        const double ep_f12 = ef12->get_element(ij);
        const double ep_mp2f12 = emp2f12->get_element(ij);
        const double ep_mp2 = ep_mp2f12 - ep_f12;
        so << indent << scprintf("  %3d     %3d     %12.9lf   %12.9lf   %12.9lf",i+1,j+1,ep_mp2,ep_f12,ep_mp2f12) << endl;
      }

    }
  }
  else {
    MPQC_ASSERT(r12eval()->dim_oo(AlphaBeta) == r12eval()->dim_gg(AlphaBeta));
    Ref<SCMatrixKit> localkit = C_[AlphaAlpha].kit();

    const Ref<OrbitalSpace> gspace = r12eval()->ggspace(Alpha);
    const int nocc = gspace->rank();
    RefSCDimension dim_oo_s = new SCDimension(nocc * (nocc + 1) / 2);
    RefSCDimension dim_oo_t = r12eval_->dim_oo(AlphaAlpha);
    RefSCVector emp2f12_0 = localkit->vector(dim_oo_s);
    RefSCVector emp2f12_1 = localkit->vector(dim_oo_t);
    RefSCVector ef12_0 = localkit->vector(dim_oo_s);
    RefSCVector ef12_1 = localkit->vector(dim_oo_t);

    // Triplet pairs are easy
    assign(emp2f12_1,emp2f12_[AlphaAlpha]);
    emp2f12_1->scale(3.0);
    assign(ef12_1,ef12_[AlphaAlpha]);
    ef12_1->scale(3.0);

    // Singlet pairs are a bit trickier
    const RefSCVector emp2f12_ab = emp2f12_[AlphaBeta];
    const RefSCVector emp2f12_aa = emp2f12_[AlphaAlpha];
    const RefSCVector ef12_ab = ef12_[AlphaBeta];
    const RefSCVector ef12_aa = ef12_[AlphaAlpha];
    SpatialMOPairIter_eq ij_iter(gspace->rank());
    int ij_s = 0;
    for(ij_iter.start(); ij_iter; ij_iter.next(), ++ij_s) {
      const int ij_ab = ij_iter.ij_ab();
      const int ij_aa = ij_iter.ij_aa();
      const int i = ij_iter.i();
      const int j = ij_iter.j();
      {
        double eab = emp2f12_ab->get_element(ij_ab);
        double eaa = 0.0;
        if (ij_aa != -1)
          eaa = emp2f12_aa->get_element(ij_aa);
        double e_s = (i != j ? 2.0 : 1.0) * eab - eaa;
        emp2f12_0->set_element(ij_s,e_s);
      }
      {
        double eab = ef12_ab->get_element(ij_ab);
        double eaa = 0.0;
        if (ij_aa != -1)
          eaa = ef12_aa->get_element(ij_aa);
        double e_s = (i != j ? 2.0 : 1.0) * eab - eaa;
        ef12_0->set_element(ij_s,e_s);
      }
    }
    // compute total singlet and triplet energies
    RefSCVector unit_0 = ef12_0.clone();
    RefSCVector unit_1 = ef12_1.clone();
    unit_0->assign(1.0);
    unit_1->assign(1.0);
    ef12tot_0 = ef12_0.dot(unit_0);
    ef12tot_1 = ef12_1.dot(unit_1);
    emp2f12tot_0 = emp2f12_0.dot(unit_0);
    emp2f12tot_1 = emp2f12_1.dot(unit_1);

    so << endl << indent << "Singlet MBPT2-R12/" << SA_str << " pair energies:" << endl;
    so << indent << scprintf("    i       j        mp2(ij)        r12(ij)      mp2-r12(ij)") << endl;
    so << indent << scprintf("  -----   -----   ------------   ------------   ------------") << endl;
    int ng = gspace->rank();
    for(int i=0,ij=0;i<ng;i++) {
      for(int j=0;j<=i;j++,ij++) {
        const double ep_f12_0 = ef12_0.get_element(ij);
        const double ep_mp2f12_0 = emp2f12_0.get_element(ij);
        so << indent << scprintf("  %3d     %3d     %12.9lf   %12.9lf   %12.9lf",
                                 i+1,j+1,ep_mp2f12_0-ep_f12_0,ep_f12_0,ep_mp2f12_0) << endl;
      }
    }

    so << endl << indent << "Triplet MBPT2-R12/" << SA_str << " pair energies:" << endl;
    so << indent << scprintf("    i       j        mp2(ij)        r12(ij)      mp2-r12(ij)") << endl;
    so << indent << scprintf("  -----   -----   ------------   ------------   ------------") << endl;
    for(int i=0,ij=0;i<ng;i++) {
      for(int j=0;j<i;j++,ij++) {
        const double ep_f12_1 = ef12_1.get_element(ij);
        const double ep_mp2f12_1 = emp2f12_1.get_element(ij);
        so << indent << scprintf("  %3d     %3d     %12.9lf   %12.9lf   %12.9lf",
                                 i+1,j+1,ep_mp2f12_1-ep_f12_1,ep_f12_1,ep_mp2f12_1) << endl;
      }
    }

  }

  const double ef12_corr_energy = ef12tot(AlphaAlpha) + ef12tot(BetaBeta) + ef12tot(AlphaBeta);
  const double emp2_obs_singles_energy = include_obs_singles_ ? r12eval()->emp2_obs_singles() : 0.0;
  const double emp2f12_corr_energy = emp2f12tot(AlphaAlpha) +
                                     emp2f12tot(BetaBeta) +
                                     emp2f12tot(AlphaBeta) +
                                     emp2_obs_singles_energy +
                                     emp2_cabs_singles_energy;

  ///////////////////////////////////////////////////////////////
  // The computation of the MP2 energy is now complete on each
  // node;
  ///////////////////////////////////////////////////////////////

  if (spinadapted) {
    so <<endl<<indent
    <<scprintf("Singlet MP2 correlation energy [au]:           %17.12lf\n", emp2f12tot_0 - ef12tot_0);
    so <<indent
    <<scprintf("Triplet MP2 correlation energy [au]:           %17.12lf\n", emp2f12tot_1 - ef12tot_1);
    so <<indent
    <<scprintf("Singlet (MP2)-R12/%3s correlation energy [au]: %17.12lf\n", SA_str.c_str(), ef12tot_0);
    so <<indent
    <<scprintf("Triplet (MP2)-R12/%3s correlation energy [au]: %17.12lf\n", SA_str.c_str(), ef12tot_1);
    so <<indent
    <<scprintf("Singlet MP2-R12/%3s correlation energy [au]:   %17.12lf\n", SA_str.c_str(),
    emp2f12tot_0);
    so <<indent
    <<scprintf("Triplet MP2-R12/%3s correlation energy [au]:   %17.12lf\n", SA_str.c_str(),
    emp2f12tot_1);
  }

  const double etotal = escf + emp2f12_corr_energy;
  so <<endl<<indent
  <<scprintf("RHF energy [au]:                               %17.12lf\n", escf);
  if (emp2_obs_singles_energy != 0.0)
  so <<indent
  <<scprintf("OBS singles MP2 correlation energy [au]:       %17.12lf\n", emp2_obs_singles_energy);
  if (emp2_cabs_singles_energy != 0.0)
  so <<indent
  <<scprintf("CABS singles MP2 correlation energy [au]:      %17.12lf\n", emp2_cabs_singles_energy);
  so <<indent
  <<scprintf("MP2 correlation energy [au]:                   %17.12lf\n", emp2f12_corr_energy - ef12_corr_energy);
  so <<indent
  <<scprintf("(MBPT2)-R12/%3s correlation energy [au]:       %17.12lf\n", SA_str.c_str(), ef12_corr_energy);
  so <<indent
  <<scprintf("MBPT2-R12/%3s correlation energy [au]:         %17.12lf\n", SA_str.c_str(),
  emp2f12_corr_energy);
  so <<indent
  <<scprintf("MBPT2-R12/%3s energy [au]:                     %17.12lf\n", SA_str.c_str(), etotal) << endl;

  so.flush();

  return;
}

#if MP2R12ENERGY_CAN_COMPUTE_PAIRFUNCTION

void
MP2R12Energy::compute_pair_function(unsigned int i, unsigned int j, SpinCase2 spincase2,
                                    const Ref<TwoBodyGrid>& tbgrid)
{
  const bool spin_polarized = r12eval()->r12world()->refwfn()->spin_polarized();
  const SpinCase2 sc2 = (!spin_polarized && spincase2 == BetaBeta ? AlphaAlpha : spincase2);
  const SpinCase1 spin1 = case1(sc2);
  const SpinCase1 spin2 = case2(sc2);
  const bool p1_neq_p2 = spin_polarized && (spincase2 == AlphaBeta);
  const bool antisymm = (spincase2 != AlphaBeta);

  // Cannot plot same-spin pairs yet
  if (spincase2 != AlphaBeta)
    return;

  // convert replicated matrix to local matrix
  RefSCMatrix C;
  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  {
    RefSCMatrix Crepl = C_[sc2];
    C = localkit->matrix(Crepl.rowdim(),Crepl.coldim());
    double* c = new double[Crepl.rowdim().n()*Crepl.coldim().n()];
    Crepl.convert(c);
    C.assign(c);
    delete[] c;
  }
  // and transpose so that row dimension is for |ij> pairs
  C = C.t();
  if (debug_ >= DefaultPrintThresholds::mostO2N2) C.print("C amplitudes");

  Ref<R12WavefunctionWorld> r12world = r12eval()->r12world();
  Ref<OrbitalSpace> vir1_act = r12eval()->vir_act(spin1);
  Ref<OrbitalSpace> vir2_act = r12eval()->vir_act(spin2);
  Ref<OrbitalSpace> occ1_act = r12eval()->occ_act(spin1);
  Ref<OrbitalSpace> occ2_act = r12eval()->occ_act(spin2);
  Ref<OrbitalSpace> occ1 = r12eval()->occ(spin1);
  Ref<OrbitalSpace> occ2 = r12eval()->occ(spin2);
  Ref<OrbitalSpace> ribs1 = r12world->cabs_space(spin1);
  Ref<OrbitalSpace> ribs2 = r12world->cabs_space(spin2);

  // Pair index
  unsigned int ij;
  switch(spincase2) {
  case AlphaBeta:  ij = i*occ2_act->rank() + j; break;
  case AlphaAlpha:
  case BetaBeta:
  {
    // Cannot violate Pauli exclusion principle
    if (i == j)
      return;
    const unsigned int ii = std::max(i,j);
    const unsigned int jj = std::min(i,j);
    ij = ii*(ii-1)/2 + jj;
  }
  break;
  default: MPQC_ASSERT(false);
  }

  const Ref<R12Technology::CorrelationFactor> corrfactor = r12world->r12tech()->corrfactor();
  const unsigned int nf12 = corrfactor->nfunctions();
  // No same-spin pairs if number of orbitals == 1
  if (spincase2 != AlphaBeta && occ1_act->rank() == 1)
    return;

  Ref<R12Amplitudes> Amps = r12eval_->amps();
  RefSCMatrix T2 = Amps->T2(sc2);  if (debug_ >= DefaultPrintThresholds::mostO2N2) T2.print("T2 amplitudes");
  const unsigned int nij = T2.rowdim().n();
  if (ij >= nij)
    return;
  RefSCMatrix Fvv = Amps->Fvv(sc2);  if (debug_ >= DefaultPrintThresholds::mostO2N2)Fvv.print("F12(vv) matrix");
  RefSCMatrix Foo = Amps->Foo(sc2);  if (debug_ >= DefaultPrintThresholds::mostO2N2)Foo.print("F12(oo) matrix");
  RefSCMatrix Fov = Amps->Fov(sc2);  if (debug_ >= DefaultPrintThresholds::mostO2N2)Fov.print("F12(ov) matrix");
  RefSCMatrix Fox = Amps->Fox(sc2);  if (debug_ >= DefaultPrintThresholds::mostO2N2)Fox.print("F12(ox) matrix");
  RefSCMatrix Fvo, Fxo;
  if (p1_neq_p2) {
    Fvo = Amps->Fvo(sc2);  if (debug_ >= DefaultPrintThresholds::mostO2N2)Fvv.print("F12(vo) matrix");
    Fxo = Amps->Fxo(sc2);  if (debug_ >= DefaultPrintThresholds::mostO2N2)Fvv.print("F12(xo) matrix");
  }

  RefSCMatrix Cvv = C * Fvv;  if (debug_ >= DefaultPrintThresholds::mostO2N2)Cvv.print("C(vv) matrix");
  RefSCMatrix Coo = C * Foo;  if (debug_ >= DefaultPrintThresholds::mostO2N2)Coo.print("C(oo) matrix");
  RefSCMatrix Cov = C * Fov;  if (debug_ >= DefaultPrintThresholds::mostO2N2)Cov.print("C(ov) matrix");
  RefSCMatrix Cox = C * Fox;  if (debug_ >= DefaultPrintThresholds::mostO2N2)Cox.print("C(ox) matrix");
  RefSCMatrix Cvo, Cxo;
  if (p1_neq_p2) {
    Cvo = C * Fvo;  if (debug_ >= DefaultPrintThresholds::mostO2N2) Cvv.print("C(vo) matrix");
    Cxo = C * Fxo;  if (debug_ >= DefaultPrintThresholds::mostO2N2) Cvv.print("C(xo) matrix");
  }

  const int nelem = tbgrid->nelem();
  std::string spinlabel;
  switch(spincase2) {
  case AlphaBeta: spinlabel = "ab"; break;
  case AlphaAlpha: spinlabel = "aa"; break;
  case BetaBeta: spinlabel = "bb"; break;
  default: MPQC_ASSERT(false);
  }
  std::stringstream output_file_name;
  output_file_name << SCFormIO::default_basename() << ".pair_function." << tbgrid->name() << "." << spinlabel << "."
                   << i << "_" << j << ".txt";
  ofstream ofile(output_file_name.str().c_str());

  // get coefficients for ij-th pair
  double* c_ij = new double[nf12*nij];
  {
    RefSCVector Cij = C.get_row(ij);
    Cij.convert(c_ij);
  }
  RefSCVector C_ij_f = localkit->vector(C.rowdim());

  for(int p=0; p<nelem; p++) {
    const SCVector3 r1 = tbgrid->xyz1(p);
    const SCVector3 r2 = tbgrid->xyz2(p);

    RefSCVector phi_aa = compute_2body_values_(antisymm,occ1_act,occ2_act,r1,r2);
    RefSCVector phi_vv = compute_2body_values_(antisymm,vir1_act,vir2_act,r1,r2);
    RefSCVector phi_oo = compute_2body_values_(antisymm,occ1,occ2,r1,r2);
    RefSCVector phi_ov = compute_2body_values_(antisymm,occ1,vir2_act,r1,r2);
    RefSCVector phi_ox = compute_2body_values_(antisymm,occ1,ribs2,r1,r2);
    RefSCVector phi_vo, phi_xo;
    if (p1_neq_p2) {
      phi_vo = compute_2body_values_(antisymm,vir1_act,occ2,r1,r2);
      phi_xo = compute_2body_values_(antisymm,ribs1,occ2,r1,r2);
    }
    else {
      phi_vo = compute_2body_values_(antisymm,occ1,vir2_act,r2,r1);
      phi_xo = compute_2body_values_(antisymm,occ1,ribs2,r2,r1);
    }

    double phi_t2 = T2.get_row(ij).dot(phi_vv);

    const double r12 = (r1-r2).norm();
    double phi_r12 = 0.0;

    for(int f=0; f<nf12; ++f) {
      C_ij_f.assign(c_ij + f*nij);
      phi_r12 += 0.5 * C_ij_f.dot(phi_aa) * corrfactor->value(f,r12);
    }
    phi_r12 -= 0.5 * Cvv.get_row(ij).dot(phi_vv);
    phi_r12 -= 0.5 * Coo.get_row(ij).dot(phi_oo);
    phi_r12 -= 0.5 * Cov.get_row(ij).dot(phi_ov);
    phi_r12 -= 0.5 * Cox.get_row(ij).dot(phi_ox);
    if (p1_neq_p2) {
      phi_r12 -= 0.5 * Cvo.get_row(ij).dot(phi_vo);
      phi_r12 -= 0.5 * Cxo.get_row(ij).dot(phi_xo);
    }
    else {
      phi_r12 -= 0.5 * Cov.get_row(ij).dot(phi_vo);
      phi_r12 -= 0.5 * Cox.get_row(ij).dot(phi_xo);
    }

    print_psi_values(ofile,r1,r2,phi_aa.get_element(ij),phi_t2,phi_r12);
  }

  ofile.close();
}

RefSCVector
MP2R12Energy::compute_2body_values_(bool equiv, const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2,
                                    const SCVector3& r1, const SCVector3& r2) const
{
  const Ref<Integral> ints = r12eval_->r12world()->integral();
  const Ref<GaussianBasisSet> bs1 = space1->basis();
  const Ref<GaussianBasisSet> bs2 = space2->basis();
  ints->set_basis(bs1,bs2);
  GaussianBasisSet::ValueData* vdata1 = new GaussianBasisSet::ValueData(bs1,ints);
  GaussianBasisSet::ValueData* vdata2 = new GaussianBasisSet::ValueData(bs2,ints);

  const bool space1_eq_space2 = (space1 == space2);
  const int nbasis1 = bs1->nbasis();
  const int nbasis2 = bs2->nbasis();
  const int rank1 = space1->rank();
  const int rank2 = space2->rank();

  const int npair = (space1_eq_space2 && equiv) ? rank1*(rank1-1)/2 : rank1*rank2;
  RefSCDimension pairdim = new SCDimension(npair);

  double* values11 = new double[nbasis1];
  double* values12 = new double[nbasis1];
  double* values21 = new double[nbasis2];
  double* values22 = new double[nbasis2];

  bs1->values(r1,vdata1,values11);
  bs1->values(r2,vdata1,values12);
  bs2->values(r1,vdata2,values21);
  bs2->values(r2,vdata2,values22);

  RefSCMatrix ao2mo_1 = space1->coefs().t();
  RefSCMatrix ao2mo_2 = space2->coefs().t();

  Ref<SCMatrixKit> kit = ao2mo_1.kit();
  RefSCVector vals11 = kit->vector(ao2mo_1.coldim());
  RefSCVector vals12 = kit->vector(ao2mo_1.coldim());
  RefSCVector vals21 = kit->vector(ao2mo_2.coldim());
  RefSCVector vals22 = kit->vector(ao2mo_2.coldim());
  vals11.assign(values11);
  vals12.assign(values12);
  vals21.assign(values21);
  vals22.assign(values22);
  delete[] values11;
  delete[] values12;
  delete[] values21;
  delete[] values22;

  RefSCVector movals11 = ao2mo_1 * vals11;
  RefSCVector movals12 = ao2mo_1 * vals12;
  RefSCVector movals21 = ao2mo_2 * vals21;
  RefSCVector movals22 = ao2mo_2 * vals22;

  kit = new LocalSCMatrixKit;
  RefSCVector vals = kit->vector(pairdim);

  MOPairIterFactory PIFactory;
  Ref<SpatialMOPairIter> ij_iter = PIFactory.mopairiter(space1,space2);
  for(ij_iter->start();int(*ij_iter.pointer());ij_iter->next()) {
    const int i = ij_iter->i();
    const int j = ij_iter->j();
    const int ij_aa = ij_iter->ij_aa();
    const int ij_ab = ij_iter->ij_ab();
    const int ij_ba = ij_iter->ij_ba();

    if (equiv) {
      if (ij_aa != -1) {
        const double value = movals11.get_element(i) * movals22.get_element(j) -
          movals12.get_element(i) * movals21.get_element(j);
        vals.set_element(ij_aa,value);
      }
    }
    else {
      const double value = movals11.get_element(i) * movals22.get_element(j);
      vals.set_element(ij_ab,value);
      if (space1_eq_space2 && ij_ab != ij_ba) {
        const double value = movals11.get_element(j) * movals22.get_element(i);
        vals.set_element(ij_ba,value);
      }
    }

  }

  vdata1->~ValueData();
  vdata2->~ValueData();

  return vals;
}

#endif  /* MP2R12ENERGY_CAN_COMPUTE_PAIRFUNCTION */

const RefSCVector&
MP2R12Energy::emp2f12(SpinCase2 s)
{
  compute();
  return emp2f12_[s];
}

const RefSCVector&
MP2R12Energy::ef12(SpinCase2 s)
{
  compute();
  return ef12_[s];
}

double MP2R12Energy::energy()
{
  const double emp2_obs_singles_energy = include_obs_singles_ ? r12eval()->emp2_obs_singles() : 0.0;
  const double value = emp2f12tot(AlphaAlpha) + emp2f12tot(BetaBeta) + emp2f12tot(AlphaBeta) + emp2_obs_singles_energy;
  return value;
}

RefSCMatrix
MP2R12Energy::C(SpinCase2 S)
{
  compute();
  return C_[static_cast<int>(S)];
}

RefSCMatrix MP2R12Energy::T2(SpinCase2 S)
{
  Ref<R12Amplitudes> Amps = r12eval_->amps();
  RefSCMatrix T2mat = Amps->T2(S);
  return(T2mat);
}

/*-------------------------
 * MP2R12Energy_SpinOrbital
  -------------------------*/
static ClassDesc MP2R12Energy_SpinOrbital_cd(
                           typeid(MP2R12Energy_SpinOrbital),"MP2R12Energy_SpinOrbital",1,"public MP2R12Energy",
                           0, 0, create<MP2R12Energy_SpinOrbital>);

MP2R12Energy_SpinOrbital::MP2R12Energy_SpinOrbital(Ref<R12EnergyIntermediates> &r12intermediates,
                                                   bool include_obs_singles,
                                                   int debug) :
  MP2R12Energy(r12intermediates,include_obs_singles,debug) {
}

MP2R12Energy_SpinOrbital::MP2R12Energy_SpinOrbital(StateIn &si) :
  MP2R12Energy(si)
{
}

MP2R12Energy_SpinOrbital::~MP2R12Energy_SpinOrbital()
{
}

void MP2R12Energy_SpinOrbital::save_data_state(StateOut &so){
  MP2R12Energy::save_data_state(so);
}

/*-------------------------
 * MP2R12Energy_Diag
  -------------------------*/
static ClassDesc MP2R12Energy_Diag_cd(
                           typeid(MP2R12Energy_Diag),"MP2R12Energy_Diag",1,"public MP2R12Energy",
                           0, 0, create<MP2R12Energy_Diag>);

MP2R12Energy_Diag::MP2R12Energy_Diag(Ref<R12EnergyIntermediates> &r12intermediates,
                                     bool include_obs_singles,
                                     int debug) :
  MP2R12Energy(r12intermediates,include_obs_singles,debug) {
}

MP2R12Energy_Diag::MP2R12Energy_Diag(StateIn &si) :
  MP2R12Energy(si)
{
}

MP2R12Energy_Diag::~MP2R12Energy_Diag()
{
}

void MP2R12Energy_Diag::save_data_state(StateOut &so){
  MP2R12Energy::save_data_state(so);
}

Ref<MP2R12Energy> sc::construct_MP2R12Energy(Ref<R12EnergyIntermediates> &r12intermediates,
                                             bool include_obs_singles,
                                             int debug,
                                             bool diag) {
  Ref<MP2R12Energy> mp2r12energy = diag
      ? static_cast<MP2R12Energy*>(new MP2R12Energy_Diag(r12intermediates,include_obs_singles,debug))
      : static_cast<MP2R12Energy*>(new MP2R12Energy_SpinOrbital(r12intermediates,include_obs_singles,debug));

  return mp2r12energy;
}

#include <chemistry/qc/mbptr12/mp2r12_energy_util.h>
#include <chemistry/qc/lcao/utils.h>
#include <chemistry/qc/lcao/utils.impl.h>
#include <chemistry/qc/mbptr12/mp2r12_energy_compute.cc>
#include <chemistry/qc/mbptr12/mp2r12_energy_diag.cc>
#include <chemistry/qc/mbptr12/mp2r12_energy_diag2.cc>
#include <chemistry/qc/mbptr12/compute_density_diag.cc>

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
