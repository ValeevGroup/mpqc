//
// mp2r12_energy.cc
//
// Copyright (C) 2003 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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
#pragma implementation
#endif

#include <ostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <util/misc/string.h>
#include <util/misc/formio.h>
#include <util/misc/timer.h>
#include <util/ref/ref.h>
#include <math/scmat/local.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbptr12/mp2r12_energy.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/transform_factory.h>
#include <chemistry/qc/mbptr12/svd.h>

using namespace std;
using namespace sc;
using namespace sc::exp;

inline int max(int a,int b) { return (a > b) ? a : b;}

#define USE_INVERT 0

/*-------------
  MP2R12Energy
 -------------*/
static ClassDesc MP2R12Energy_cd(
  typeid(MP2R12Energy),"MP2R12Energy",2,"virtual public SavableState",
  0, 0, create<MP2R12Energy>);

MP2R12Energy::MP2R12Energy(Ref<R12IntEval>& r12eval, LinearR12::StandardApproximation stdapp, int debug)
{
  r12eval_ = r12eval;
  stdapprox_ = stdapp;
  if (debug >= 0)
    debug_ = debug;
  else
    debug_ = 0;
  evaluated_ = false;

  init_();
}

void
MP2R12Energy::init_()
{
  const Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
  RefSCDimension dim_oo_aa = r12eval_->dim_oo_aa();
  RefSCDimension dim_oo_ab = r12eval_->dim_oo_ab();
  Ref<SCMatrixKit> kit = r12info->matrixkit();

  RefSCDimension dim_vv_aa = r12eval_->dim_vv_aa();
  RefSCDimension dim_vv_ab = r12eval_->dim_vv_ab();
  for(int s=0; s<NSpinCases2; s++) {
    const bool spin_polarized = r12info->refinfo()->ref()->spin_polarized();
    if (spin_polarized || s != BetaBeta) {
      RefSCDimension dim_oo = r12eval()->dim_oo(static_cast<SpinCase2>(s));
      RefSCDimension dim_f12 = r12eval()->dim_f12(static_cast<SpinCase2>(s));
      C_[s] = kit->matrix(dim_f12,dim_oo);
      ef12_[s] = kit->vector(dim_oo);
      emp2f12_[s] = kit->vector(dim_oo);
    }
    else {
      C_[BetaBeta] = C_[AlphaAlpha];
      ef12_[BetaBeta] = ef12_[AlphaAlpha];
      emp2f12_[BetaBeta] = emp2f12_[AlphaAlpha];
    }
  }
} 

MP2R12Energy::MP2R12Energy(StateIn& si) : SavableState(si)
{
  r12eval_ << SavableState::restore_state(si);

  init_();

  if (si.version(::class_desc<MP2R12Energy>()) >= 2) {
    for(int s=0; s<NSpinCases2; s++) {
      ef12_[s].restore(si);
      emp2f12_[s].restore(si);
      C_[s].restore(si);
    }
  }
  
  int stdapprox;
  si.get(stdapprox);
  stdapprox_ = (LinearR12::StandardApproximation) stdapprox;
  si.get(debug_);
  int evaluated;
  si.get(evaluated);
  evaluated_ = (bool) evaluated;
}

MP2R12Energy::~MP2R12Energy()
{
  r12eval_ = 0;
}

void MP2R12Energy::save_data_state(StateOut& so)
{
  SavableState::save_state(r12eval_.pointer(),so);

  for(int s=0; s<NSpinCases2; s++) {
    ef12_[s].save(so);
    emp2f12_[s].save(so);
    C_[s].save(so);
  }
  
  so.put((int)stdapprox_);
  so.put(debug_);
  so.put((int)evaluated_);
}

void MP2R12Energy::obsolete()
{
  evaluated_ = false;
}

Ref<R12IntEval> MP2R12Energy::r12eval() const { return r12eval_; };
const bool MP2R12Energy::ebc() const { return ebc_; };
const bool MP2R12Energy::gbc() const { return r12eval_->gbc(); };
LinearR12::StandardApproximation MP2R12Energy::stdapp() const { return stdapprox_; };
void MP2R12Energy::set_debug(int debug) { debug_ = debug; };
int MP2R12Energy::get_debug() const { return debug_; };

double MP2R12Energy::energy()
{
  double value = emp2f12tot(AlphaAlpha) + emp2f12tot(BetaBeta) + emp2f12tot(AlphaBeta);
  return value;
}

double MP2R12Energy::emp2f12tot(SpinCase2 s) const
{
  RefSCVector unit = emp2f12_[s].clone();
  unit.assign(1.0);
  return emp2f12_[s].dot(unit);
}

double MP2R12Energy::ef12tot(SpinCase2 s) const
{
  RefSCVector unit = ef12_[s].clone();
  unit.assign(1.0);
  return ef12_[s].dot(unit);
}

static void print_psi_values(std::ostream& fout, const SCVector3& r1, const SCVector3& r2, double phi_0, double phi_1_mp2, double phi_1_r12)
{
  fout << scprintf("%9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %12.8lf %12.8lf %12.8lf",
                   r1.x(),r1.y(),r1.z(),r2.x(),r2.y(),r2.z(),phi_0,phi_1_mp2,phi_1_r12) << endl;
}

#if MP2R12ENERGY_CAN_COMPUTE_PAIRFUNCTION
double
MP2R12Energy::compute_pair_function_aa(int ij, const SCVector3& r1, const SCVector3& r2)
{
  Ref<R12Amplitudes> Amps = r12eval_->amps();
  RefSCMatrix T2aa = Amps->T2_aa();
  RefSCMatrix Rvv_aa = Amps->Rvv_aa();
  RefSCMatrix Roo_aa = Amps->Roo_aa();
  RefSCMatrix Rvo_aa = Amps->Rvo_aa();
  RefSCMatrix Rxo_aa = Amps->Rxo_aa();

  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  RefSCMatrix Caa = localkit->matrix(Caa_.rowdim(),Caa_.coldim());
  double* caa = new double[Caa_.rowdim().n()*Caa_.coldim().n()];
  Caa_.convert(caa);
  Caa.assign(caa);
  delete[] caa;
  RefSCMatrix Cvv = Caa * Rvv_aa;
  RefSCMatrix Coo = Caa * Roo_aa;
  RefSCMatrix Cov = Caa * Rvo_aa;
  RefSCMatrix Cox = Caa * Rxo_aa;

  Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
  Ref<MOIndexSpace> act_vir_space = r12info->vir_act();
  Ref<MOIndexSpace> act_occ_space = r12info->refinfo()->docc_act();
  Ref<MOIndexSpace> occ_space = r12info->refinfo()->docc();
  Ref<MOIndexSpace> ribs_space = r12info->ribs_space();

  RefSCVector phi_aa = compute_2body_values_(true,act_occ_space,act_occ_space,r1,r2);
  RefSCVector phi_vv = compute_2body_values_(true,act_vir_space,act_vir_space,r1,r2);
  RefSCVector phi_oo = compute_2body_values_(true,occ_space,occ_space,r1,r2);
  RefSCVector phi_ov = compute_2body_values_(true,occ_space,act_vir_space,r1,r2);
  RefSCVector phi_ox = compute_2body_values_(true,occ_space,ribs_space,r1,r2);

  double phi_t2 = T2aa.get_row(ij).dot(phi_vv);

  SCVector3 r12 = r1 - r2;
  const double dist12 = r12.norm();
  double phi_r12;
  phi_r12 = 0.5 * Caa.get_row(ij).dot(phi_aa) * dist12;
  phi_r12 -= 0.5 * Cvv.get_row(ij).dot(phi_vv);
  phi_r12 -= 0.5 * Coo.get_row(ij).dot(phi_oo);
  phi_r12 -= 1.0 * Cov.get_row(ij).dot(phi_ov);
  phi_r12 -= 1.0 * Cox.get_row(ij).dot(phi_ox);

  print_psi_values(ExEnv::out0(),r1,r2,phi_aa.get_element(ij),phi_t2,phi_r12);

}

void
MP2R12Energy::compute_pair_function_aa(int ij, const Ref<TwoBodyGrid>& tbgrid)
{
  Ref<R12Amplitudes> Amps = r12eval_->amps();
  RefSCMatrix T2aa = Amps->T2_aa();
  RefSCMatrix Rvv_aa = Amps->Rvv_aa();
  RefSCMatrix Roo_aa = Amps->Roo_aa();
  RefSCMatrix Rvo_aa = Amps->Rvo_aa();
  RefSCMatrix Rxo_aa = Amps->Rxo_aa();

  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  RefSCMatrix Caa = localkit->matrix(Caa_.rowdim(),Caa_.coldim());
  double* caa = new double[Caa_.rowdim().n()*Caa_.coldim().n()];
  Caa_.convert(caa);
  Caa.assign(caa);
  delete[] caa;
  RefSCMatrix Cvv = Caa * Rvv_aa;
  RefSCMatrix Coo = Caa * Roo_aa;
  RefSCMatrix Cov = Caa * Rvo_aa;
  RefSCMatrix Cox = Caa * Rxo_aa;

  Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
  Ref<MOIndexSpace> act_vir_space = r12info->vir_act();
  Ref<MOIndexSpace> act_occ_space = r12info->refinfo()->docc_act();
  Ref<MOIndexSpace> occ_space = r12info->refinfo()->docc();
  Ref<MOIndexSpace> ribs_space = r12info->ribs_space();

  const int nelem = tbgrid->nelem();
  std::stringstream output_file_name;
  output_file_name << tbgrid->name() << ".ab.pair"
                   << ij << ".txt";
  ofstream ofile(output_file_name.str().c_str());

  for(int i=0; i<nelem; i++) {
    RefSCVector phi_aa = compute_2body_values_(true,act_occ_space,act_occ_space,tbgrid->xyz1(i),tbgrid->xyz2(i));
    RefSCVector phi_vv = compute_2body_values_(true,act_vir_space,act_vir_space,tbgrid->xyz1(i),tbgrid->xyz2(i));
    RefSCVector phi_oo = compute_2body_values_(true,occ_space,occ_space,tbgrid->xyz1(i),tbgrid->xyz2(i));
    RefSCVector phi_ov = compute_2body_values_(true,occ_space,act_vir_space,tbgrid->xyz1(i),tbgrid->xyz2(i));
    RefSCVector phi_ox = compute_2body_values_(true,occ_space,ribs_space,tbgrid->xyz1(i),tbgrid->xyz2(i));

    double phi_t2 = T2aa.get_row(ij).dot(phi_vv);

    SCVector3 r12 = tbgrid->xyz1(i) - tbgrid->xyz2(i);
    const double dist12 = r12.norm();
    double phi_r12;
    phi_r12 = 0.5 * Caa.get_row(ij).dot(phi_aa) * dist12;
    phi_r12 -= 0.5 * Cvv.get_row(ij).dot(phi_vv);
    phi_r12 -= 0.5 * Coo.get_row(ij).dot(phi_oo);
    phi_r12 -= 1.0 * Cov.get_row(ij).dot(phi_ov);
    phi_r12 -= 1.0 * Cox.get_row(ij).dot(phi_ox);

    print_psi_values(ofile,tbgrid->xyz1(i),tbgrid->xyz2(i),phi_aa.get_element(ij),phi_t2,phi_r12);
  }
}

double
MP2R12Energy::compute_pair_function_ab(int ij, const SCVector3& r1, const SCVector3& r2)
{
  Ref<R12Amplitudes> Amps = r12eval_->amps();
  RefSCMatrix T2ab = Amps->T2_ab();
  RefSCMatrix Rvv_ab = Amps->Rvv_ab();
  RefSCMatrix Roo_ab = Amps->Roo_ab();
  RefSCMatrix Rvo_ab = Amps->Rvo_ab();
  RefSCMatrix Rxo_ab = Amps->Rxo_ab();

  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  RefSCMatrix Cab = localkit->matrix(Cab_.rowdim(),Cab_.coldim());
  double* cab = new double[Cab_.rowdim().n()*Cab_.coldim().n()];
  Cab_.convert(cab);
  Cab.assign(cab);
  delete[] cab;
  RefSCMatrix Cvv = Cab * Rvv_ab;
  RefSCMatrix Coo = Cab * Roo_ab;
  RefSCMatrix Cov = Cab * Rvo_ab;
  RefSCMatrix Cox = Cab * Rxo_ab;

  Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
  Ref<MOIndexSpace> act_vir_space = r12info->vir_act();
  Ref<MOIndexSpace> act_occ_space = r12info->refinfo()->docc_act();
  Ref<MOIndexSpace> occ_space = r12info->refinfo()->docc();
  Ref<MOIndexSpace> ribs_space = r12info->ribs_space();

  RefSCVector phi_aa = compute_2body_values_(false,act_occ_space,act_occ_space,r1,r2);
  RefSCVector phi_vv = compute_2body_values_(false,act_vir_space,act_vir_space,r1,r2);
  RefSCVector phi_oo = compute_2body_values_(false,occ_space,occ_space,r1,r2);
  RefSCVector phi_ov = compute_2body_values_(false,occ_space,act_vir_space,r1,r2);
  RefSCVector phi_ox = compute_2body_values_(false,occ_space,ribs_space,r1,r2);

  double phi_t2 = T2ab.get_row(ij).dot(phi_vv);

  SCVector3 r12 = r1 - r2;
  const double dist12 = r12.norm();
  double phi_r12;
  phi_r12 = 0.5*Cab.get_row(ij).dot(phi_aa) * dist12;
  phi_r12 -= 0.5 * Cvv.get_row(ij).dot(phi_vv);
  phi_r12 -= 0.5 * Coo.get_row(ij).dot(phi_oo);
  phi_r12 -= 1.0 * Cov.get_row(ij).dot(phi_ov);
  phi_r12 -= 1.0 * Cox.get_row(ij).dot(phi_ox);

  print_psi_values(ExEnv::out0(),r1,r2,phi_aa.get_element(ij),phi_t2,phi_r12);
}

void
MP2R12Energy::compute_pair_function_ab(int ij, const Ref<TwoBodyGrid>& tbgrid)
{
  Ref<R12Amplitudes> Amps = r12eval_->amps();
  RefSCMatrix T2ab = Amps->T2_ab();
  RefSCMatrix Rvv_ab = Amps->Rvv_ab();
  RefSCMatrix Roo_ab = Amps->Roo_ab();
  RefSCMatrix Rvo_ab = Amps->Rvo_ab();
  RefSCMatrix Rxo_ab = Amps->Rxo_ab();

  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  RefSCMatrix Cab = localkit->matrix(Cab_.rowdim(),Cab_.coldim());
  double* cab = new double[Cab_.rowdim().n()*Cab_.coldim().n()];
  Cab_.convert(cab);
  Cab.assign(cab);
  delete[] cab;
  RefSCMatrix Cvv = Cab * Rvv_ab;
  RefSCMatrix Coo = Cab * Roo_ab;
  RefSCMatrix Cov = Cab * Rvo_ab;
  RefSCMatrix Cox = Cab * Rxo_ab;

  Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
  Ref<MOIndexSpace> act_vir_space = r12info->vir_act();
  Ref<MOIndexSpace> act_occ_space = r12info->refinfo()->docc_act();
  Ref<MOIndexSpace> occ_space = r12info->refinfo()->docc();
  Ref<MOIndexSpace> ribs_space = r12info->ribs_space();

  const int nelem = tbgrid->nelem();
  std::stringstream output_file_name;
  output_file_name << tbgrid->name() << ".ab.pair"
                   << ij << ".txt";
  ofstream ofile(output_file_name.str().c_str());

  for(int i=0; i<nelem; i++) {
    RefSCVector phi_aa = compute_2body_values_(false,act_occ_space,act_occ_space,tbgrid->xyz1(i),tbgrid->xyz2(i));
    RefSCVector phi_vv = compute_2body_values_(false,act_vir_space,act_vir_space,tbgrid->xyz1(i),tbgrid->xyz2(i));
    RefSCVector phi_oo = compute_2body_values_(false,occ_space,occ_space,tbgrid->xyz1(i),tbgrid->xyz2(i));
    RefSCVector phi_ov = compute_2body_values_(false,occ_space,act_vir_space,tbgrid->xyz1(i),tbgrid->xyz2(i));
    RefSCVector phi_ox = compute_2body_values_(false,occ_space,ribs_space,tbgrid->xyz1(i),tbgrid->xyz2(i));
    
    double phi_t2 = T2ab.get_row(ij).dot(phi_vv);
    
    SCVector3 r12 = tbgrid->xyz1(i) - tbgrid->xyz2(i);
    const double dist12 = r12.norm();
    double phi_r12;
    phi_r12 = 0.5*Cab.get_row(ij).dot(phi_aa) * dist12;
    phi_r12 -= 0.5 * Cvv.get_row(ij).dot(phi_vv);
    phi_r12 -= 0.5 * Coo.get_row(ij).dot(phi_oo);
    phi_r12 -= 1.0 * Cov.get_row(ij).dot(phi_ov);
    phi_r12 -= 1.0 * Cox.get_row(ij).dot(phi_ox);

    print_psi_values(ofile,tbgrid->xyz1(i),tbgrid->xyz2(i),phi_aa.get_element(ij),phi_t2,phi_r12);
  }
}

RefSCVector
MP2R12Energy::compute_2body_values_(bool equiv, const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2,
                                    const SCVector3& r1, const SCVector3& r2) const
{
  const Ref<Integral> ints = r12eval_->r12info()->integral();
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
#endif

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
};

void MP2R12Energy::print_pair_energies(bool spinadapted, std::ostream& so)
{
  compute();
  
  std::string SA_str;
  switch (stdapprox_) {
    case LinearR12::StdApprox_A:
      SA_str = "A";
      break;

    case LinearR12::StdApprox_Ap:
      SA_str = "A'";
      break;

    case LinearR12::StdApprox_B:
      SA_str = "B";
      break;

    default:
      throw std::runtime_error("MP2R12Energy::print_pair_energies -- stdapprox_ is not valid");
  }

  const Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
  const double escf = r12info->refinfo()->ref()->energy();
  // WARNING assuming only RHF and ROHF
  const bool spin_polarized = r12info->refinfo()->ref()->spin_polarized();
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
      const Ref<MOIndexSpace> occ1_act = r12eval()->occ_act(case1(spincase2));
      const Ref<MOIndexSpace> occ2_act = r12eval()->occ_act(case2(spincase2));
      SpinMOPairIter ij_iter(occ1_act, occ2_act, spincase2);
      
      so << endl << indent << prepend_spincase(spincase2,"MBPT2-F12/") << SA_str << " pair energies:" << endl;
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
    
    Ref<SCMatrixKit> localkit = C_[AlphaAlpha].kit();
    RefSCVector emp2f12_0 = localkit->vector(r12eval_->dim_oo_s());
    RefSCVector emp2f12_1 = localkit->vector(r12eval_->dim_oo_t());
    RefSCVector ef12_0 = localkit->vector(r12eval_->dim_oo_s());
    RefSCVector ef12_1 = localkit->vector(r12eval_->dim_oo_t());
    
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
    const Ref<MOIndexSpace> occ_act = r12eval()->occ_act(Alpha);
    SpatialMOPairIter_eq ij_iter(occ_act);
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
    
    so << endl << indent << "Singlet MBPT2-F12/" << SA_str << " pair energies:" << endl;
    so << indent << scprintf("    i       j        mp2(ij)        r12(ij)      mp2-r12(ij)") << endl;
    so << indent << scprintf("  -----   -----   ------------   ------------   ------------") << endl;
    const int nocc_act = occ_act->rank();
    for(int i=0,ij=0;i<nocc_act;i++) {
      for(int j=0;j<=i;j++,ij++) {
        const double ep_f12_0 = ef12_0.get_element(ij);
        const double ep_mp2f12_0 = emp2f12_0.get_element(ij);
        so << indent << scprintf("  %3d     %3d     %12.9lf   %12.9lf   %12.9lf",
                                 i+1,j+1,ep_mp2f12_0-ep_f12_0,ep_f12_0,ep_mp2f12_0) << endl;
      }
    }
    
    so << endl << indent << "Triplet MBPT2-F12/" << SA_str << " pair energies:" << endl;
    so << indent << scprintf("    i       j        mp2(ij)        r12(ij)      mp2-r12(ij)") << endl;
    so << indent << scprintf("  -----   -----   ------------   ------------   ------------") << endl;
    for(int i=0,ij=0;i<nocc_act;i++) {
      for(int j=0;j<i;j++,ij++) {
        const double ep_f12_1 = ef12_1.get_element(ij);
        const double ep_mp2f12_1 = emp2f12_1.get_element(ij);
        so << indent << scprintf("  %3d     %3d     %12.9lf   %12.9lf   %12.9lf",
                                 i+1,j+1,ep_mp2f12_1-ep_f12_1,ep_f12_1,ep_mp2f12_1) << endl;
      }
    }
    
  }
  
  const double ef12_corr_energy = ef12tot(AlphaAlpha) + ef12tot(BetaBeta) + ef12tot(AlphaBeta);
  const double emp2f12_corr_energy = emp2f12tot(AlphaAlpha) + emp2f12tot(BetaBeta) + emp2f12tot(AlphaBeta);
  
  ///////////////////////////////////////////////////////////////
  // The computation of the MP2 energy is now complete on each
  // node;
  ///////////////////////////////////////////////////////////////
  
  if (spinadapted) {
    so <<endl<<indent
    <<scprintf("Singlet MP2 correlation energy [au]:          %17.12lf\n", emp2f12tot_0 - ef12tot_0);
    so <<indent
    <<scprintf("Triplet MP2 correlation energy [au]:          %17.12lf\n", emp2f12tot_1 - ef12tot_1);
    so <<indent
    <<scprintf("Singlet (MP2)-F12/%2s correlation energy [au]: %17.12lf\n", SA_str.c_str(), ef12tot_0);
    so <<indent
    <<scprintf("Triplet (MP2)-F12/%2s correlation energy [au]: %17.12lf\n", SA_str.c_str(), ef12tot_1);
    so <<indent
    <<scprintf("Singlet MP2-F12/%2s correlation energy [au]:   %17.12lf\n", SA_str.c_str(),
    emp2f12tot_0);
    so <<indent
    <<scprintf("Triplet MP2-F12/%2s correlation energy [au]:   %17.12lf\n", SA_str.c_str(),
    emp2f12tot_1);
  }
  
  double etotal = escf + emp2f12_corr_energy;
  so <<endl<<indent
  <<scprintf("RHF energy [au]:                           %17.12lf\n", escf);
  so <<indent
  <<scprintf("MP2 correlation energy [au]:               %17.12lf\n", emp2f12_corr_energy - ef12_corr_energy);
  so <<indent
  <<scprintf("(MBPT2)-F12/%2s correlation energy [au]:    %17.12lf\n", SA_str.c_str(), ef12_corr_energy);
  so <<indent
  <<scprintf("MBPT2-F12/%2s correlation energy [au]:      %17.12lf\n", SA_str.c_str(),
  emp2f12_corr_energy);
  so <<indent
  <<scprintf("MBPT2-F12/%2s energy [au]:                  %17.12lf\n", SA_str.c_str(), etotal) << endl;
  
  so.flush();
  
  return;
}


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
