//
// compute_DC_energy_GenRefansatz2.cc
//
// Copyright (C) 2009 Martin Torheyden
//
// Author: Martin Torheyden <mtorhey@vt.edu>
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

#include <stdexcept>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <scconfig.h>
#include <util/misc/formio.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/state_text.h>
#include <util/state/state_bin.h>
#include <math/scmat/local.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/mbptr12/blas.h>
#include <chemistry/qc/mbptr12/distarray4.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/creator.h>
#include <chemistry/qc/mbptr12/container.h>
#include <chemistry/qc/mbptr12/compute_tbint_tensor.h>
#include <chemistry/qc/mbptr12/contract_tbint_tensor.h>
#include <chemistry/qc/mbptr12/contract_tbint_tensors_to_obtensor.h>
#include <chemistry/qc/mbptr12/twoparticlecontraction.h>
#include <chemistry/qc/mbptr12/utils.h>
#include <chemistry/qc/mbptr12/utils.impl.h>
#include <chemistry/qc/mbptr12/print.h>
#include <chemistry/qc/psi/psiwfn.h>

using namespace std;
using namespace sc;

double PsiCorrWavefunction_PT2R12::compute_DC_energy_GenRefansatz2() {
  using LinearR12::TwoParticleContraction;
  using LinearR12::Direct_Contraction;
  Ref<R12IntEval> thisref(r12eval_);
  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;

  // initialize intermediates I and J
  RefSCMatrix I_alpha_q[NSpinCases1];
  RefSCMatrix J_alpha_q[NSpinCases1];
  for(int s1=0; s1<NSpinCases1; s1++) {
    SpinCase1 singlespin = static_cast<SpinCase1>(s1);
    Ref<OrbitalSpace> cabs_sgsp = r12eval_->r12world()->cabs_space(singlespin);
    Ref<OrbitalSpace> orbs_sgsp = r12eval_->orbs(singlespin);
    RefSCDimension cabs_sgsp_dim = cabs_sgsp->dim();
    RefSCDimension orbs_sgsp_dim = orbs_sgsp->dim();
    unsigned int ncabs_sgsp = cabs_sgsp_dim.n();
    unsigned int norbs_sgsp = orbs_sgsp_dim.n();
    I_alpha_q[singlespin] = r12eval_->B(AlphaBeta).kit()->matrix(cabs_sgsp_dim,orbs_sgsp_dim);
    J_alpha_q[singlespin] = r12eval_->B(AlphaBeta).kit()->matrix(cabs_sgsp_dim,orbs_sgsp_dim);
    I_alpha_q[singlespin].assign(0.0);
    J_alpha_q[singlespin].assign(0.0);
  }

  // compute intermediates I and J
  Timer IJtimer("computation of intermediates I and J for the D contribution");
  ExEnv::out0() << "Starting computation of intermediates I and J for the D contribution." << endl;
  for(int s2=0; s2<NSpinCases2; s2++) {
    SpinCase2 pairspin = static_cast<SpinCase2>(s2);
    SpinCase1 spin1 = case1(pairspin);
    SpinCase1 spin2 = case2(pairspin);

    Ref<R12WavefunctionWorld> r12world = r12eval_->r12world();
    Ref<R12RefWavefunction> ref = r12world->ref();
    Ref<OrbitalSpace> gg1space = r12eval_->ggspace(spin1);
    Ref<OrbitalSpace> gg2space = r12eval_->ggspace(spin2);
    Ref<OrbitalSpace> GG1space = r12eval_->GGspace(spin1);
    Ref<OrbitalSpace> GG2space = r12eval_->GGspace(spin2);
    Ref<OrbitalSpace> orbs1 = r12eval_->orbs(spin1);
    Ref<OrbitalSpace> orbs2 = r12eval_->orbs(spin2);
    Ref<OrbitalSpace> cabs1 = r12eval_->r12world()->cabs_space(spin1);
    Ref<OrbitalSpace> cabs2 = r12eval_->r12world()->cabs_space(spin2);

    // computing I contribution
    {
      Ref<OrbitalSpace> gammaFgamma_p_p1 = r12eval_->gammaFgamma_p_p(spin1);
      Ref<TwoParticleContraction> dircontract_p_p = new Direct_Contraction(gg2space->rank(),gg1space->rank(),1.0);
      std::vector<std::string> tforms_bra_p_p;
      {
        R12TwoBodyIntKeyCreator tform_creator(
                            r12world->world()->moints_runtime4(),
                            cabs2,GG2space,gammaFgamma_p_p1,GG1space,
                            r12world->r12tech()->corrfactor(),true
                            );
        fill_container(tform_creator,tforms_bra_p_p);
      }
      std::vector<std::string> tforms_ket_p_p;
      {
        R12TwoBodyIntKeyCreator tform_creator(
                            r12world->world()->moints_runtime4(),
                            gg2space,GG2space,gg1space,GG1space,
                            r12world->r12tech()->corrfactor()
                            );
        fill_container(tform_creator,tforms_ket_p_p);
      }
      TwoBodyTensorInfo tbtensor_info_bra(r12eval_->corrfactor()->tbint_type_f12());
      TwoBodyTensorInfo tbtensor_info_ket(TwoBodyTensorInfo::geminalcoeff);
      r12eval_->contract_tbint_tensors_to_obtensor<ManyBodyTensors::I_to_T,
                                         ManyBodyTensors::I_to_T,
                                         true,false,false>(
                                         I_alpha_q[spin2],
                                         pairspin,
                                         tbtensor_info_bra,tbtensor_info_ket,
                                         cabs2,gammaFgamma_p_p1,
                                         GG2space,GG1space,
                                         gg2space,gg1space,
                                         GG2space,GG1space,
                                         dircontract_p_p,tforms_bra_p_p,tforms_ket_p_p);
    }
    if(pairspin==AlphaBeta) {
      Ref<OrbitalSpace> gammaFgamma_p_p2 = r12eval_->gammaFgamma_p_p(spin2);
      Ref<TwoParticleContraction> dircontract_p_p = new Direct_Contraction(gg1space->rank(),gg2space->rank(),1.0);
      std::vector<std::string> tforms_bra_p_p;
      {
        R12TwoBodyIntKeyCreator tform_creator(
                            r12world->world()->moints_runtime4(),
                            cabs1,GG1space,gammaFgamma_p_p2,GG2space,
                            r12world->r12tech()->corrfactor(),true
                            );
        fill_container(tform_creator,tforms_bra_p_p);
      }
      std::vector<std::string> tforms_ket_p_p;
      {
        R12TwoBodyIntKeyCreator tform_creator(
                            r12world->world()->moints_runtime4(),
                            gg1space,GG1space,gg2space,GG2space,
                            r12world->r12tech()->corrfactor()
                            );
        fill_container(tform_creator,tforms_ket_p_p);
      }
      TwoBodyTensorInfo tbtensor_info_bra(r12eval_->corrfactor()->tbint_type_f12());
      TwoBodyTensorInfo tbtensor_info_ket(TwoBodyTensorInfo::geminalcoeff);
      r12eval_->contract_tbint_tensors_to_obtensor<ManyBodyTensors::I_to_T,
                                         ManyBodyTensors::I_to_T,
                                         true,false,false>(
                                         I_alpha_q[spin1],
                                         pairspin,
                                         tbtensor_info_bra,tbtensor_info_ket,
                                         cabs1,gammaFgamma_p_p2,
                                         GG1space,GG2space,
                                         gg1space,gg2space,
                                         GG1space,GG2space,
                                         dircontract_p_p,tforms_bra_p_p,tforms_ket_p_p);
    }

    // computing J contribution
    {
      Ref<OrbitalSpace> Fgamma_p_P1 = r12eval_->Fgamma_p_P(spin1);
      Ref<TwoParticleContraction> dircontract_p_p = new Direct_Contraction(gg2space->rank(),gg1space->rank(),1.0);
      std::vector<std::string> tforms_bra_p_p;
      {
        R12TwoBodyIntKeyCreator tform_creator(
                            r12world->world()->moints_runtime4(),
                            cabs2,GG2space,Fgamma_p_P1,GG1space,
                            r12world->r12tech()->corrfactor(),true
                            );
        fill_container(tform_creator,tforms_bra_p_p);
      }
      std::vector<std::string> tforms_ket_p_p;
      {
        R12TwoBodyIntKeyCreator tform_creator(
                            r12world->world()->moints_runtime4(),
                            gg2space,GG2space,gg1space,GG1space,
                            r12world->r12tech()->corrfactor()
                            );
        fill_container(tform_creator,tforms_ket_p_p);
      }
      TwoBodyTensorInfo tbtensor_info_bra(r12eval_->corrfactor()->tbint_type_f12());
      TwoBodyTensorInfo tbtensor_info_ket(TwoBodyTensorInfo::geminalcoeff);
      r12eval_->contract_tbint_tensors_to_obtensor<ManyBodyTensors::I_to_T,
                                         ManyBodyTensors::I_to_T,
                                         true,false,false>(
                                         J_alpha_q[spin2],
                                         pairspin,
                                         tbtensor_info_bra,tbtensor_info_ket,
                                         cabs2,Fgamma_p_P1,
                                         GG2space,GG1space,
                                         gg2space,gg1space,
                                         GG2space,GG1space,
                                         dircontract_p_p,tforms_bra_p_p,tforms_ket_p_p);
    }
    if(pairspin==AlphaBeta) {
      Ref<OrbitalSpace> Fgamma_p_P2 = r12eval_->Fgamma_p_P(spin2);
      Ref<TwoParticleContraction> dircontract_p_p = new Direct_Contraction(gg1space->rank(),gg2space->rank(),1.0);
      std::vector<std::string> tforms_bra_p_p;
      {
        R12TwoBodyIntKeyCreator tform_creator(
                            r12world->world()->moints_runtime4(),
                            cabs1,GG1space,Fgamma_p_P2,GG2space,
                            r12world->r12tech()->corrfactor(),true
                            );
        fill_container(tform_creator,tforms_bra_p_p);
      }
      std::vector<std::string> tforms_ket_p_p;
      {
        R12TwoBodyIntKeyCreator tform_creator(
                            r12world->world()->moints_runtime4(),
                            gg1space,GG1space,gg2space,GG2space,
                            r12world->r12tech()->corrfactor()
                            );
        fill_container(tform_creator,tforms_ket_p_p);
      }
      TwoBodyTensorInfo tbtensor_info_bra(r12eval_->corrfactor()->tbint_type_f12());
      TwoBodyTensorInfo tbtensor_info_ket(TwoBodyTensorInfo::geminalcoeff);
      r12eval_->contract_tbint_tensors_to_obtensor<ManyBodyTensors::I_to_T,
                                         ManyBodyTensors::I_to_T,
                                         true,false,false>(
                                         J_alpha_q[spin1],
                                         pairspin,
                                         tbtensor_info_bra,tbtensor_info_ket,
                                         cabs1,Fgamma_p_P2,
                                         GG1space,GG2space,
                                         gg1space,gg2space,
                                         GG1space,GG2space,
                                         dircontract_p_p,tforms_bra_p_p,tforms_ket_p_p);
    }
  }
  ExEnv::out0() << "Computation of intermediates I and J for the D contribution finished." << endl;
  IJtimer.exit();

  // Add up intermediates I and J
  RefSCMatrix IminusJ[NSpinCases1];
  for(int s1=0; s1<NSpinCases1; s1++) {
    SpinCase1 singlespin = static_cast<SpinCase1>(s1);
    I_alpha_q[singlespin].print(prepend_spincase(singlespin,"compute_DC_energy_GenRefansatz2 I_alpha_q").c_str());
    J_alpha_q[singlespin].print(prepend_spincase(singlespin,"compute_DC_energy_GenRefansatz2 J_alpha_q").c_str());
    IminusJ[singlespin] = I_alpha_q[singlespin] - J_alpha_q[singlespin];
    //IminusJ[singlespin].print(prepend_spincase(singlespin,"IminusJ").c_str());
  }

  double energy_pairspin[NSpinCases2];
  double energy = 0.0;
  // Evaluate r tensor
  for(int s2=0; s2<NSpinCases2; s2++) {
    RefSCMatrix intermed1;    /// lambda * R
    SpinCase2 pairspin = static_cast<SpinCase2>(s2);
    SpinCase1 spin1 = case1(pairspin);
    SpinCase1 spin2 = case2(pairspin);
    const bool antisymmetrize = pairspin!=AlphaBeta;

    Ref<R12WavefunctionWorld> r12world = r12eval_->r12world();
    Ref<R12RefWavefunction> ref = r12world->ref();
    Ref<OrbitalSpace> gg1space = r12eval_->ggspace(spin1);
    Ref<OrbitalSpace> gg2space = r12eval_->ggspace(spin2);
    Ref<OrbitalSpace> GG1space = r12eval_->GGspace(spin1);
    Ref<OrbitalSpace> GG2space = r12eval_->GGspace(spin2);
    Ref<OrbitalSpace> orbs1 = r12eval_->orbs(spin1);
    Ref<OrbitalSpace> orbs2 = r12eval_->orbs(spin2);
    Ref<OrbitalSpace> cabs1 = r12eval_->r12world()->cabs_space(spin1);
    Ref<OrbitalSpace> cabs2 = r12eval_->r12world()->cabs_space(spin2);
    Ref<OrbitalSpace> IminusJ_p_A1 = r12eval_->obtensor_p_A(IminusJ[spin1],spin1);
    Ref<OrbitalSpace> IminusJ_p_A2 = r12eval_->obtensor_p_A(IminusJ[spin2],spin2);

    // choose pairspin=AlphaBeta, since the indices are not equivalent.
    SpinMOPairIter rbra(GG1space,GG2space,pairspin);
    //SpinMOPairIter rket(orbs1,IminusJ_p_A2,pairspin);
    SpinMOPairIter rket(gg1space,gg2space,pairspin);
    RefSCDimension bradim = new SCDimension(rbra.nij());
    RefSCDimension ketdim = new SCDimension(rket.nij());
    RefSCMatrix R = localkit->matrix(bradim,ketdim);
    R.assign(0.0);

    {
      std::vector<std::string> tforms;
      {
        R12TwoBodyIntKeyCreator tform_creator(
                            r12world->world()->moints_runtime4(),
                            GG1space,orbs1,GG2space,IminusJ_p_A2,
                            r12world->r12tech()->corrfactor(),true
                            );
        fill_container(tform_creator,tforms);
      }
      r12eval_->compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
                        R,r12eval_->corrfactor()->tbint_type_f12(),
                        GG1space,orbs1,
                        GG2space,IminusJ_p_A2,true,tforms);
    }

    RefSymmSCMatrix lambda = lambda_refmo(pairspin);
    intermed1 = R*lambda;
    intermed1.scale(2.0);
    lambda = 0;
    R = 0;
    RefSCMatrix T = C(pairspin);

    RefSCMatrix intermed2 = T*intermed1;
    T = 0;
    intermed1 = 0;

    ExEnv::out0() << "projected contribution to pair energies:" << endl;
    SpinMOPairIter gg_iter(gg1space,gg2space,pairspin);
    for(gg_iter.start(); int(gg_iter); gg_iter.next()) {
      int i = gg_iter.i();
      int j = gg_iter.j();
      int ij = gg_iter.ij();
      ExEnv::out0() << setw(6) << i << setw(6) << j
                    << setw(20) << setprecision(12) << intermed2.get_element(ij,ij) << endl;
    }

    energy_pairspin[pairspin] = 0.0;
    if(nfzc_==0) {
      ExEnv::out0() << "Computing projected contribution to energy as trace of intermediate matrix." << endl;
      energy_pairspin[pairspin] = intermed2.trace();
    }
    else {
      ExEnv::out0() << "Computing projected contribution to energy as sum of pair energies." << endl;
      for(gg_iter.start(); int(gg_iter); gg_iter.next()) {
        int i = gg_iter.i();
        int j = gg_iter.j();
        int ij = gg_iter.ij();
        if((i>=nfzc_) && (j>=nfzc_)) {
          energy_pairspin[pairspin]+=intermed2.get_element(ij,ij);
        }
      }
    }
    energy += energy_pairspin[pairspin];
    ExEnv::out0() << "projected contribution to energy " << setprecision(12) << energy_pairspin[pairspin] << endl;
  }

  return(energy);
}
