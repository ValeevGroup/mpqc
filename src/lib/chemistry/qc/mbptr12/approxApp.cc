//
// approxApp.cc
//
// Copyright (C) 2006 Edward Valeev
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

#include <stdexcept>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <mpqc_config.h>
#include <util/misc/formio.h>
#include <util/misc/regtime.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/state_text.h>
#include <util/state/state_bin.h>
#include <math/scmat/local.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <math/distarray4/distarray4.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <math/mmisc/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/creator.h>
#include <chemistry/qc/mbptr12/container.h>
#include <chemistry/qc/mbptr12/compute_tbint_tensor.h>
#include <chemistry/qc/mbptr12/twoparticlecontraction.h>
#include <chemistry/qc/lcao/utils.h>
#include <chemistry/qc/lcao/utils.impl.h>
#include <util/misc/print.h>

using namespace std;
using namespace sc;

#define INCLUDE_Q 1

void
R12IntEval::compute_BApp_()
{
  if (evaluated_)
    return;

  const bool vbs_eq_obs = r12world()->obs_eq_vbs();
  const bool abs_eq_obs = r12world()->obs_eq_ribs();
  const unsigned int maxnabs = r12world()->r12tech()->maxnabs();

  Timer tim_B_app_App("B(app. A'') intermediate");
  ExEnv::out0() << endl << indent
  << "Entered B(app. A'') intermediate evaluator" << endl;
  ExEnv::out0() << incindent;

  //
  // if I'm relativistic, all relativistic terms will either be treated as K(i.e. dropped from F giving h+J)
  // OR kept if user asks for that. The latter assumes that all higher-order terms in
  // DK hamiltonian commute with Q12 f12
  //
  // to accomodate the former option (drop reltivistic terms from h+J) do this
  // 1) make new space hJnr = h(nr) + J = (h(r) + J) - (H(r) - H(nr)) = hJ - dH
  //
  Ref<OrbitalSpace> hJnr[NSpinCases1];
  if (this->dk() > 0) {

    const R12Technology::H0_dk_approx_pauli H0_dk_approx_pauli = r12world()->r12tech()->H0_dk_approx();
    const bool H0_dk_keep = r12world()->r12tech()->H0_dk_keep();

    const int nspins1 = this->nspincases1();
    for (int s = 0; s < nspins1; ++s) {
      const SpinCase1 spin = static_cast<SpinCase1> (s);

      // is this supported?
      if (!vbs_eq_obs && maxnabs < 2)
        throw FeatureNotImplemented("Relativistic R12/A'' computations are not supported with maxnabs==1 && VBS!=OBS",__FILE__,__LINE__);
      if (abs_eq_obs)
        throw FeatureNotImplemented("Relativistic R12/A'' computations are not supported with ABS==OBS",__FILE__,__LINE__);

      // which space is used as RIBS?
      Ref<OrbitalSpace> ribs = (maxnabs < 2) ? this->orbs(spin) : r12world()->ribs_space();
      // get AO space for RIBS
      const Ref<OrbitalSpace>& aoribs =
        this->ao_registry()->value(ribs->basis());

      Ref<OrbitalSpace> hJ_x_P = (maxnabs < 2) ? hj_x_p(spin) : hj_x_P(spin);
      if (H0_dk_approx_pauli == R12Technology::H0_dk_approx_pauli_false && !H0_dk_keep) { // use nonrelativistic hamiltonian in h+J
        const Ref<OrbitalSpace>& x = GGspace(spin);
        const Ref<OrbitalSpace>& aox = this->ao_registry()->value(x->basis());
        // compute dH = H(rel) - H(nonrel) in AO basis
        RefSCMatrix Hr = this->fock(aox, aoribs, spin, 0.0, 0.0, 1.0);
        const std::string nonrel_hkey =
          ParsedOneBodyIntKey::key(aox->id(), aoribs->id(),
                                   std::string("H"));
        RefSCMatrix Hnr = fockbuild_runtime()->get(nonrel_hkey);
        RefSCMatrix dH = Hnr.clone();
        dH->convert(Hr);
        Hnr.scale(-1.0);
        dH.accumulate(Hnr);
        Hr = 0;
        Hnr = 0;
        // transform to MO basis, transpose so that x dimension is coldim
        RefSCMatrix dH_mo = x->coefs().t() * dH * ribs->coefs();
        dH_mo = dH_mo.t();

        // Make new space hJ - dH
        std::string id = x->id();
        id += "_hJnr(";
        id += ribs->id();
        id += ")";
        ExEnv::out0() << "id = " << id << endl;
        std::string name = "(hJnr)-weighted space";
#if 1
        hJnr[spin] = new OrbitalSpace(id, name, hJ_x_P, hJ_x_P->coefs() - ribs->coefs() * dH_mo, hJ_x_P->basis());
#else
        hJnr[spin] = new OrbitalSpace(id, name, hJ_x_P, hJ_x_P->coefs(), hJ_x_P->basis());
#endif

        this->r12world()->world()->tfactory()->orbital_registry()->add(make_keyspace_pair(hJnr[s]));
      }
      else { // use pure h+J
        hJnr[spin] = hJ_x_P;
      }

    }
    if (hJnr[Beta] == 0) hJnr[Beta] = hJnr[Alpha];
  }

  for(int s=0; s<nspincases2(); s++) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);
    Ref<OrbitalSpace> GGspace1 = GGspace(spin1);
    Ref<OrbitalSpace> GGspace2 = GGspace(spin2);
    const bool x1_eq_x2 = (GGspace1 == GGspace2);

    if (dim_oo(spincase2).n() == 0)
      continue;

#if INCLUDE_Q

    std::string Qlabel = prepend_spincase(spincase2,"Q(A'') intermediate");
    Timer tim_Q(Qlabel);
    ExEnv::out0() << endl << indent
                  << "Entered " << Qlabel << " evaluator" << endl;
    ExEnv::out0() << incindent;

    // compute Q = X_{xy}^{xy_{hj}}
    RefSCMatrix Q;
    if (maxnabs > 1) { // if can only use 2 RI index, h+J can be resolved by the RIBS
      Ref<OrbitalSpace> hj_x2 = hj_x_P(spin2);
      // if I'm relativistic, use the hJ adjusted for user options
      if (this->dk() > 0) { hj_x2 = hJnr[spin2]; }
      compute_X_(Q,spincase2,GGspace1,GGspace2,
                 GGspace1,hj_x2);
    }
    else { // else do RI in orbital basis...
      if (vbs_eq_obs) { // which is just p if VBS == OBS
        Ref<OrbitalSpace> hj_x2 = hj_x_p(spin2);
        // if I'm relativistic, use the hJ adjusted for user options
        if (this->dk() > 0) { hj_x2 = hJnr[spin2]; }
        compute_X_(Q,spincase2,GGspace1,GGspace2,
                   GGspace1,hj_x2);
      }
      else { // if VBS != OBS, p = m + a
        Ref<OrbitalSpace> hj_x2 = hj_x_m(spin2);
        compute_X_(Q,spincase2,GGspace1,GGspace2,
                   GGspace1,hj_x2);
	    hj_x2 = hj_x_a(spin2);
	    compute_X_(Q,spincase2,GGspace1,GGspace2,
	               GGspace1,hj_x2);
      }
    }

    if (x1_eq_x2) {
      if (maxnabs > 1) { // if can only use 2 RI index, h+J can be resolved by the RIBS
        Ref<OrbitalSpace> hj_x1 = hj_x_P(spin1);
        // if I'm relativistic, use the hJ adjusted for user options
        if (this->dk() > 0) { hj_x1 = hJnr[spin1]; }
        compute_X_(Q,spincase2,GGspace1,GGspace2,
                   hj_x1,GGspace2);
      }
      else { // else do RI in orbital basis...
	    if (vbs_eq_obs) { // which is just p if VBS == OBS
	      Ref<OrbitalSpace> hj_x1 = hj_x_p(spin1);
	      // if I'm relativistic, use the hJ adjusted for user options
	      if (this->dk() > 0) { hj_x1 = hJnr[spin1]; }
	      compute_X_(Q,spincase2,GGspace1,GGspace2,
	                 hj_x1,GGspace2);
	    }
	    else { // if VBS != OBS, p = m + a
	      Ref<OrbitalSpace> hj_x1 = hj_x_m(spin1);
	      compute_X_(Q,spincase2,GGspace1,GGspace2,
	                 hj_x1,GGspace2);
	      hj_x1 = hj_x_a(spin1);
	      compute_X_(Q,spincase2,GGspace1,GGspace2,
	                 hj_x1,GGspace2);
	    }
      }
    }
    else {
      Q.scale(2.0);
      if (spincase2 == AlphaBeta) {
        symmetrize<false>(Q,Q,GGspace1,GGspace1);
      }
    }

    ExEnv::out0() << decindent;
    ExEnv::out0() << indent << "Exited " << Qlabel << " evaluator" << endl;
    tim_Q.exit();

    if (debug_ >= DefaultPrintThresholds::mostO4) {
      ExEnv::out0() << indent << __FILE__ << ": "<<__LINE__<<"\n";
      std::string label = prepend_spincase(spincase2,"Q(A'') contribution");
      Q.print(label.c_str());
    }
    B_[s].accumulate(Q); Q = 0;
#endif // INCLUDE_Q

    // Bra-Ket symmetrize the B(A'') contribution
    B_[s].scale(0.5);
    RefSCMatrix B_t = B_[s].t();
    B_[s].accumulate(B_t);
  }

  globally_sum_intermeds_();

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited B(app. A'') intermediate evaluator" << endl;

  tim_B_app_App.exit();
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
