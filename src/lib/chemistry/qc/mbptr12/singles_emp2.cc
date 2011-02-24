//
// singles_emp2.cc
//
// Copyright (C) 2007 Edward Valeev
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
#include <math.h>
#include <limits.h>

#include <scconfig.h>
#include <util/misc/formio.h>
#include <util/misc/regtime.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/state_text.h>
#include <util/state/state_bin.h>
#include <math/scmat/matrix.h>
#include <math/scmat/blocked.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbptr12/distarray4.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>

using namespace std;
using namespace sc;

double
R12IntEval::compute_emp2_obs_singles(bool obs_singles)
{
  Ref<MessageGrp> msg = r12world()->world()->msg();
  Ref<MemoryGrp> mem = r12world()->world()->mem();
  Ref<ThreadGrp> thr = r12world()->world()->thr();

  Timer tim("OBS singles MP2 energy");
  std::string evalname;
  {
    std::ostringstream oss;
    oss << (obs_singles ? "OBS" : "VBS") << " singles MP2 energy evaluator";
    evalname = oss.str();
  }
  ExEnv::out0() << endl << indent << "Entered " << evalname << endl;
  ExEnv::out0() << incindent;

  int me = msg->me();
  int nproc = msg->n();

  double result = 0.0;
  for(int s=0; s<nspincases1(); s++) {
    const SpinCase1 spin = static_cast<SpinCase1>(s);

    Ref<OrbitalSpace> occ_act = this->occ_act(spin);
    Ref<OrbitalSpace> vir_act = this->vir_act(spin);
    RefSCMatrix Fia = this->fock(occ_act,vir_act,spin);

#define DEBUG_EMP2_SINGLES 0
#if DEBUG_EMP2_SINGLES
    Fia.print(prepend_spincase(spin,"occ/vir Fock matrix").c_str());

    Ref<OrbitalSpace> orbs = orbs(spin);
    RefSCMatrix Fpp = fock(orbs,orbs,spin);
    orbs->print_detail(ExEnv::out0());
    Fpp.print(prepend_spincase(spin,"OBS/OBS Fock matrix").c_str());
#endif

    const int ni = occ_act->rank();
    const int na = vir_act->rank();
    const RefDiagSCMatrix& ievals = occ_act->evals();
    const RefDiagSCMatrix& aevals = vir_act->evals();

    for(int i=0; i<ni; i++) {
      for(int a=0; a<na; a++) {
        const double fia = Fia(i,a);
        result -= fia*fia/(-ievals(i)+aevals(a));
      }
    }

    if (! spin_polarized()) result *= 2.0;
  }

  ExEnv::out0() << indent << "E(MP2 OBS singles) = " << scprintf("%25.15lf",result) << endl;
  ExEnv::out0() << decindent << indent << "Exited " << evalname << endl;

  return result;
}

double
R12IntEval::compute_emp2_cabs_singles()
{
  Ref<MessageGrp> msg = r12world()->world()->msg();
  Ref<MemoryGrp> mem = r12world()->world()->mem();
  Ref<ThreadGrp> thr = r12world()->world()->thr();

  Timer tim("CABS singles MP2 energy");
  std::string evalname("CABS singles MP2 energy evaluator");
  ExEnv::out0() << endl << indent << "Entered " << evalname << endl;
  ExEnv::out0() << incindent;

  double result = 0.0;
  for(int s=0; s<nspincases1(); s++) {
    const SpinCase1 spin = static_cast<SpinCase1>(s);

    Ref<OrbitalSpace> occ = this->occ(spin);
    if (occ->rank() == 0) continue;
    Ref<OrbitalSpace> vir = this->vir(spin);
#if 0
    Ref<OrbitalSpace> cabs = cabs_space_canonical(spin);
#endif
    Ref<OrbitalSpace> cabs = r12world()->cabs_space(spin);
    if (cabs->rank() == 0) continue;
    RefSCMatrix FiA = fock(occ,cabs,spin);
    RefSCMatrix Fia = fock(occ,vir,spin);
    RefSCMatrix FaA = fock(vir,cabs,spin);
    RefSCMatrix FAA = fock(cabs,cabs,spin);

    const int ni = occ->rank();
    const int na = vir->rank();
    const int nA = cabs->rank();
    const RefDiagSCMatrix& ievals = occ->evals();
    const RefDiagSCMatrix& aevals = vir->evals();

    double* FiA_ptr = new double[ni * nA];
    FiA.convert(FiA_ptr);
    double* Fia_ptr = new double[ni * na];
    Fia.convert(Fia_ptr);

    for(int i=0; i<ni; i++) {
      const double ieval = ievals(i);

      // zeroth-order Hamiltonian in CABS basis includes the standard contribution
      // H0(i)_AB = F_AB - \delta_AB F_ii
      RefSCMatrix H0_i = FAA.copy();
      for(int a=0; a<nA; a++) {
        H0_i.accumulate_element(a, a, -ieval);
      }
      // H1(i)_A = F_iA
      RefSCVector H1_i = H0_i.kit()->vector(H0_i.coldim());  H1_i.assign(FiA_ptr + i*nA);

#define INCLUDE_EBC_BLOCK_ZEROTH_ORDER 1
#if INCLUDE_EBC_BLOCK_ZEROTH_ORDER
      // if including EBC block of the Fock matrix (F_aA), then it's easy to block diagonalize it
      // so that I only need to invert F_AA
      // this is done in the same manner as the coupling is handled in MP2-R12
      // H0(i)_AB -= F_Aa * (H0(i)_ab)^{-1} * F_bB, where H0(i)_ab is diagonal for canonicalized references
      // H1(i)_A -= F_Ac * (H0(i)_ab)^{-1} * H1(i)_b   <--- this term is nonzero only for non-Brillouin references
      RefDiagSCMatrix H0_i_inv_vv = FaA.kit()->diagmatrix(FaA.rowdim());
      for(int c=0; c<na; c++) {
        H0_i_inv_vv.set_element(c, 1.0/(aevals(c) - ieval));
      }

      // and the contribution from the EBC block to H0
      RefSCMatrix H1_H0_inv = FaA.t() * H0_i_inv_vv;  H1_H0_inv.scale(-1.0);
      H0_i.accumulate(H1_H0_inv * FaA);
      // ... and to H1
      RefSCVector fia = H1_H0_inv.kit()->vector(H1_H0_inv.coldim());  fia.assign(Fia_ptr + i*na);
      H1_i.accumulate( H1_H0_inv * fia );
#endif

      RefSCVector C1_i = H0_i.gi() * H1_i; C1_i.scale(-1.0);
      result += H1_i.dot(C1_i);

    }

    if (! spin_polarized()) result *= 2.0;

    delete[] FiA_ptr;
  }

  ExEnv::out0() << indent << "E(MP2 CABS singles) = " << scprintf("%25.15lf",result) << endl;
  ExEnv::out0() << decindent << indent << "Exited " << evalname << endl;

  return result;
}

// rotates CABS so that the Fock matrix is diagonal
const Ref<OrbitalSpace>&
R12IntEval::cabs_space_canonical(SpinCase1 spin)
{
  static Ref<OrbitalSpace> cabs_canonical[] = {0, 0};

  assert(r12world()->obs_eq_ribs() == false);

  if (cabs_canonical[spin] != 0)
    return cabs_canonical[spin];

  Ref<OrbitalSpace> cabs = r12world()->cabs_space(spin);
  const int ncabs = cabs->rank();

  // note that I'm overriding pauli flag here -- true Fock matrix must always be used
  const double scale_J = 1.0;
  const double scale_K = 1.0;
  const double scale_H = 1.0;
  const int pauli = 0;
  RefSCMatrix F_cabs = fock(cabs,cabs,spin,scale_J,scale_K,scale_H,pauli);
  RefSymmSCMatrix F_cabs_lt(F_cabs.rowdim(),F_cabs->kit());
  F_cabs_lt.assign(0.0);
  // CABS is a blocked space, convert to a symmetric matrix block by block
  {
    BlockedSCMatrix* bF = require_dynamic_cast<BlockedSCMatrix*>(F_cabs.pointer(),"expected a blocked matrix");
    BlockedSymmSCMatrix* bF_lt = require_dynamic_cast<BlockedSymmSCMatrix*>(F_cabs_lt.pointer(),"expected a blocked matrix");
    Ref<SCBlockInfo> blocks = F_cabs.rowdim()->blocks();
    const int nb = blocks->nblock();
    for(int b=0; b<nb; ++b) {
      const int bsize = blocks->size(b);
      bF_lt->block(b).assign_subblock(bF->block(b),0,bsize-1,0,bsize-1);
    }
  }
#define TEST_CABS_CANONICAL 0
#if TEST_CABS_CANONICAL
  F_cabs.print("Fock matrix in the original CABS basis");
  F_cabs_lt.print("Fock matrix in the original CABS basis");
#endif

  std::string id_sb, id;
  if (spin_polarized()) {
    id = ParsedOrbitalSpaceKey::key(std::string("c'"),Alpha);
  }
  else {
    id = "c'";
  }
  std::ostringstream oss;  oss << "CABS (" << to_string(spin) << ";canonicalized)";
  cabs_canonical[spin] = new OrbitalSpace(id, oss.str(),
                                          cabs->coefs()*F_cabs_lt.eigvecs(),
                                          cabs->basis(),
                                          cabs->integral(),
                                          F_cabs_lt.eigvals(),
                                          0, 0,
                                          OrbitalSpace::symmetry);

  const Ref<OrbitalSpaceRegistry> idxreg = this->orbital_registry();
  idxreg->add(make_keyspace_pair(cabs_canonical[spin]));

#if TEST_CABS_CANONICAL
  {
    RefSCMatrix F_CC = fock(cabs_canonical[spin],cabs_canonical[spin],spin);
    F_CC.print("Fock matrix in the CABS canonical basis");
    cabs_canonical[spin]->evals().print("eigenvalues of CABS canonical basis");
  }
#endif

  return cabs_canonical[spin];
}


////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
