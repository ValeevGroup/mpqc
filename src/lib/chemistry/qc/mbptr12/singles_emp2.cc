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
#include <cassert>

#include <mpqc_config.h>
#include <util/misc/formio.h>
#include <util/misc/regtime.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/state_text.h>
#include <util/state/state_bin.h>
#include <math/scmat/matrix.h>
#include <math/scmat/blocked.h>
#include <math/optimize/conjgrad.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <math/distarray4/distarray4.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <math/mmisc/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <math/scmat/svd.h>

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

    // !!! because frozen core/vir orbitals are rotated in semicanonicalization, must use full spaces here!
    Ref<OrbitalSpace> occ = this->occ(spin);
    Ref<OrbitalSpace> vir = this->vir(spin);
    RefSCMatrix Fia = this->fock(occ,vir,spin);

    const int ni = occ->rank();
    const int na = vir->rank();
    const RefDiagSCMatrix& ievals = occ->evals();
    const RefDiagSCMatrix& aevals = vir->evals();

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

#if 0
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
#endif

namespace {
  // this functor helps to implement steepest descent CABS singles solver
  struct CABS_singles_residual {
      CABS_singles_residual(const RefSCMatrix& h1,
                            const RefSCMatrix& h0_ab,
                            const RefSCMatrix& h0_ij) : H1(h1), H0_AB(h0_ab), H0_IJ(h0_ij)
      {
      }

      const RefSCMatrix& H1;
      const RefSCMatrix& H0_AB;
      const RefSCMatrix& H0_IJ;

      void operator()(const RefSCMatrix& T1, RefSCMatrix& R1) {
        R1.assign(H1);
        R1.accumulate_product(T1, H0_AB);
        RefSCMatrix tmp = H0_IJ * T1; tmp.scale(-1.0);
        R1.accumulate( tmp );
      }
  };

  /** this functor helps to implement conjugate gradient CABS singles solver
   */
  struct CABS_singles_h0t1 {
      /**
       * @param h0_AB allvirt/allvirt Fock operator
       * @param h0_ij occ/occ Fock operator
       */
      CABS_singles_h0t1(const RefSCMatrix& h0_AB,
                        const RefSCMatrix& h0_ij) : H0_AB(h0_AB), H0_IJ(h0_ij)
      {
      }

      const RefSCMatrix& H0_AB;
      const RefSCMatrix& H0_IJ;

      /**
       * @param[in] T1 t_i^A
       * @param[out] R1 R_i^A
       */
      void operator()(const RefSCMatrix& T1, RefSCMatrix& R1) {
        R1.assign(0.0);
        R1.accumulate_product(T1, H0_AB);
        RefSCMatrix tmp = H0_IJ * T1; tmp.scale(-1.0);
        R1.accumulate( tmp );
      }
  };
}

double
R12IntEval::compute_emp2_cabs_singles_noncanonical(bool vir_cabs_coupling) {

  // prerequsites:
  // H1_A^i = F_A^i where A is all virtual orbitals (allvirt = virt + cabs)
  // H0_A^B = F_A^B
  // H0_i^j = F_i^j
  // compute:
  // R_A^i =  F_A^B t_A^i - F_j^i t^j_A + H1_A^i
  // update:
  // t_A^i -= R_A^i / (F_A^A - F_i^i)
  // energy correction:
  // E2 = t_A^i H1^A_i

  Ref<MessageGrp> msg = r12world()->world()->msg();
  Ref<MemoryGrp> mem = r12world()->world()->mem();
  Ref<ThreadGrp> thr = r12world()->world()->thr();

  Timer tim("(2)_S energy");
  std::string evalname("(2)_S energy evaluator");
  ExEnv::out0() << endl << indent << "Entered " << evalname << endl;
  ExEnv::out0() << incindent;

  double result = 0.0;
  for(int s=0; s<nspincases1(); s++) {
    const SpinCase1 spin = static_cast<SpinCase1>(s);

    Ref<OrbitalSpace> occ = this->occ(spin);
    if (occ->rank() == 0) continue;
    Ref<OrbitalSpace> vir = this->vir(spin);
    Ref<OrbitalSpace> cabs = r12world()->cabs_space(spin);
    // test for CABS Singles contribution to dipole moment
    //Ref<OrbitalSpace> cabs = cabs_space_canonical(spin);
    if (cabs->rank() == 0) continue;
    Ref<OrbitalSpace> aspace = cabs;
    if (vir_cabs_coupling) {
      std::string allvirt_key = ParsedOrbitalSpaceKey::key(std::string("A'"),this->spin_polarized() ? spin : AnySpinCase1);
      if (not this->orbital_registry()->key_exists(allvirt_key)) { // this should be the first time this is created
        Ref<OrbitalSpaceUnion> allvirt = new OrbitalSpaceUnion(allvirt_key,
                                                               prepend_spincase(spin, std::string("all virtual orbitals")),
                                                               *vir, *cabs, true);
        this->orbital_registry()->add(allvirt->id(), allvirt);
        aspace = allvirt;
        // this space potentially uses a union of OBS (or VBS) and ABS
        // make sure that this AO basis is known
        if (this->ao_registry()->key_exists(allvirt->basis()) == false) {
          Ref<OrbitalSpace> muprime = new AtomicOrbitalSpace("mu''", "VBS(AO)+CABS(AO)", allvirt->basis(), allvirt->integral());
          this->orbital_registry()->add(make_keyspace_pair(muprime));
          this->ao_registry()->add(muprime->basis(),muprime);
        }
      }
      else
        aspace = orbital_registry()->value(allvirt_key);
    }

    RefSCMatrix FiA = fock(occ,aspace,spin);
    RefSCMatrix Fii = fock(occ,occ,spin);
    RefSCMatrix FAA = fock(aspace,aspace,spin);
    T1_cabs_[s] = FiA.clone(); T1_cabs_[s].assign(0.0);

    // pre-compute preconditioner: PC(A,i) = 1/(FAA-Fii)
    RefSCMatrix PC = T1_cabs_[s].clone();
    const unsigned ni = occ->rank();
    const unsigned nA = aspace->rank();
    for (unsigned i = 0; i < ni; ++i) {
      for (unsigned A = 0; A < nA; ++A) {
        PC(i, A) = 1.0 / (FAA(A, A) - Fii(i, i));
      }
    }

    bool converged = true;
    { // use conjugate gradient
      CABS_singles_h0t1 h0t1(FAA, Fii);
      RefSCMatrix rhs = FiA.copy();
      rhs.scale(-1.0);
      double Rnorm2;
      try {
        Rnorm2 = linsolv_conjugate_gradient(h0t1, rhs, T1_cabs_[s], PC,
                                            1e-10);
        //ConjugateGradientSolver<RefSCMatrix, CABS_singles_h0t1> cg;
        //Rnorm2 = cg(h0t1, rhs, T1_cabs_[s], PC,
        //            1e-10);
      }
      catch (...) {
        // if failed for some reason, at least compute the energy to help with troubleshooting
        Ref<SCElementScalarProduct> dotprodop = new SCElementScalarProduct;
        T1_cabs_[s].element_op(dotprodop, FiA);
        const double E2 = dotprodop->result();
        std::cout << indent << "WARNING: CG solver for (2)_S ("
                  << prepend_spincase(spin, std::string("spin")) << ") did not converge, using direct solver" << std::endl;
        converged = false; // do not rethrow ...
      }

      // ... try a direct solver as a last resort
      if (not converged) {
        const int niA = ni*nA;
        RefSCDimension iA_dim = new SCDimension(niA, 1, &niA);
        iA_dim->blocks()->set_subdim(0, new SCDimension(niA));

        RefSymmSCMatrix C = T1_cabs_[s].kit()->symmmatrix(iA_dim);
        C.assign(0.0);
        for (unsigned i = 0, iA = 0; i < ni; ++i) {
          for (unsigned A = 0; A < nA; ++A, ++iA) {
            unsigned int iB = i*nA;
            for (unsigned B = 0; B <= A; ++B, ++iB) {
              C.set_element(iA, iB, FAA(A, B));
            }
          }
        }
        for (unsigned i = 0, iA = 0; i < ni; ++i) {
          for (unsigned j = 0; j <= i; ++j) {
            unsigned int iA = i*nA;
            unsigned int jA = j*nA;
            for (unsigned A = 0; A < nA; ++A, ++iA, ++jA) {
              const double newval = C(iA, jA) - Fii(i, j);
              C(iA, jA) = newval;
            }
          }
        }

        RefSCVector B = C.kit()->vector(C.dim());
        for (unsigned i = 0, iA = 0; i < ni; ++i) {
          for (unsigned A = 0; A < nA; ++A, ++iA) {
            B(iA) = rhs(i, A);
          }
        }

        RefSCVector X = B.clone();
        lapack_linsolv_symmnondef(C, X, B);

        for (unsigned i = 0, iA = 0; i < ni; ++i) {
          for (unsigned A = 0; A < nA; ++A, ++iA) {
            T1_cabs_[s](i, A) = X(iA);
          }
        }
      }

      Ref<SCElementScalarProduct> dotprodop = new SCElementScalarProduct;
      T1_cabs_[s].element_op(dotprodop, FiA);
      result += dotprodop->result();
    }

    if (! spin_polarized()) {
      result *= 2.0;
      T1_cabs_[Beta] = T1_cabs_[Alpha];
    }

  } // end of spin loop

  ExEnv::out0() << indent << "E(MP2 " << (vir_cabs_coupling ? "OBS+CABS" : "CABS") << " singles = (2)_S) = " << scprintf("%25.15lf",result) << endl;
  ExEnv::out0() << decindent << indent << "Exited " << evalname << endl;

  return result;
}

double
R12IntEval::compute_emp2_cabs_singles_noncanonical_ccsd(const RefSCMatrix& T1_ia_alpha,
                                                        const RefSCMatrix& T1_ia_beta) {

  // prerequsites:
  // H1_A^i = F_a'^i where a' is cabs
  // H1_a^B = F_a^b'
  // H0_A^B = F_a'^b'
  // H0_i^j = F_i^j
  // compute:
  // R_A^i =  F_A^B t_A^i - F_j^i t^j_A + F_A^a t^i_a + H1_A^i
  // update:
  // t_A^i -= R_A^i / (F_A^A - F_i^i)
  // energy correction:
  // E2 = t_A^i H1^A_i + t_A^i H1^A_a t^a_i

  Ref<MessageGrp> msg = r12world()->world()->msg();
  Ref<MemoryGrp> mem = r12world()->world()->mem();
  Ref<ThreadGrp> thr = r12world()->world()->thr();

  Timer tim("(2)_S energy");
  std::string evalname("(2)_S energy evaluator");
  ExEnv::out0() << endl << indent << "Entered " << evalname << endl;
  ExEnv::out0() << incindent;

  double result = 0.0;
  for(int s=0; s<nspincases1(); s++) {
    const SpinCase1 spin = static_cast<SpinCase1>(s);

    Ref<OrbitalSpace> occ = this->occ(spin);
    if (occ->rank() == 0) continue;
    Ref<OrbitalSpace> vir = this->vir(spin);
    Ref<OrbitalSpace> cabs = r12world()->cabs_space(spin);
    if (cabs->rank() == 0) continue;

    RefSCMatrix FiA = fock(occ,cabs,spin);
    RefSCMatrix Fii = fock(occ,occ,spin);
    RefSCMatrix FaA = fock(vir,cabs,spin);
    RefSCMatrix FAA = fock(cabs,cabs,spin);

    const RefSCMatrix& T1_ia_ref = (spin == Alpha) ? T1_ia_alpha : T1_ia_beta;
    RefSCMatrix Fia = fock(occ,vir,spin);
    RefSCMatrix T1_ia = Fia.clone();
    // doesn't work because T1_ia_ref is local, T1_ia is replicated
    //T1_ia->convert(T1_ia_ref);
    T1_ia->assign(0.0);
    const int nocc_act = T1_ia_ref.nrow();
    const int nvir_act = T1_ia_ref.ncol();
    Ref<OrbitalSpace> occ_act = this->occ_act(spin);
    Ref<OrbitalSpace> vir_act = this->vir_act(spin);
    MOIndexMap map_act2occ(*occ << *occ_act);
    MOIndexMap map_act2vir(*vir << *vir_act);
    for(int i=0; i<nocc_act; ++i) {
      for(int a=0; a<nvir_act; ++a) {
        T1_ia(map_act2occ[i], map_act2vir[a]) = T1_ia_ref(i,a);
      }
    }

    T1_cabs_[s] = FiA.clone(); T1_cabs_[s].assign(0.0);

    // pre-compute preconditioner: PC(A,i) = 1/(FAA-Fij)
    RefSCMatrix PC = T1_cabs_[s].clone();
    const unsigned ni = occ->rank();
    const unsigned nA = cabs->rank();
    for (unsigned i = 0; i < ni; ++i) {
      for (unsigned A = 0; A < nA; ++A) {
        PC(i, A) = 1.0 / (FAA(A, A) - Fii(i, i));
      }
    }

    bool converged = false;
    { // use conjugate gradient
      CABS_singles_h0t1 h0t1(FAA, Fii);
      RefSCMatrix rhs = FiA.copy();
      rhs.accumulate_product(T1_ia, FaA);
      rhs.scale(-1.0);  // rhs = - (F_i^a' + t_i^a F_a^a')
      double Rnorm2;
      try {
        Rnorm2 = linsolv_conjugate_gradient(h0t1, rhs, T1_cabs_[s], PC,
                                            1e-10);
      }
      catch (...) {
        // if failed for some reason, at least compute the energy to help with troubleshooting
        Ref<SCElementScalarProduct> dotprodop = new SCElementScalarProduct;
        T1_cabs_[s].element_op(dotprodop, FiA);
        const double E2 = dotprodop->result();
        std::cout << "nonconverged (2)_S energy = " << E2
            << " (elapsed spincase " << (spin == Alpha ? "alpha" : "alpha+beta") << ")"
            << std::endl;
        throw;
      }
      Ref<SCElementScalarProduct> dotprodop = new SCElementScalarProduct;
      T1_cabs_[s].element_op(dotprodop, rhs);
      // - because rhs included a - sign
      result -= dotprodop->result();
    }

    if (! spin_polarized()) {
      result *= 2.0;
      T1_cabs_[Beta] = T1_cabs_[Alpha];
    }

  } // end of spin loop

  ExEnv::out0() << indent << "E(MP2 CABS singles = (2)_S) = " << scprintf("%25.15lf",result) << endl;
  ExEnv::out0() << decindent << indent << "Exited " << evalname << endl;

  return result;
}

// rotates CABS so that the Fock matrix is diagonal
const Ref<OrbitalSpace>&
R12IntEval::cabs_space_canonical(SpinCase1 spin)
{
  static Ref<OrbitalSpace> cabs_canonical[] = {0, 0};

  if (!spin_polarized() && spin == Beta)
    return cabs_space_canonical(Alpha);
  if (cabs_canonical[spin] == 0)
    cabs_canonical[spin] = this->cabs_space_fockcanonical(spin,1.0,1.0,1.0);

  return cabs_canonical[spin];
}

// rotates CABS so that the Fock matrix is diagonal
const Ref<OrbitalSpace>&
R12IntEval::cabs_space_hcanonical(SpinCase1 spin)
{
  static Ref<OrbitalSpace> cabs_hcanonical[] = {0, 0};

  if (!spin_polarized() && spin == Beta)
    return cabs_space_hcanonical(Alpha);
  if (cabs_hcanonical[spin] == 0)
    cabs_hcanonical[spin] = this->cabs_space_fockcanonical(spin,1.0,0.0,0.0);

  return cabs_hcanonical[spin];
}

// rotates CABS so that the Fock matrix is diagonal
Ref<OrbitalSpace>
R12IntEval::cabs_space_fockcanonical(SpinCase1 spin,
                                     double scale_H,
                                     double scale_J,
                                     double scale_K)
{
  MPQC_ASSERT(r12world()->obs_eq_ribs() == false);

  Ref<OrbitalSpace> cabs = r12world()->cabs_space(spin);
  const int ncabs = cabs->rank();

  // note that I'm overriding pauli flag here -- true Fock matrix must always be used
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
      if (bsize != 0)
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
    id = ParsedOrbitalSpaceKey::key(std::string("c'"), spin);
  }
  else {
    id = "c'";
  }
  std::ostringstream oss;  oss << "CABS (" << to_string(spin) << ";canonicalized)";
  Ref<OrbitalSpace> cabs_canonical = new OrbitalSpace(id, oss.str(),
                                                      cabs->coefs()*F_cabs_lt.eigvecs(),
                                                      cabs->basis(),
                                                      cabs->integral(),
                                                      F_cabs_lt.eigvals(),
                                                      0, 0,
                                                      OrbitalSpace::symmetry);

  const Ref<OrbitalSpaceRegistry> idxreg = this->orbital_registry();
  idxreg->add(make_keyspace_pair(cabs_canonical));

#if TEST_CABS_CANONICAL
  {
    RefSCMatrix F_CC = fock(cabs_canonical[spin],cabs_canonical[spin],spin);
    F_CC.print("Fock matrix in the CABS canonical basis");
    cabs_canonical[spin]->evals().print("eigenvalues of CABS canonical basis");
  }
#endif

  return cabs_canonical;
}


////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
