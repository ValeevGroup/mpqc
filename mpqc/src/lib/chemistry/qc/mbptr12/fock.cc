//
// fock.cc
//
// Copyright (C) 2004 Edward Valeev
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdexcept>
#include <util/class/scexception.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/symmint.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/print.h>

using namespace std;
using namespace sc;

#define TEST_FOCKBUILD 0
#define NEW_HCORE 0

RefSCMatrix
R12IntEval::fock_(const Ref<MOIndexSpace>& bra_space,
                  const Ref<MOIndexSpace>& ket_space,
                  SpinCase1 spin,
                  double scale_J, double scale_K)
{
  Ref<SingleRefInfo> refinfo = r12info()->refinfo();
  const Ref<GaussianBasisSet> bs1 = bra_space->basis();
  const Ref<GaussianBasisSet> bs2 = ket_space->basis();
  const bool bs1_eq_bs2 = (bs1 == bs2);
  int nshell1 = bs1->nshell();
  int nshell2 = bs2->nshell();

  RefSCMatrix vec1t = bra_space->coefs().t();
  RefSCMatrix vec2 = ket_space->coefs();

#if NEW_HCORE
  RefSCDimension aodim1 = vec1t.coldim();
  RefSCDimension aodim2 = vec2.rowdim();
  Ref<SCMatrixKit> sokit = bs1->so_matrixkit();

  // cast Wavefunction to MBPT2_R12
  MBPT2_R12 *r12wfn = dynamic_cast<MBPT2_R12*>(r12info()->wfn());
  if (r12wfn == 0) {
      throw ProgrammingError(
          "r12info()->wfn() was not an MBPT2_R12 object",
          __FILE__, __LINE__, class_desc());
  }
  Ref<SCF> ref = r12wfn->ref();
  // Form the DK correction in the current basis using the momentum
  // basis of the reference wavefunction.  The momentum basis in the
  // reference should be a superset of hcore_basis
  Ref<GaussianBasisSet> p_basis = ref->momentum_basis();
  Ref<GaussianBasisSet> hcore_basis;
  if (bs1_eq_bs2) {
      hcore_basis = bs1;
    }
  else {
      hcore_basis = bs1 + bs2;
    }

  RefSymmSCMatrix hsymm
      = ref->core_hamiltonian_for_basis(hcore_basis,p_basis);


  // convert hsymm to the AO basis
  Ref<Integral> localints = r12info_->integral()->clone();
  localints->set_basis(hcore_basis,hcore_basis);
  Ref<PetiteList> hcore_pl = localints->petite_list();
  RefSymmSCMatrix hsymm_ao = hcore_pl->to_AO_basis(hsymm);
  hsymm = 0;

  RefSCMatrix h(aodim1, aodim2, sokit);
  if (bs1_eq_bs2) {
      h.assign(0.0);
      h.accumulate(hsymm_ao);
    }
  else {
      int br = 0;
      int er = bs1->nbasis() - 1;
      int bc = 0;
      int ec = bs2->nbasis() - 1;
      int source_br = 0;
      int source_bc = bs1->nbasis();
      RefSCMatrix hrect_ao(hsymm_ao.dim(), hsymm_ao.dim(), sokit);
      hrect_ao.assign(0.0);
      hrect_ao.accumulate(hsymm_ao);
      h.assign_subblock(hrect_ao, br, er, bc, ec, source_br, source_bc);
    }

#else // ! NEW_HCORE

  if (r12info()->wfn()->dk() > 0) {
      throw ProgrammingError(
	  "cannot use reference with dk > 0 with old mp2r12 fock algorithm",
	  __FILE__, __LINE__, class_desc());
  }

  Ref<Integral> localints = r12info_->integral()->clone();
  localints->set_basis(bs1,bs2);

  Ref<OneBodyInt> h_ints = localints->hcore();

  // form AO moment matrices
  RefSCDimension aodim1 = vec1t.coldim();
  RefSCDimension aodim2 = vec2.rowdim();
  Ref<SCMatrixKit> aokit = bs1->so_matrixkit();
  RefSCMatrix h(aodim1, aodim2, aokit);
  h.assign(0.0);
  
  for(int sh1=0; sh1<nshell1; sh1++) {
    int bf1_offset = bs1->shell_to_function(sh1);
    int nbf1 = bs1->shell(sh1).nfunction();

    int sh2max;
    if (bs1_eq_bs2)
      sh2max = sh1;
    else
      sh2max = nshell2-1;
    
    for(int sh2=0; sh2<=sh2max; sh2++) {
      int bf2_offset = bs2->shell_to_function(sh2);
      int nbf2 = bs2->shell(sh2).nfunction();
      
      h_ints->compute_shell(sh1,sh2);
      const double *hintsptr = h_ints->buffer();

      int bf1_index = bf1_offset;
      for(int bf1=0; bf1<nbf1; bf1++, bf1_index++, hintsptr+=nbf2) {
	int bf2_index = bf2_offset;
	const double *ptr = hintsptr;
	int bf2max;
        if (bs1_eq_bs2 && sh1 == sh2)
          bf2max = bf1;
        else
	  bf2max = nbf2-1;
	for(int bf2=0; bf2<=bf2max; bf2++, bf2_index++, ptr++) {

	  h.set_element(bf1_index, bf2_index, *ptr);

        }
      }
    }
  }

  // Symmetrize matrices, if necessary
  if (bs1_eq_bs2) {
    const int nbasis = bs1->nbasis();
    for(int bf1=0; bf1<nbasis; bf1++)
      for(int bf2=0; bf2<=bf1; bf2++) {
        h(bf2,bf1) = h(bf1,bf2);
      }
  }

  h_ints = 0;

#endif // NEW_HCORE

  // finally, transform
  RefSCMatrix F = vec1t * h * vec2;
  if (debug_ >= DefaultPrintThresholds::allN2)
    F.print("Core Hamiltonian contribution");
  // and clean up a bit
  h = 0;
  
  const bool spin_unrestricted = spin_polarized();
  
#if TEST_FOCKBUILD
  RefSCMatrix Ftest;
  if (!spin_unrestricted)
    Ftest = F.clone();
#endif
  
  //
  // add coulomb and exchange parts
  //
  const double J_prefactor = (spin_unrestricted ? 1.0 : 2.0);
  for(int s=0; s<nspincases1(); s++) {
    const SpinCase1 sc = static_cast<SpinCase1>(s);
    if (occ(sc)->rank() == 0)
      continue;
    if (scale_J != 0.0) {
      RefSCMatrix J = coulomb_(occ(sc),bra_space,ket_space);
      J.scale(J_prefactor*scale_J);
      if (debug_ >= DefaultPrintThresholds::allN2)
        J.print("Coulomb contribution");
      F.accumulate(J); J = 0;
    }
  }
  if (scale_K != 0.0) {
    RefSCMatrix K = exchange_(occ(spin),bra_space,ket_space);
    K.scale(-1.0*scale_K);
    if (debug_ >= DefaultPrintThresholds::allN2)
      K.print("Exchange contribution");
    F.accumulate(K); K = 0;
  }
  
  if (debug_ >= DefaultPrintThresholds::allN2) {
    F.print("Fock matrix");
  }
  
#if TEST_FOCKBUILD
  if (scale_J == 1.0 && scale_K == 1.0) {
    Ref<GaussianBasisSet> obs = occ(Alpha)->basis();
    Ref<FockContribution> fc = new CLHFContribution(bs1,bs2,obs,obs);
    
    // transform the density matrix to the AO basis
    RefSymmSCMatrix D_ao = reference()->ao_density();
    RefSCMatrix G_ao = Ftest.kit()->matrix(aodim1,aodim2);
    fc->set_fmat(0, G_ao);
    fc->set_pmat(0, D_ao);
    
    const Ref<Integral>& integral = r12info()->integral();
    const unsigned int nthreads = r12info()->thr()->nthread();
    Ref<TwoBodyInt>* tbints = new Ref<TwoBodyInt>[nthreads];
    for(unsigned int t=0; t<nthreads; ++t)
      tbints[t] = integral->electron_repulsion();
    
    // WARNING should replace with something reasonable
    const double accuracy = 1.0e-15;
    Ref<FockBuild> fb = new FockBuild(fc,accuracy,tbints,
                                      bs1, bs2, obs, obs,
                                      r12info()->msg(), r12info()->thr(),
                                      integral);
    fb->build();
    ExEnv::out0() << indent << scprintf("%20.0f integrals\n",
    fb->contrib()->nint());
    fb = 0;
    
    // WARNING Need to symmetrize G
    
    for(unsigned int t=0; t<nthreads; ++t)
      tbints[t] = 0;
    delete[] tbints;
  }
#endif
  
  return F;
}

///////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
