//
// etrain.cc
//
// Copyright (C) 2011 Edward Valeev
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

#include <cassert>
#include <util/keyval/keyval.h>
#include <util/group/message.h>
#include <util/group/pregtime.h>
#include <util/group/thread.h>
#include <util/group/memory.h>
#include <util/misc/bug.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <util/state/state_bin.h>
#include <util/options/GetLongOpt.h>
#include <math/scmat/repl.h>
#include <math/scmat/dist.h>
#include <chemistry/molecule/molshape.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/wfn/orbital.h>
#include <chemistry/qc/wfn/orbitalspace_utils.h>
#include <chemistry/qc/etrain/etrain.h>
#include <chemistry/qc/lcao/fockbuilder.h>

using namespace sc;

// This object maintains information about class ETraIn.
// NOTE: The last three arguments register default, KeyVal and StateIn constructors,
//       respectively. Hence only the second of the three arguments is not zero
static ClassDesc ETraIn_cd(typeid(ETraIn), "ETraIn", 1,
                                   "public Function",
                                   0, create<ETraIn>, 0);


ETraIn::ETraIn(const Ref<KeyVal>& keyval): Function(keyval)
{
  obwfn12_ << keyval->describedclassvalue("wfn12");
  obwfn1_ << keyval->describedclassvalue("wfn1");
  obwfn2_ << keyval->describedclassvalue("wfn2");

  //
  // Check wave functions. Must be for closed-shell systems and derived from OneBodyWavefunction
  //
  if (obwfn12_.null())
    throw InputError("wfn12 keyword not specified or has wrong type (must be derived from OneBodyWavefunction)", __FILE__, __LINE__);
  if (obwfn1_.null())
    throw InputError("wfn1 keyword not specified or has wrong type (must be derived from OneBodyWavefunction)", __FILE__, __LINE__);
  if (obwfn2_.null())
    throw InputError("wfn2 keyword not specified or has wrong type (must be derived from OneBodyWavefunction)", __FILE__, __LINE__);
  if (obwfn12_->nelectron()%2)
    throw InputError("wfn12 wave function must be of closed-shell type");
  if (obwfn1_->nelectron()%2)
    throw InputError("wfn1 wave function must be of closed-shell type");
  if (obwfn2_->nelectron()%2)
    throw InputError("wfn2 wave function must be of closed-shell type");

  // Must use canonical orthogonalization method since symmetric does not drop redundant combinations
  if (obwfn12_->orthog_method() != OverlapOrthog::Canonical ||
      obwfn1_->orthog_method() != OverlapOrthog::Canonical ||
      obwfn2_->orthog_method() != OverlapOrthog::Canonical) {
    throw InputError("all Wavefunctions must use canonical orthogonalization method",__FILE__,__LINE__);
  }

  // known monomer ionization potentials?
  if (keyval->exists("ip1")) {
    read_ip(keyval, "ip1", ip1_, obwfn1_->nelectron()/2, ip1_orbs_);
  }
  if (keyval->exists("ip2")) {
    read_ip(keyval, "ip2", ip2_, obwfn2_->nelectron()/2, ip2_orbs_);
  }

  if (keyval->exists("grid")) {
    grid_ << keyval->describedclassvalue("grid");
    if (grid_.null()) { // check if grid = auto was used
      const std::string make_grid = keyval->stringvalue("grid", KeyValValuestring(""));
      if (make_grid == "auto") { // construct the grid automatically
        const Ref<VDWShape> vdwshape = new VDWShape(obwfn12_->molecule());
        SCVector3 gmin, gmax;
        vdwshape->boundingbox(-1.0, 1.0, gmin, gmax);
        SCVector3 gorigin = gmin;
        SCVector3 gsize = gmax - gmin;
        const double resolution = 0.2;
        int n[3]; for(int i=0; i<3; ++i) n[i] = int(std::ceil(gsize[i] / 0.2));
        SCVector3 axis0(gsize[0]/n[0], 0.0, 0.0);
        SCVector3 axis1(0.0, gsize[1]/n[1], 0.0);
        SCVector3 axis2(0.0, 0.0, gsize[2]/n[2]);
        grid_ = new Grid(n[0], n[1], n[2], gorigin, axis0, axis1, axis2);
      }
    }
  }

  // nocc_ and nuocc_ default to -1, which means to use all monomer orbitals
  nocc_  = keyval->intvalue("nocc",KeyValValueint(-1));
  nuocc_  = keyval->intvalue("nuocc",KeyValValueint(-1));
  debug_  = keyval->intvalue("debug",KeyValValueint(0));

  // check that monomer atoms are subsets of the n-mer and construct atom maps
  // since Molecule does not rotate atoms, the difference in coordinates can be only due to different
  // origins
  Ref<Molecule> mol1 = obwfn1_->molecule();
  Ref<Molecule> mol2 = obwfn2_->molecule();
  Ref<Molecule> mol12 = obwfn12_->molecule();
  // monomer 1
  {
    SCVector3 shift1 = mol12->ref_origin() - mol1->ref_origin();
    const int natom1 = mol1->natom();
    atom_map1_.resize(natom1);
    for (int a1 = 0; a1 < natom1; ++a1) {
      double r[3];
      for (int xyz = 0; xyz < 3; ++xyz)
        r[xyz] = mol1->r(a1, xyz) + shift1[xyz];
      int a12 = mol12->atom_at_position(r, 1e-4);
      if (a12 < 0) {
        mol1->print();
        mol12->print();
        std::ostringstream oss;
        oss << "atom " << a1 << " in monomer 1 not found in the n-mer";
        throw InputError(oss.str().c_str(), __FILE__, __LINE__);
      }
      atom_map1_[a1] = a12;
    }
  }
  // monomer 2
  {
    SCVector3 shift2 = mol12->ref_origin() - mol2->ref_origin();
    const int natom2 = mol2->natom();
    atom_map2_.resize(natom2);
    for (int a2 = 0; a2 < natom2; ++a2) {
      double r[3];
      for (int xyz = 0; xyz < 3; ++xyz)
        r[xyz] = mol2->r(a2, xyz) + shift2[xyz];
      int a12 = mol12->atom_at_position(r, 1e-4);
      if (a12 < 0) {
        mol2->print();
        mol12->print();
        std::ostringstream oss;
        oss << "atom " << a2 << " in monomer 2 not found in the n-mer";
        throw InputError(oss.str().c_str(), __FILE__, __LINE__);
      }
      atom_map2_[a2] = a12;
    }
  }
}

void
ETraIn::obsolete() {
  Function::obsolete();
  obwfn12_->obsolete();
  obwfn1_->obsolete();
  obwfn2_->obsolete();
}

void
ETraIn::run(void) {
  this->compute();
}

void
ETraIn::compute(void)
{
  if(gradient_needed()) {
    ExEnv::out0() << "No gradients can be provided for this object" << std::endl;
    abort();
  }

  obwfn12_->set_desired_value_accuracy(desired_value_accuracy());
  obwfn1_->set_desired_value_accuracy(desired_value_accuracy());
  obwfn2_->set_desired_value_accuracy(desired_value_accuracy());
  double energy12 = obwfn12_->energy();
  double energy1 = obwfn1_->energy();
  double energy2 = obwfn2_->energy();
  compute_train();

  set_value(0.0);
  set_actual_value_accuracy(0.0);
}

// based on LocalDiagSCMatrix::vprint
void
ip_print(const char *title1, const char *title2, std::ostream& os, int prec, double conv, RefDiagSCMatrix m1, RefDiagSCMatrix m2)
{
  MPQC_ASSERT(m1.n() == m2.n()); // meant to print 2 matrices side by side

  double max = std::max(m1->maxabs()*conv, m2->maxabs()*conv);
  max = (max==0.0) ? 1.0 : log10(max);
  if (max < 0.0) max=1.0;

  const int twidth = std::max(strlen(title1), strlen(title2));
  const int lwidth = std::max(prec + 5 + (int) max, twidth);

  os << std::endl << indent << scprintf("      %*s %*s\n", lwidth, title1, lwidth, title2);

  if (m1.n()==0) {
    os << indent << "empty matrices\n";
    return;
  }

  for (int i=0; i<m1.n(); i++)
    os << indent
       << scprintf("%5d %*.*f %*.*f\n",
                   i+1,
                   lwidth,prec,m1.get_element(i)*conv,
                   lwidth,prec,m2.get_element(i)*conv
                  );
  os << std::endl;

  os.flush();
}

namespace sc {
  template <typename T>
  T max(T t1, T t2, T t3) {
    return std::max(std::max(t1,t2),t3);
  }
}

void
ip_print(const char *title1, const char *title2, const char *title3,
       std::ostream& os, int prec, double conv,
       RefDiagSCMatrix m1, RefDiagSCMatrix m2, RefDiagSCMatrix m3
      )
{
  MPQC_ASSERT(m1.n() == m2.n()); // meant to print 3 matrices side by side
  MPQC_ASSERT(m1.n() == m3.n());

  double max = sc::max(m1->maxabs()*conv, m2->maxabs()*conv, m3->maxabs()*conv);
  max = (max==0.0) ? 1.0 : log10(max);
  if (max < 0.0) max=1.0;

  const int twidth = sc::max(strlen(title1), strlen(title2), strlen(title3));
  const int lwidth = std::max(prec + 5 + (int) max, twidth);

  os << std::endl << indent << scprintf("      %*s %*s %*s\n", lwidth, title1, lwidth, title2, lwidth, title3);

  if (m1.n()==0) {
    os << indent << "empty matrices\n";
    return;
  }

  for (int i=m1.n()-1; i >= 0; i--) {
    const double ip1 = m1.get_element(i)*conv;
    const double ip2 = m2.get_element(i)*conv;
    const double ip3 = m3.get_element(i)*conv;
    if (ip1 >= 0.0 && ip2 >= 0.0 && ip3 >= 0.0)
      os << indent
         << scprintf("%5d %*.*f %*.*f %*.*f\n",
                     m1.n()-i,
                     lwidth,prec,ip1,
                     lwidth,prec,ip2,
                     lwidth,prec,ip3
                    );
  }
  os << std::endl;

  os.flush();
}

void
ETraIn::compute_train()
{
  // Effective Hailtonian for SCF and other (ExtendedHuckel) wave functions is recovered via different members
  // thus need to determine which wave function we have
  Ref<CLSCF> obwfn12_clscf; obwfn12_clscf << obwfn12_;

  if (! obwfn12_->integral()->equiv(obwfn1_->integral()) ||
      ! obwfn12_->integral()->equiv(obwfn2_->integral()) )
    throw InputError("Integral factories must match for all calculations",__FILE__,__LINE__);
  const Ref<Integral>& integral = obwfn12_->integral();

  RefSCMatrix vec12_so = obwfn12_->eigenvectors();
  RefSCMatrix vec1_so = obwfn1_->eigenvectors();
  RefSCMatrix vec2_so = obwfn2_->eigenvectors();
  Ref<PetiteList> plist12 = obwfn12_->integral()->petite_list();
  Ref<PetiteList> plist1 = obwfn1_->integral()->petite_list();
  Ref<PetiteList> plist2 = obwfn2_->integral()->petite_list();
  RefSCMatrix vec12 = plist12->evecs_to_AO_basis(vec12_so);
  RefSCMatrix vec1 = plist1->evecs_to_AO_basis(vec1_so);
  RefSCMatrix vec2 = plist2->evecs_to_AO_basis(vec2_so);

  // map vec1 and vec2 to the n-mer basis. Because Molecule does not rotate the frame, as long as
  // monomers are given in the same frame as the n-mer, there is no need to rotate the basis functions.
  // we thank ye, MPQC gods!
  //
  // because frames may be shifted, however, can't just say bs12 << bs1
  Ref<GaussianBasisSet> bs1 = obwfn1_->basis();
  Ref<GaussianBasisSet> bs2 = obwfn2_->basis();
  Ref<GaussianBasisSet> bs12 = obwfn12_->basis();
  std::vector<unsigned int> basis_map1(bs1->nbasis());  // maps basis functions of bs1 to bs12
  std::vector<unsigned int> basis_map2(bs2->nbasis());  // etc.
  {  // map bs1 to bs12
    const int natom = atom_map1_.size();
    for(int a=0; a<natom; ++a) {
      const int a12 = atom_map1_[a];

      const int nshell = bs1->nshell_on_center(a);
      for(int s=0; s<nshell; ++s) {
        const int S = bs1->shell_on_center(a, s);
        int s12;
        try {
          s12 = ishell_on_center(a12, bs12, bs1->shell(S));
          const int nbf = bs1->shell(S).nfunction();
          const int offset = bs1->shell_to_function(S);
          const int offset12 = bs12->shell_to_function(s12);
          for(int f=0; f<nbf; ++f)
            basis_map1[f + offset] = f + offset12;
        }
        catch (ProgrammingError&) {
          std::ostringstream oss; oss << "shell " << S << " in basis set of scf1 is not found in that of scf12";
          throw InputError(oss.str().c_str(),
                           __FILE__,__LINE__);
        }
      }
    }
  }
  {  // map bs2 to bs12
    const int natom = atom_map2_.size();
    for(int a=0; a<natom; ++a) {
      const int a12 = atom_map2_[a];

      const int nshell = bs2->nshell_on_center(a);
      for(int s=0; s<nshell; ++s) {
        const int S = bs2->shell_on_center(a, s);
        int s12;
        try {
          s12 = ishell_on_center(a12, bs12, bs2->shell(S));
          const int nbf = bs2->shell(S).nfunction();
          const int offset = bs2->shell_to_function(S);
          const int offset12 = bs12->shell_to_function(s12);
          for(int f=0; f<nbf; ++f)
            basis_map2[f + offset] = f + offset12;
        }
        catch (ProgrammingError&) {
          std::ostringstream oss; oss << "shell " << S << " in basis set of scf2 is not found in that of scf12";
          throw InputError(oss.str().c_str(),
                           __FILE__,__LINE__);
        }
      }
    }
  }
  // use basis function maps to map coefficients to the n-mer basis
  RefSCMatrix vec1_12 = vec1.kit()->matrix(vec12.rowdim(), vec1.coldim());  vec1_12.assign(0.0);
  RefSCMatrix vec2_12 = vec2.kit()->matrix(vec12.rowdim(), vec2.coldim());  vec2_12.assign(0.0);
  { // vec1 -> vec1_12
    const int nbasis = vec1.rowdim().n();
    const int nmo = vec1.coldim().n();
    for(int f=0; f<nbasis; ++f) {
      const int f12 = basis_map1[f];
      for(int mo=0; mo<nmo; ++mo)
        vec1_12.set_element(f12, mo, vec1.get_element(f, mo));
    }
  }
  { // vec2 -> vec2_12
    const int nbasis = vec2.rowdim().n();
    const int nmo = vec2.coldim().n();
    for(int f=0; f<nbasis; ++f) {
      const int f12 = basis_map2[f];
      for(int mo=0; mo<nmo; ++mo)
        vec2_12.set_element(f12, mo, vec2.get_element(f, mo));
    }
  }
  vec1 = 0;
  vec2 = 0;

  // Compute how many monomer orbitals to consider for the transfer matrices
  int  nomo_omit1, nomo_omit2;
  int  numo_omit1, numo_omit2;
  const int nocc1 = obwfn1_->nelectron()/2;
  const int nocc2 = obwfn2_->nelectron()/2;
  if (nocc_ == -1) {
    nomo_omit1 = 0;
    nomo_omit2 = 0;
  }
  else {
    nomo_omit1 = nocc1 - nocc_;
    nomo_omit2 = nocc2 - nocc_;
  }
  if (nuocc_ == -1) {
    numo_omit1 = 0;
    numo_omit2 = 0;
  }
  else {
    numo_omit1 = vec1_12.coldim().n() - nocc1 - nuocc_;
    numo_omit2 = vec2_12.coldim().n() - nocc2 - nuocc_;
  }
  if (nomo_omit1 < 0) nomo_omit1 = 0;
  if (nomo_omit2 < 0) nomo_omit2 = 0;
  if (numo_omit1 < 0) numo_omit1 = 0;
  if (numo_omit2 < 0) numo_omit2 = 0;
  // select only the requested HOMOs and LUMOs
  Ref<OrbitalSpace> dspace =  new OrbitalSpace("D", "n-mer basis set space", vec12, obwfn12_->basis(), obwfn12_->integral(), obwfn12_->eigenvalues(), 0, 0);
  Ref<OrbitalSpace> m1space = new OrbitalSpace("m1", "Monomer 1 active MO space", vec1_12, obwfn12_->basis(), obwfn1_->integral(), obwfn1_->eigenvalues(), nomo_omit1, numo_omit1);
  Ref<OrbitalSpace> m2space = new OrbitalSpace("m2", "Monomer 2 active MO space", vec2_12, obwfn12_->basis(), obwfn2_->integral(), obwfn2_->eigenvalues(), nomo_omit2, numo_omit2);

  if (debug_ > 0) {
    RefSCMatrix S11 = compute_overlap_ints(m1space,m1space);
    S11.print("Original S11 matrix");
    Ref<OrbitalSpace> proj11 = gen_project(m1space,m1space,"m1->m1", "Monomer 1 MO space projected on itself",1e-10);
    S11 = compute_overlap_ints(m1space,proj11);
    S11.print("S11 matrix after projection");
  }

  if (debug_ > 0) {
    RefSCMatrix S22 = compute_overlap_ints(m2space,m2space);
    S22.print("Original S22 matrix");
    Ref<OrbitalSpace> proj22 = gen_project(m2space,m2space,"m2->m2", "Monomer 2 MO space projected on itself",1e-10);
    S22 = compute_overlap_ints(m2space,proj22);
    S22.print("S22 matrix after projection");
  }

  RefSCMatrix C1projT = m1space->coefs().t();
  RefSCMatrix C2proj  = m2space->coefs();

  // Compute the n-mer Fock matrix in SO basis -- must transform from MO basis to SO basis
  RefSymmSCMatrix fock12_so(obwfn12_->so_dimension(), obwfn12_->basis_matrixkit());
  fock12_so.assign(0.0);
  RefSCMatrix sobymo = obwfn12_->mo_to_so();
  // SCF Hamiltonian is obtained with effective_fock()
#if 1
  if (obwfn12_clscf.nonnull()) {
    fock12_so.accumulate_transform(obwfn12_->mo_to_so(), obwfn12_clscf->effective_fock());
  }
  // Other (incl. Huckel) Hamiltonian are assumed to be simply the eigenvalues
  else {
#else
  {
#endif
    RefSymmSCMatrix fock12_mo(obwfn12_->oso_dimension(), obwfn12_->basis_matrixkit());
    const int nmo = obwfn12_->oso_dimension().n();
    RefDiagSCMatrix fock12_evals = obwfn12_->eigenvalues();
    fock12_mo.assign(0.0);
    for(int i=0; i<nmo; i++)
      fock12_mo.set_element(i,i,fock12_evals(i));
    fock12_so.accumulate_transform(obwfn12_->mo_to_so(), fock12_mo);
  }

  // now transform the n-mer Fock matrix to AO basis
  RefSymmSCMatrix fock12_ao = plist12->to_AO_basis(fock12_so);

  // Compute n-mer Fock matrix between M1d and M2d
  RefSCMatrix F12 = C1projT * fock12_ao * C2proj;
#if 0 // need to look at other one-electron operators?
  F12 = C1projT * plist12->to_AO_basis(sc::detail::onebodyint<&Integral::hcore>(m1space->basis(), Integral::get_default_integral())) * C2proj;
#endif
  F12.print("Transfer Fock matrix");

  // Compute overlap matrix between M1d and M2d
#if 0
  RefSCMatrix S12 = compute_overlap_ints(proj1,proj2);
#else
  RefSCMatrix S12 = compute_overlap_ints(m1space, m2space);
#endif
  S12.print("Overlap matrix");

  // print out monomer MO energies
  m1space->evals().print("Monomer 1 orbital energies");
  m2space->evals().print("Monomer 2 orbital energies");
  dspace->evals().print("n-mer orbital energies");

  // Compute n-mer Fock matrix in monomer MO bases
  RefSCMatrix F11 = C1projT * fock12_ao * C1projT.t();
  RefSCMatrix F22 = C2proj.t() * fock12_ao * C2proj;
  F11.print("n-mer Fock matrix in monomer 1 basis");
  F22.print("n-mer Fock matrix in monomer 2 basis");

  // print out overlap between monomer and n-mer orbitals to help map monomer HOMOs/LUMOs to the n-mer orbitals
  if (debug_ > 0) {
    RefSCMatrix S1D = compute_overlap_ints(m1space,dspace);  S1D.print("Overlap between monomer 1 and n-mer orbitals");
    RefSCMatrix S2D = compute_overlap_ints(m2space,dspace);  S2D.print("Overlap between monomer 1 and n-mer orbitals");
  }

  //
  // if given ip1 or ip2, use them to recompute adiabatic states
  //
  if (ip1_.empty() == false ||
      ip2_.empty() == false) {

    const int nuocc1 = vec1_12.coldim().n() - nocc1;
    const int nuocc2 = vec2_12.coldim().n() - nocc2;
    Ref<OrbitalSpace> ospace1 = new OrbitalSpace("om1", "Monomer 1 occupied subspace", vec1_12, obwfn12_->basis(), obwfn1_->integral(), obwfn1_->eigenvalues(), 0, nuocc1);
    Ref<OrbitalSpace> ospace2 = new OrbitalSpace("om2", "Monomer 2 occupied subspace", vec2_12, obwfn12_->basis(), obwfn2_->integral(), obwfn2_->eigenvalues(), 0, nuocc2);
    Ref<OrbitalSpace> space1 = new OrbitalSpace("m1", "Monomer 1 subspace", vec1_12, obwfn12_->basis(), obwfn1_->integral(), obwfn1_->eigenvalues(), 0, 0);
    Ref<OrbitalSpace> space2 = new OrbitalSpace("m2", "Monomer 2 subspace", vec2_12, obwfn12_->basis(), obwfn2_->integral(), obwfn2_->eigenvalues(), 0, 0);

    // use 2 methods
    for(int method=0; method<=1; ++method) {
      // method 1: use occupied space only
      // method 2: use all orbitals
      Ref<OrbitalSpace> s1 = (method == 0) ? ospace1 : space1;
      Ref<OrbitalSpace> s2 = (method == 0) ? ospace2 : space2;

      // make union of spaces
      Ref<OrbitalSpace> s12 = new OrbitalSpaceUnion("m1+m2", "Dimer subspace",
                                                    *s1, *s2, true);
      // map smaller to larger spaces
      MOIndexMap map_1_to_12;  map_1_to_12 = *s12 << *s1;  // will throw in case of error
      MOIndexMap map_2_to_12;  map_2_to_12 = *s12 << *s2;
      MOIndexMap map_o1_to_1;  map_o1_to_1 = *s1 << *ospace1;  // will throw in case of error
      MOIndexMap map_o2_to_2;  map_o2_to_2 = *s2 << *ospace2;

      // prepare dimer fock and overlap matrices in s12
      RefSymmSCMatrix F12 = s12->coefs()->kit()->symmmatrix(s12->coefs()->coldim());
      F12.assign(0.0);
      RefSymmSCMatrix S12 = F12.copy();
      RefSCMatrix s12_coefs;
      if (s12->nblocks() == 1) // c1 symmetry? coefficient matrix row dimension is blocked by shells,
                               // but fock12_ao has 1 block (see PetiteList::AO_basisdim())
                               // convert to compatible shape
      {
        s12_coefs = s12->coefs()->kit()->matrix(fock12_ao.dim(), s12->coefs().coldim());
        s12_coefs->convert(s12->coefs());
      }
      else {
        s12_coefs = s12->coefs();
      }
      F12.accumulate_transform(s12_coefs, fock12_ao, SCMatrix::TransposeTransform);
      RefSymmSCMatrix S12_so = sc::detail::overlap(s12->basis(), integral);
      RefSymmSCMatrix S12_ao = plist12->to_AO_basis(S12_so);
      S12.accumulate_transform(s12_coefs, S12_ao, SCMatrix::TransposeTransform);
      RefDiagSCMatrix e12 = s12->evals();  // eigenvalues of *monomer* fock operators (i.e. not including intermonomer polarization)

      //F12.print("F12");
      //S12.print("S12");

      // solve unperturbed eigensystem
      RefDiagSCMatrix evals_0 = F12->kit()->diagmatrix(F12.dim());
      RefSCMatrix evecs_0 = F12->kit()->matrix(F12.dim(), F12.dim());
      try {
        F12->eigensystem(S12.pointer(), evals_0.pointer(), evecs_0.pointer());
      }
      catch (AlgorithmException& e) {
        // failure is likely due to monomer orbitals being strongly nonorthogonal!
        // very likely an error in geometry, or something else
        throw AlgorithmException("generalized eigensolver failed, possible reason could be monomers too close together, etc.",
                                 __FILE__, __LINE__, this->class_desc());
      }

      // modify F using provided IPs
      // for monomer 1
      typedef IPs::iterator iter;
      Ref<Units> eV = new Units("eV");
      const double eV_to_au = eV->to_atomic_units();
      const double au_to_eV = 1.0 / eV_to_au;
      for(iter o1=ip1_.begin(); o1!=ip1_.end(); ++o1) {
        //
        // map occ index to s12
        //
        const unsigned int i1 = map_o1_to_1[ nocc1 - o1->first ];
        const unsigned int i12 = map_1_to_12[ i1 ];
        const double e0 = F12.get_element(i12, i12);
        const double delta_e0 = -o1->second * eV_to_au - e12.get_element(i12);
        const double e1 = e0 + delta_e0;
        F12.set_element(i12, i12, e1);
      }
      // and monomer 2
      for(iter o2=ip2_.begin(); o2!=ip2_.end(); ++o2) {
        // map occ index to s12
        const unsigned int i2 = map_o2_to_2[ nocc2 - o2->first ];
        const unsigned int i12 = map_2_to_12[ i2 ];
        const double e0 = F12.get_element(i12, i12);
        const double delta_e0 = -o2->second * eV_to_au - e12.get_element(i12);
        const double e1 = e0 + delta_e0;
        F12.set_element(i12, i12, e1);
      }

      // and solve perturbed eigensystem
      RefDiagSCMatrix evals_1 = F12->kit()->diagmatrix(F12.dim());
      RefSCMatrix evecs_1 = F12->kit()->matrix(F12.dim(), F12.dim());
      try {
        F12->eigensystem(S12.pointer(), evals_1.pointer(), evecs_1.pointer());

        // test normalization
        //(evecs_1.t() * S12 * evecs_1).print("V^t . S . V: should be one");
      }
      catch (AlgorithmException& e) {
        // failure is likely due to monomer orbitals being strongly nonorthogonal!
        // very likely an error in geometry, or something else
        throw AlgorithmException("generalized eigensolver failed, possible reason could be monomers too close together, etc.",
                                 __FILE__, __LINE__, this->class_desc());
      }

      // print out IPs side-by-side
      ExEnv::out0() << indent << "Ionization potentials (n-mer IPs obtained with "
                    << (method == 0 ? "occ" : "full") << " Fock operator)" << std::endl;
      {
        // sort monomer IPs
        RefDiagSCMatrix evals_m = s12->evals().clone();
        std::vector<double> evals_m_vec(evals_m.n());
        s12->evals().convert(&evals_m_vec[0]);
        std::sort(evals_m_vec.begin(), evals_m_vec.end());
        evals_m.assign(&evals_m_vec[0]);

        ip_print("1-mer", "n-mer", "n-mer+extIP", ExEnv::out0(), 6, -au_to_eV, evals_m, evals_0, evals_1);
      }

      // optional: print adiabatic occupied dimer orbitals on the grid
      if (grid_.nonnull() && method == 0) {
        // compute adiabatic orbitals
        Ref<OrbitalSpace> as12 = new OrbitalSpace("d12", "Dimer active MO space", s12->coefs() * evecs_1,
                                                  s12->basis(), s12->integral(), evals_1, 0, 0);
        // create labels for mos
        std::vector<int> labels;
        for(int o=0; o<as12->rank(); ++o)
          labels.push_back(o + 1);
        Ref<WriteOrbitals> wrtorbs = new WriteOrbitals(as12, labels,
                                                       grid_,
                                                       std::string("gaussian_cube"),
                                                       SCFormIO::fileext_to_filename(".dimer_mo.cube"));
        wrtorbs->run();
      }

    }

  }

  // optional: print monomer orbitals on the grid
  if (grid_.nonnull()) {
    // merge the two spaces together
    Ref<OrbitalSpace> m12space = new OrbitalSpaceUnion("m1+m2", "Monomers 1+2 active MO space",
                                                       *m1space, *m2space, true);
    // create labels for mos: fragment 1 will be negative, fragment 2 positive
    std::vector<int> labels;
    for(int o=0; o<m1space->rank(); ++o)
      labels.push_back(- (o + nomo_omit1 + 1));
    for(int o=0; o<m2space->rank(); ++o)
      labels.push_back(  (o + nomo_omit2 + 1));
    Ref<WriteOrbitals> wrtorbs = new WriteOrbitals(m12space, labels,
                                                   grid_,
                                                   std::string("gaussian_cube"),
                                                   SCFormIO::fileext_to_filename(".mo.cube"));
    wrtorbs->run();
  }

  obwfn12_->print();
  obwfn1_->print();
  obwfn2_->print();
}

void
ETraIn::read_ip(const Ref<KeyVal> & kv, const std::string & ip_key, IPs& ip,
                unsigned int norbs, Ref<OrbitalSpace>& ip_orbs)
{
  Keyword kw(kv, ip_key);
  kw >> ip;
  // validate orbital indices and IP signs
  typedef IPs::iterator iter;
  for(iter i=ip.begin(); i!=ip.end(); ++i) {
    if (i->first < 1 || i->first > norbs) {
      std::ostringstream oss;
      oss << ip_key << ":" << i->first;
      throw InputError("invalid orbital index",
                       __FILE__, __LINE__,
                       oss.str().c_str(), "", this->class_desc());
    }
    if (i->second < 0)
      ExEnv::out0() << indent << "WARNING: negative IP # " << i->first << " in " << ip_key << ", did you mean positive?"<< std::endl;
  }

  // ip_orbs is optional
  const std::string ip_orbs_key = ip_key + "_orbs";
  if (kv->exists(ip_orbs_key.c_str())) {
    Ref<OneBodyWavefunction> ip_orbs_wfn;
    ip_orbs_wfn << kv->describedclassvalue(ip_orbs_key.c_str());
    if (ip_orbs_wfn.null())
      ExEnv::out0() << indent << "WARNING: could not understand value for the optional keyword " << ip_orbs_key << ", will ignore";
    else {
      const int nocc = ip_orbs_wfn->nelectron() / 2;
      MPQC_ASSERT(ip_orbs_wfn->nelectron() % 2 == 0);
      const int norbs = ip_orbs_wfn->so_dimension().n();
      const int nuocc = norbs - nocc;

      RefSCMatrix vec_so = ip_orbs_wfn->eigenvectors();
      Ref<PetiteList> plist = ip_orbs_wfn->integral()->petite_list();
      RefSCMatrix vec = plist->evecs_to_AO_basis(vec_so);

      ip_orbs = new OrbitalSpace(ip_key.c_str(), ip_key.c_str(), vec, ip_orbs_wfn->basis(), ip_orbs_wfn->integral(),
                                 ip_orbs_wfn->eigenvalues(), 0, nuocc);
      MPQC_ASSERT(false); // not implemented yet
    }
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
