//
// soad.cc
//
// Copyright (C) 2013 Edward Valeev
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
#pragma implementation
#endif

#include <cassert>
#include <mpqc_config.h>
#include <chemistry/qc/lcao/soad.h>
#include <chemistry/qc/wfn/femo.h>
#include <chemistry/qc/basis/split.h>
#ifdef HAVE_LIBINT2
#include <chemistry/qc/libint2/libint2.h>
#endif
#include <math/scmat/abstract.h>
#include <chemistry/qc/basis/petite.h>

using namespace sc;

ClassDesc
SuperpositionOfAtomicDensities::class_desc_(typeid(SuperpositionOfAtomicDensities),
                     "SuperpositionOfAtomicDensities",
                     1,               // version
                     "public OneBodyWavefunction", // must match parent
                     0,               // change to create<SuperpositionOfAtomicDensities> if this class is DefaultConstructible
                     create<SuperpositionOfAtomicDensities>, // change to 0 if this class is not KeyValConstructible
                     create<SuperpositionOfAtomicDensities>  // change to 0 if this class is not StateInConstructible
                     );

SuperpositionOfAtomicDensities::SuperpositionOfAtomicDensities(const Ref<KeyVal>& kv) : OneBodyWavefunction(kv)
{
  if (molecule()->max_z() > 86) // for Rb-Rn (Z=86) will use WTBS, for H-Kr will use STO-6G
    throw FeatureNotImplemented("SuperpositionOfAtomicDensities only implemented for up to Z=86",
                                __FILE__,__LINE__,this->class_desc());

  relax_ = kv->booleanvalue("relax", KeyValValueboolean(true));
}

SuperpositionOfAtomicDensities::SuperpositionOfAtomicDensities(StateIn& si) : OneBodyWavefunction(si) {
  si.get(relax_);
}

SuperpositionOfAtomicDensities::~SuperpositionOfAtomicDensities() {
  // this may be necessary if this is a templated class
  const bool make_sure_class_desc_initialized = (&class_desc_ != 0);
}

void
SuperpositionOfAtomicDensities::save_data_state(StateOut& so) {
  so.put(relax_);
}

RefSymmSCMatrix sc::SuperpositionOfAtomicDensities::density() {
  if (not density_.computed())
    compute();
  return density_;
}

double sc::SuperpositionOfAtomicDensities::magnetic_moment() const {
  // count unpaired electrons
  const int norbs = occs_.n();
  double magmom = 0.0;
  for(int o=0; o<norbs; ++o)
    if (occs_.get_element(o) == 1.0)
      magmom += 1.0;
  return magmom;
}

int sc::SuperpositionOfAtomicDensities::nelectron() {
  return this->molecule()->total_Z();
}

int sc::SuperpositionOfAtomicDensities::value_implemented() const {
  return 1;
}

void sc::SuperpositionOfAtomicDensities::compute() {

  Ref<Integral> intf = integral()->clone();
  intf->set_storage(integral()->storage_used() + integral()->storage_unused());

    if (not relax_) {
      // get guess AO density in this basis
      // and transform it to SO basis
      RefSymmSCMatrix density_ao =
          SuperpositionOfAtomicDensities::guess_density(basis(),
                                                        intf);
      RefSymmSCMatrix dens = basis_matrixkit()->symmmatrix(so_dimension());
      {
        Ref<PetiteList> pl = integral()->petite_list();
        // copy density_ao into a blocked matrix of the desired dimensions
        RefSymmSCMatrix density_ao_blk = basis_matrixkit()->symmmatrix(
            pl->AO_basisdim());
        density_ao_blk->convert(density_ao);

        dens->assign(0.0);
        dens->accumulate_transform(pl->sotoao(), density_ao_blk);
      }
      if (debug_ > 1)
        dens->print("projected SO density");
      density_ = dens;
    }
    else { // relax_ == true

      // this will take quite a bit of work
      // 0) make minimal basis
      // 1) make a wave function based on a minbasis
      // 2) make a WavefunctionWorld around this wfn; optionally tell it to use density fitting, and such
      // 3) register basis()
      // 4) use world to get a Fock matrix
      // 5) diagonalize
      // 6) repopulate (using FEMO?)

      // make minimal basis
      Ref<GaussianBasisSet> minbasis = SuperpositionOfAtomicDensities::minimal_basis_set(molecule());

      // make a SOAD based on minbasis
      // this is purely for the sake of serving as the leader of the world since
      // this wavefunction's basis() is not the basis in which we will have the density matrix
      Ref<Wavefunction> minbasis_wfn;
      {
        Ref<AssignedKeyVal> akv = new AssignedKeyVal;
        akv->assign("molecule", Ref<DescribedClass>(minbasis->molecule()));
        akv->assign("basis", Ref<DescribedClass>(minbasis));
        akv->assign("relax", false);
        minbasis_wfn = new SuperpositionOfAtomicDensities(akv);
      }

      // make WavefunctionWorld
      Ref<WavefunctionWorld> world;
      {
        Ref<AssignedKeyVal> akv = new AssignedKeyVal;
        akv->assign("wfn", Ref<DescribedClass>(minbasis_wfn));
        akv->assign("ints_precision", 1e-4);
        world = new WavefunctionWorld(akv);
      }

      // we have to use 2 (or 3) basis sets here, make sure the world knows about them
      Ref<OrbitalSpaceRegistry> oreg = world->moints_runtime()->factory()->orbital_registry();
      Ref<AOSpaceRegistry> aoreg = world->moints_runtime()->factory()->ao_registry();
      { // register the minimal basis
        std::string aospace_id = new_unique_key(oreg);
        if (aoreg->key_exists(minbasis) == false) {
          Ref<OrbitalSpace> aospace = new AtomicOrbitalSpace(
              aospace_id, "SuperpositionOfAtomicDensities minimal AO basis set",
              minbasis, integral());
          aoreg->add(minbasis, aospace);
          MPQC_ASSERT(oreg->key_exists(aospace_id) == false);
          // should be ensured by using new_unique_key
          oreg->add(make_keyspace_pair(aospace));
        }
      }
      { // register the full basis
        std::string aospace_id = new_unique_key(oreg);
        if (aoreg->key_exists(basis()) == false) {
          Ref<OrbitalSpace> aospace = new AtomicOrbitalSpace(
              aospace_id, "SuperpositionOfAtomicDensities AO basis set",
              basis(), integral());
          aoreg->add(basis(), aospace);
          MPQC_ASSERT(oreg->key_exists(aospace_id) == false);
          // should be ensured by using new_unique_key
          oreg->add(make_keyspace_pair(aospace));
        }
      }

      // get guess AO density in the minimal basis and feed it to the world's fock builder
      RefSymmSCMatrix density_ao =
          SuperpositionOfAtomicDensities::guess_minimal_density(minbasis,
                                                                intf);
      density_ao->scale(0.5); // assuming density spin averaged
      // copy density_ao into a blocked matrix of the desired dimensions
      RefSymmSCMatrix density_ao_blk = basis_matrixkit()->symmmatrix(density_ao.dim());
      density_ao_blk->convert(density_ao);
      world->fockbuild_runtime()->set_densities(density_ao_blk, density_ao_blk);

      // now fun part: build the Fock matrix in the given basis with the density expressed in the minimal basis
      Ref<OrbitalSpace> aospace = aoreg->value(basis());
      RefSCMatrix F;
      {
        const std::string fkey = ParsedOneBodyIntKey::key(aospace->id(),aospace->id(),std::string("F"));
        RefSCMatrix Floc = world->fockbuild_runtime()->get(fkey);
        F = Floc.copy();
      }
#if 0
      {
        const std::string jkey = ParsedOneBodyIntKey::key(aospace->id(),aospace->id(),std::string("J"));
        RefSCMatrix J = world->fockbuild_runtime()->get(jkey);
        G = J.copy();
      }
      {
        const std::string kkey = ParsedOneBodyIntKey::key(aospace->id(),aospace->id(),std::string("K"),AnySpinCase1);
        RefSCMatrix K = world->fockbuild_runtime()->get(kkey);
        G.accumulate( -1.0 * K);
      }
#endif
      Ref<SCElementOp> accum_F_op = new SCElementAccumulateSCMatrix(F.pointer());
      RefSymmSCMatrix F_symm = F.kit()->symmmatrix(F.coldim()); F_symm.assign(0.0);
      F_symm.element_op(accum_F_op); F = 0;
      Ref<PetiteList> pl = integral()->petite_list();
      F_symm = pl->to_SO_basis(F_symm);

      // transform to OSO basis and diagonalize the Fock matrix
      RefSCMatrix ortho = so_to_orthog_so();
      RefSymmSCMatrix F_oso(oso_dimension(), basis_matrixkit());
      F_oso.assign(0.0);
      F_oso.accumulate_transform(ortho,F_symm);
      RefDiagSCMatrix evals = F_oso.kit()->diagmatrix(F_oso.dim());
      RefSCMatrix evecs = F_oso.kit()->matrix(F_oso.dim(),F_oso.dim());
      F_oso.diagonalize(evals, evecs);
      eigenvalues_ = evals;  eigenvalues_.computed() = 1;
      oso_eigenvectors_ = evecs;  oso_eigenvectors_.computed() = 1;

      // compute occupations using Hunds
      //                            # electron,       Etol, allow closed-shell?
      HundsFEMOSeeker femoseeker(nelectron(), HundsFEMOSeeker::tolerance, true,
                                 evals, evals);
      Ref<FEMO> femo = femoseeker.result();
      occs_ = evals.clone();

      BlockedDiagSCMatrix *occs_blk =
          require_dynamic_cast<BlockedDiagSCMatrix*>(occs_.pointer(),
              "SuperpositionOfAtomicDensities::compute: val"
              );
      for(int i=0; i<occs_blk->nblocks(); ++i) {
        if (occs_blk->block(i)) {
          const int nso_blk = occs_blk->block(i).n();
          const int nalpha = femo->nalpha(i);
          const int nbeta = femo->nbeta(i);
          for(int o=0; o<nso_blk; ++o) {
            double occ = 0.0;
            if (o < nalpha) occ += 1.0;
            if (o < nbeta) occ += 1.0;
            occs_blk->block(i)->set_element(o, occ);
          }
        }
      }

      // transform oso_eigenvector to so_eigenvector
      RefSCMatrix evecs_so = so_to_orthog_so().t() * evecs;

      // form the SO density
      RefSymmSCMatrix density_so = F_oso.kit()->symmmatrix(evecs_so.rowdim());
      density_so.assign(0.0);
      density_so.accumulate_transform(evecs_so, occs_);
      if (debug_ > 1)
        density_so->print("relaxed SO density");
      density_ = density_so;
    }
    density_.computed() = 1;

    // test idempotency of this density and McWeeny-purify if necessary
#if 0
    {
      RefSymmSCMatrix S = overlap(); // in SO basis
      RefSymmSCMatrix R = density_->copy();
      R.scale(0.5);

      // Purify an approximate density matrix so that it meets the requirements
      // of a density matirx mainly that
      // R = R^T
      //  Trace(S.R) = Number of electrons * 0.5
      // All eigenvalues of the matrix are 1.

      ExEnv::out0() << indent << "Purifying SuperpositionOfAtomicDensities\n";

      // This routine is essentially R = 3*R^2 - 2 * R^3 it only converges if R is
      // close enough to the correct answer.
      double rdiff = 0.0;
      do {

        { // compute natural occupation numbers and natural orbitals

          // transform the density into an orthogonal basis
          RefSCMatrix ortho = so_to_orthog_so();

          RefSymmSCMatrix densortho(oso_dimension(), basis_matrixkit());
          densortho.assign(0.0);
          densortho.accumulate_transform(so_to_orthog_so_inverse().t(),R);
          densortho.eigvals().print("natural occupancies");
          //densortho.eigvecs().print("natural orbitals");

        }

        RefSCMatrix RS = R * S;

        /// Peform the actual iteration step.
        RefSCMatrix RSR = RS * R;
        RefSCMatrix Rnew_sq = 3.0 * RSR - 2 * RS * RSR;
        RefSymmSCMatrix Rnew = R.clone();
        Rnew.assign(0.0);
        Rnew.accumulate_symmetric_sum(Rnew_sq);
        Rnew.scale(0.5);

        // norm2 of the lower triangle!
        Ref<SCElementKNorm> norm2_op = new SCElementKNorm(2);
        (Rnew - R).element_op(norm2_op);
        rdiff = norm2_op->result();
        ExEnv::out0() << indent << "res = " << scprintf("%e",rdiff) << std::endl;
        R = Rnew;

      }while(rdiff > 1e-8);

      density_ = 2.0 * R;
    }
#endif
}

void sc::SuperpositionOfAtomicDensities::obsolete() {
  OneBodyWavefunction::obsolete();
}

RefSCMatrix sc::SuperpositionOfAtomicDensities::oso_eigenvectors() {
  if (not relax_) {
    if (!natural_orbitals_.computed()) {
      compute();
      density_.result_noupdate()->scale(-1.0); // little trickeration to avoid recomputation of density in natural_orbitals()
      this->natural_orbitals(); // this will set *negatives* of occupation numbers as eigenvalues to order them properly
      density_.result_noupdate()->scale(-1.0);

      oso_eigenvectors_ = so_to_orthog_so_inverse() * natural_orbitals();
      oso_eigenvectors_.computed() = 1;
      //oso_eigenvectors_->print("SOAD OSO eigenvectors");
    }
  }
  else {
    if (not oso_eigenvectors_.computed())
      compute();
  }

  return oso_eigenvectors_.result_noupdate();
}

RefDiagSCMatrix sc::SuperpositionOfAtomicDensities::eigenvalues() {
  if (not relax_) {
    if (natural_density_.computed() == false) {
      this->oso_eigenvectors();
    }
    return natural_density();
  }
  else {
    if (not eigenvalues_.computed())
      compute();
    return eigenvalues_.result_noupdate();
  }
}

Ref<GaussianBasisSet>
SuperpositionOfAtomicDensities::minimal_basis_set(const Ref<Molecule>& mol) {
  // make mother minimal basis
  Ref<GaussianBasisSet> mother;
  try {
    Ref<AssignedKeyVal> akv = new AssignedKeyVal;
    akv->assign("molecule", Ref<DescribedClass>(mol));
    // mix STO-6G (for H-Kr) and WTB (for Rb-Rn)
    for (int a = 0; a < mol->natom(); ++a) {
      std::ostringstream oss;
      oss << "basis:" << a;
      const std::string keyword = oss.str();
      if (mol->Z(a) <= 38)
        akv->assign(keyword.c_str(), "STO-3G");
      else
        akv->assign(keyword.c_str(), "WTBS");
    }
    mother = new GaussianBasisSet(akv);
  }
  catch(...) {}
  if (mother.null())
    throw ProgrammingError("could not construct a minimal basis for SuperpositionOfAtomicDensities",
                           __FILE__, __LINE__);

  // and split (because some integrals need it, and because the logic of atomic aufbau does not tolerate sp shells)
  Ref<AssignedKeyVal> akv1 = new AssignedKeyVal;
  akv1->assign("molecule", Ref<DescribedClass>(mol));
  akv1->assign("basis", Ref<DescribedClass>(mother));
  Ref<GaussianBasisSet> result = new SplitBasisSet(Ref<KeyVal>(akv1));

  return result;
}

RefSymmSCMatrix
SuperpositionOfAtomicDensities::guess_minimal_density(const Ref<GaussianBasisSet>& minimal_basis_set,
                                                      const Ref<Integral>& intf) {

  Ref<Molecule> mol = minimal_basis_set->molecule();
  intf->set_basis(minimal_basis_set);
  Ref<PetiteList> pl = intf->petite_list();
  RefSymmSCMatrix minbasis_density = minimal_basis_set->matrixkit()->symmmatrix(pl->AO_basisdim());

  // populate atomic orbitals with electrons, smear density across partially filled subshells
  minbasis_density.assign(0.0);

  const int ncenters = mol->natom();
  for (int center = 0; center < ncenters; ++center) {
    // using charge() instead of Z() to account for ghost atoms
    // skipping point charges also
    if (mol->charge(center) != 0.0 && not mol->is_Q(center)) {
      const int Z = mol->Z(center);
      int nelectrons = Z;

      std::map<int, std::vector<int> > shell_occs; // maps l -> occupations of subshells (subshell of chemistry = shell of QC),
      // in aufbau order
      int aufbau_l_order[] = { 0, // 1s
        0, // 2s
        1, // 2p
        0, // 3s
        1, // 3p
        0, // 4s
        2, // 3d
        1, // 4p
        0, // 5s
        2, // 4d
        1, // 5p
        0, // 6s
        3, // 4f
        2, // 5d
        1, // 6p
        0, // 7s
        3, // 5f
        2, // 6d
        1 // 7p
        };
      const int aufbau_size = sizeof(aufbau_l_order) / sizeof(int);

      int k = 0;
      while (nelectrons > 0) {
        const int nelectrons_to_subshell = std::min(nelectrons,
                                                    4 * aufbau_l_order[k] + 2);
        shell_occs[aufbau_l_order[k]].push_back(nelectrons_to_subshell);
        nelectrons -= nelectrons_to_subshell;
        ++k;
      }

      // sort shells to groups of l blocks
      // since can't assume that they appear in aufbau or even l order;
      // however we do assume that within same-l set they do appear in aufbau order!
      std::map<int, std::vector<int> > l_to_shells;
      const int nshells = minimal_basis_set->nshell_on_center(center);
      for (int s = 0; s < nshells; s++) {
        const int ss = minimal_basis_set->shell_on_center(center, s);
        // this is a SplitBasisSet, hence Shells have 1 l
        l_to_shells[minimal_basis_set->shell(ss).max_am()].push_back(ss);
      }

      // fill up the density:
      // 1) for each l
      typedef std::map<int, std::vector<int> >::const_iterator citer;
      for (citer i = l_to_shells.begin(); i != l_to_shells.end(); ++i) {
        const int l = i->first;
        const std::vector<int>& shells = i->second;
        // loop over shells
        const int nocc_shells_of_this_l = shell_occs[l].size();
        // if this fails, aufbau produced more occupied shells than in the basis
        // either aufbau algorithm is broken or the basis is broken
        MPQC_ASSERT(nocc_shells_of_this_l <= shells.size());
        for (int s = 0; s < nocc_shells_of_this_l; ++s) {
          const int nelectrons_in_shell = shell_occs[l].at(s);
          const double nelectrons_per_bf =
              static_cast<double>(nelectrons_in_shell) / (2 * l + 1);

          // place electrons into each function
          const int shell = shells[s];
          const int nbf = minimal_basis_set->shell(shell).nfunction();
          for (int bf = 0; bf < nbf; ++bf) {
            const int ao = minimal_basis_set->shell_to_function(shell) + bf;
            minbasis_density(ao, ao) = nelectrons_per_bf;
          }
        }
      }

    }
  }
  //minbasis_density.print("minimal basis AO guess density");

  return minbasis_density;
}

RefSymmSCMatrix SuperpositionOfAtomicDensities::guess_density(
    const Ref<GaussianBasisSet>& basis,
    const Ref<Integral>& intf) {

  // generate density in minimal basis
  Ref<GaussianBasisSet> minbasis = SuperpositionOfAtomicDensities::minimal_basis_set(basis->molecule());
  RefSymmSCMatrix minbasis_density = SuperpositionOfAtomicDensities::guess_minimal_density(minbasis, intf);

  // make full basis density by projecting minbasis density

  // make Smg = min AO x given AO overlap
  intf->set_basis(minbasis, basis);
  RefSCMatrix Smg(minbasis_density.dim(), basis->basisdim(),
                  basis->matrixkit());
  {
    Ref<SCElementOp> op = new OneBodyIntOp(intf->overlap());
    Smg.assign(0.0);
    Smg.element_op(op);
  }
  //Smg.print("(min|given) AO overlap");

  // compute given AO overlap, and its inverse
  RefSymmSCMatrix Sgg(basis->basisdim(),
                      basis->matrixkit());
  intf->set_basis(basis, basis);
  {
    Ref<SCElementOp> op = new OneBodyIntOp(intf->overlap());
    Sgg.assign(0.0);
    Sgg.element_op(op);
  }
  RefSymmSCMatrix Sgg_inv = Sgg.gi(1e12);

  // P(given) = U^t P(min) U, where U = Smg Sgg^-1
  RefSCMatrix U = Smg * Sgg_inv;
  RefSymmSCMatrix density_ao = minbasis_density.kit()->symmmatrix(basis->basisdim());
  density_ao.assign(0.0);
  density_ao.accumulate_transform(U, minbasis_density, SCMatrix::TransposeTransform);
  //density_ao.print("projected AO density");
  //ExEnv::out0() << indent << "trace of the projected AO density = " << (density_ao * Sgg).trace() << std::endl;

  return density_ao;
}

int sc::SuperpositionOfAtomicDensities::spin_unrestricted() {
  return 0;
}

double sc::SuperpositionOfAtomicDensities::occupation(int irrep,
                                                      int vectornum) {
  if (!eigenvalues_.computed())
    compute();
  RefDiagSCMatrix occs = relax_ ? occs_ : eigenvalues();
  BlockedDiagSCMatrix *occs_blk =
      require_dynamic_cast<BlockedDiagSCMatrix*>(occs.pointer(),
          "OneBodyWavefunction::projected_eigenvalues: val"
          );
  MPQC_ASSERT(irrep >= 0 && irrep < occs_blk->nblocks());
  double occ = 0.0;
  if (occs_blk->block(irrep)) {
    if (vectornum < occs_blk->block(irrep)->n()) {
      occ = occs_blk->block(irrep)->get_element(vectornum);
      if (not relax_) occ *= -1.0; // eigenvalues = -occupancies
    }
  }
  return occ;
}



/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
