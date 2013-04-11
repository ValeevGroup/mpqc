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

#include <chemistry/qc/wfn/soad.h>
#include <chemistry/qc/basis/split.h>
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
  total_charge_ = kv->intvalue("total_charge", KeyValValueint(0));
  assert(total_charge_ == 0); // nonzero charge not implemented yet
  spin_unrestricted_ = kv->booleanvalue("spin_unrestricted", KeyValValueboolean(false));
  assert(spin_unrestricted_ == false); // unrestricted not implemented yet
}

SuperpositionOfAtomicDensities::SuperpositionOfAtomicDensities(StateIn& si) : OneBodyWavefunction(si) {
  si.get(total_charge_);
}

SuperpositionOfAtomicDensities::~SuperpositionOfAtomicDensities() {
  // this may be necessary if this is a templated class
  const bool make_sure_class_desc_initialized = (&class_desc_ != 0);
}

void
SuperpositionOfAtomicDensities::save_data_state(StateOut& so) {
  so.put(total_charge_);
}

RefSymmSCMatrix sc::SuperpositionOfAtomicDensities::density() {
  compute();
  return density_;
}

int sc::SuperpositionOfAtomicDensities::spin_polarized() {
  return this->nelectron() % 2;
}

int sc::SuperpositionOfAtomicDensities::nelectron() {
  return static_cast<int>(this->molecule()->nuclear_charge()) - total_charge_;
}

int sc::SuperpositionOfAtomicDensities::value_implemented() const {
  return 1;
}

void sc::SuperpositionOfAtomicDensities::compute() {
  if (minbasis_.null()) {

    { // make minimal basis

      // make mother STO-6G basis
      Ref<AssignedKeyVal> akv = new AssignedKeyVal;
      akv->assign("molecule", Ref<DescribedClass>(molecule()));
      akv->assign("name", "STO-6G");
      Ref<GaussianBasisSet> mother = new GaussianBasisSet(akv);

      // and split
      Ref<AssignedKeyVal> akv1 = new AssignedKeyVal;
      akv1->assign("molecule", Ref<DescribedClass>(molecule()));
      akv1->assign("basis", Ref<DescribedClass>(mother));
      minbasis_ = new SplitBasisSet(Ref<KeyVal>(akv1));
    }

    // make minimal basis density
    Ref<Integral> localintfactory = integral()->clone();
    localintfactory->set_basis(minbasis_);
    Ref<PetiteList> pl = localintfactory->petite_list();
    Ref<PetiteList> plg = integral()->petite_list();
    minbasis_density_ = minbasis_->matrixkit()->symmmatrix(pl->AO_basisdim());

    { // populate atomic orbitals with electrons, smear density across partially filled subshells
      minbasis_density_.assign(0.0);

      const int ncenters = molecule()->natom();
      for(int center=0; center<ncenters; ++center) {
        if (molecule()->charge(center) != 0.0) { // using charge() instead of Z() to account for ghost atoms
          const int Z = molecule()->Z(center);
          int nelectrons = Z;

          std::map<int, std::vector<int> > shell_occs; // maps l -> occupations of subshells (subshell of chemistry = shell of QC),
                                                 // in aufbau order
          int aufbau_l_order[] = {
            0, // 1s
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
            1  // 7p
          };
          const int aufbau_size = sizeof(aufbau_l_order)/sizeof(int);

          int k=0;
          while (nelectrons > 0) {
            const int nelectrons_to_subshell = std::min(nelectrons, 4*aufbau_l_order[k] + 2);
            shell_occs[aufbau_l_order[k]].push_back(nelectrons_to_subshell);
            nelectrons -= nelectrons_to_subshell;
            ++k;
          }

          // sort shells to groups of l blocks
          // since can't assume that they appear in aufbau or even l order;
          // however we do assume that within same-l set they do appear in aufbau order!
          std::map<int, std::vector<int> > l_to_shells;
          const int nshells = minbasis_->nshell_on_center(center);
          for(int s=0; s<nshells; s++) {
            const int ss = minbasis_->shell_on_center(center, s);
            // this is a SplitBasisSet, hence Shells have 1 l
            l_to_shells[minbasis_->shell(ss).max_am()].push_back(ss);
          }

          // fill up the density:
          // 1) for each l
          typedef std::map<int, std::vector<int> >::const_iterator citer;
          for(citer i=l_to_shells.begin();
              i!=l_to_shells.end();
              ++i) {
            const int l = i->first;
            const std::vector<int>& shells = i->second;
            // loop over shells
            for(int s=0; s<shells.size(); ++s) {
              const int nelectrons_in_shell = shell_occs[l].at(s);
              const double nelectrons_per_bf = static_cast<double>(nelectrons_in_shell) / (2*l + 1);

              // place electrons into each function
              const int shell = shells[s];
              const int nbf = minbasis_->shell(shell).nfunction();
              for(int bf=0; bf<nbf; ++bf) {
                const int ao = minbasis_->shell_to_function(shell) + bf;
                minbasis_density_(ao,ao) = nelectrons_per_bf;
              }

            }
          }

        }
      }
    }

    //minbasis_density_.print("minimal basis AO guess density");

    { // make full basis density by projecting minbasis density

      // make Smg = min AO x given AO overlap
      localintfactory->set_basis(minbasis_, basis());
      RefSCMatrix Smg(minbasis_density_.dim(), basis()->basisdim(),
                      basis()->matrixkit());
      {
        Ref<SCElementOp> op = new OneBodyIntOp(localintfactory->overlap());
        Smg.assign(0.0);
        Smg.element_op(op);
      }
      if (debug_ > 1) Smg.print("(min|given) AO overlap");

      // compute given AO overlap, and its inverse
      RefSymmSCMatrix Sgg(basis()->basisdim(),
                          basis()->matrixkit());
      {
        Ref<SCElementOp> op = new OneBodyIntOp(integral()->overlap());
        Sgg.assign(0.0);
        Sgg.element_op(op);
      }
      RefSymmSCMatrix Sgg_inv = Sgg.gi(1e12);

      // P(given) = U^t P(min) U, where U = Smg Sgg^-1
      RefSCMatrix U = Smg * Sgg_inv;
      RefSymmSCMatrix density_ao = matrixkit()->symmmatrix(basis()->basisdim());
      density_ao.assign(0.0);
      density_ao.accumulate_transform(U, minbasis_density_, SCMatrix::TransposeTransform);
      if (debug_ > 1) density_ao.print("projected AO density");
      if (debug_ > 0)
        ExEnv::out0() << indent << "trace of the projected AO density = " << (density_ao * Sgg).trace() << std::endl;

      // tranform P(given) from AO to SO basis
      RefSymmSCMatrix dens = basis_matrixkit()->symmmatrix(so_dimension());
      {
        // copy density_ao into a blocked matrix of the desired dimensions
        RefSymmSCMatrix density_ao_blk = basis_matrixkit()->symmmatrix(plg->AO_basisdim());
        density_ao_blk->convert(density_ao);

        dens->assign(0.0);
        dens->accumulate_transform(plg->sotoao(), density_ao_blk);
      }
      if (debug_ > 1) dens->print("projected SO density");
      density_ = dens;
      density_.computed() = 1;

    }

  }
}

void sc::SuperpositionOfAtomicDensities::obsolete() {
  OneBodyWavefunction::obsolete();
  minbasis_ = 0;
  minbasis_density_ = 0;
}

int sc::SuperpositionOfAtomicDensities::spin_unrestricted() {
  return spin_unrestricted_ ? 1 : 0;
}

RefSCMatrix sc::SuperpositionOfAtomicDensities::oso_eigenvectors() {
  if (!natural_orbitals_.computed()) {
    density();
    density_->scale(-1.0);
    this->natural_orbitals(); // this will set
    density_->scale(-1.0);

    oso_eigenvectors_ = so_to_orthog_so_inverse() * natural_orbitals();
    oso_eigenvectors_.computed() = 1;
  }
  return oso_eigenvectors_.result_noupdate();
}

RefDiagSCMatrix sc::SuperpositionOfAtomicDensities::eigenvalues() {
  if (natural_density_.computed() == false) {
    this->oso_eigenvectors();
  }
  return natural_density();
}

double sc::SuperpositionOfAtomicDensities::occupation(int irrep,
                                                      int vectornum) {
  // make sure natural occupancies are computed
  RefDiagSCMatrix occs = eigenvalues();
  BlockedDiagSCMatrix *occs_blk =
      require_dynamic_cast<BlockedDiagSCMatrix*>(occs.pointer(),
          "OneBodyWavefunction::projected_eigenvalues: val"
          );
  assert(irrep >= 0 && irrep < occs_blk->nblocks());
  double occ = 0.0;
  if (occs_blk->block(irrep).nonnull()) {
    if (vectornum < occs_blk->block(irrep)->n())
      occ = -1.0 * occs_blk->block(irrep)->get_element(vectornum);  // eigenvalues = -occupancies
  }
  return occ;
}



/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
