
/*
 * Copyright 2009 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This file is a part of the MPQC LMP2 library.
 *
 * The MPQC LMP2 library is free software: you can redistribute it
 * and/or modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 */

// Pipek-Mezey localization of the occupied molecular orbitals
// Reference: J. Pipek and P. G. Mezey, JCP v.90, p.4916 (1989)

#include <util/misc/regtime.h>
#include <math/scmat/matrix.h>
#include <math/scmat/vector3.h>
#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/molecule/energy.h>

#include <set>
#include <map>
#include <vector>
#include <stdexcept>

#include <chemistry/qc/lmp2/pop_local.h>

using namespace std;

namespace sc {

RefSCMatrix
convert_complete_to_occupied_vector(const Ref<OneBodyWavefunction> &wfn,
                                    const RefSCMatrix &vec)
{
  
  const double occthresh = 1.99999;

  // Compute occ orb offset array (since symmetry causes occupied orbitals
  // to be spread out in vector).
  int nblock = wfn->so_dimension()->blocks()->nblock();
  int *n_in_block = new int[nblock];
  int nocc = 0;
  for (int i=0; i<nblock; i++) {
    n_in_block[i] = 0;
    for (int j=0; j<wfn->so_dimension()->blocks()->size(i); j++) {
      if (wfn->occupation(i,j) > occthresh) {
        n_in_block[i]++;
        nocc++;
        }
      }
    }

  // Construct the dimensions of the occupied vector
  Ref<SCBlockInfo> blockinfo = new SCBlockInfo(nocc, nblock, n_in_block);
  for (int i=0; i<nblock; i++) {
      blockinfo->set_subdim(i, new SCDimension(n_in_block[i]));
    }
  RefSCDimension coldim = new SCDimension(blockinfo,"occ");
  RefSCDimension rowdim = vec->rowdim();

  RefSCMatrix occvec(rowdim,coldim,vec.kit());
  for (int i=0; i<nblock; i++) {
    int norb_in_block = wfn->so_dimension()->blocks()->size(i);
    int iocc_in_block = 0;
    int iorb_in_block = 0;
    RefSCMatrix occvecb
        = dynamic_cast<BlockedSCMatrix*>(occvec.pointer())->block(i);
    if (occvecb == 0) throw runtime_error("vector not blocked");
    RefSCMatrix vecb
        = dynamic_cast<BlockedSCMatrix*>(vec.pointer())->block(i);
    if (occvecb == 0) throw runtime_error("new vector not blocked");
    for (int j=0; j<norb_in_block; j++) {
      if (wfn->occupation(i,j) > occthresh) {
          occvecb.assign_column(vecb.get_column(iorb_in_block),iocc_in_block);
          iocc_in_block++;
        }
      iorb_in_block++;
      }
    }

  delete[] n_in_block;

  return occvec;
}

RefSCMatrix
convert_complete_to_occupied_vector_nosymm(const Ref<OneBodyWavefunction> &wfn,
                                           int nfzc, const RefSCMatrix &vec)
{
  // Use this for C1 symmetry only - this works also with frozen
  // core orbitals. 
  // Use this for now until we get the routine
  // convert_complete_to_occupied_vector to work with frozen core orbitals

  int nocc = wfn->nelectron()/2;

  RefSCDimension coldim = new SCDimension(nocc-nfzc);
  coldim->blocks()->set_subdim(0, new SCDimension(nocc-nfzc));
  RefSCDimension rowdim = vec.rowdim();

  RefSCMatrix occvec(rowdim,coldim,vec.kit());
  for (int i=0; i<rowdim.n(); i++) {
      for (int j=0; j<nocc-nfzc; j++) {
          occvec(i,j) = vec(i,j+nfzc);
        }
    }

  return occvec;
}

// Returns scf_vector (AO->MO transformation matrix, rows:AO's, columns:MO's)
static RefSCMatrix
ao_to_mo(const Ref<OneBodyWavefunction> &wfn)
{
  return (wfn->so_to_mo()*(wfn->integral()->petite_list()->aotoso()).t()).t();
}

// Returns scf_vector, including only active occupied MO's
// (AO->MO transformation matrix, rows:AO's, columns:active occupied MO's)
static RefSCMatrix
ao_to_occact_mo(const Ref<OneBodyWavefunction> &wfn, int nfzc)
{
  RefSCMatrix scf_vector = (wfn->so_to_mo()*(wfn->integral()->petite_list()->aotoso()).t()).t();
  return convert_complete_to_occupied_vector_nosymm(wfn, nfzc, scf_vector);
}

class RotInfo {
  public:
    int mo1;
    int mo2;
    double delta;
    double gamma;
    RotInfo(): mo1(0), mo2(0), delta(0), gamma(0) {}
    RotInfo(int mo1a,int mo2a,double deltaa,double gammaa):
      mo1(mo1a), mo2(mo2a), delta(deltaa), gamma(gammaa) {}
};

class DeltaGreater {
  public:
    bool operator() (const RotInfo &r1, const RotInfo &r2) {
      return r1.delta > r2.delta;
    }
};

void
compute_delta(int natoms, const Ref<GaussianBasisSet> &basis, int nocc_act,
              double *scf_vector_dat, double *S_half_trans_dat,
              int mo1, int mo2, double *delta, double *gamma)
{
  // First compute the matrix elements <mo1|P_atom|mo2>
  // of the atomic population operator P_atom and the
  // intermediates A_mo1_mo2, B_mo1_mo2
  double A_mo1_mo2 = 0.0;
  double B_mo1_mo2 = 0.0;
  int offset = 0;
  int nbasis = basis->nbasis();
  for (int atom=0; atom<natoms; atom++) {
      double P_mo1_mo2 = 0.0;
      double P_mo1_mo1 = 0.0;
      double P_mo2_mo2 = 0.0;
      double * RESTRICT tmp_ptr1 = &scf_vector_dat[offset + mo1*nbasis];
      double * RESTRICT tmp_ptr2 = &scf_vector_dat[offset + mo2*nbasis];
      double * RESTRICT tmp_ptr3 = &S_half_trans_dat[offset + mo1*nbasis];
      double * RESTRICT tmp_ptr4 = &S_half_trans_dat[offset + mo2*nbasis];
      int index = basis->nbasis_on_center(atom);
      for (int ao1=0; ao1<index; ao1++) {
          double tmp1 = *tmp_ptr1;
          double tmp2 = *tmp_ptr2;
          double tmp3 = *tmp_ptr3;
          double tmp4 = *tmp_ptr4;
          P_mo1_mo2 += tmp1*tmp4 + tmp2*tmp3;
          P_mo1_mo1 += tmp1*tmp3;
          P_mo2_mo2 += tmp2*tmp4;
          tmp_ptr1 ++;
          tmp_ptr2 ++;
          tmp_ptr3 ++;
          tmp_ptr4 ++;
        }
      offset += index;
      P_mo1_mo2 *= 0.5;
      // Compute intermediates A_mo1_mo2, B_mo1_mo2
      A_mo1_mo2 += P_mo1_mo2*P_mo1_mo2 - 0.25*(P_mo1_mo1 - P_mo2_mo2)*(P_mo1_mo1 - P_mo2_mo2);
      B_mo1_mo2 += P_mo1_mo2*(P_mo1_mo1 - P_mo2_mo2);
    }
        
  // Compute rotation angle gamma and the corresponding change, delta,
  // in the functional
  double tmp = sqrt(A_mo1_mo2*A_mo1_mo2 + B_mo1_mo2*B_mo1_mo2);
  if (tmp > DBL_EPSILON) {
      double cos_zeta = -A_mo1_mo2/tmp;
      double sin_zeta =  B_mo1_mo2/tmp;
      double zeta = acos(cos_zeta)*((sin_zeta < 0) ? -1:1);
      *gamma = 0.25*zeta;
      *delta = A_mo1_mo2 + tmp;
    }
  else {
      *delta = 0.0;
      *gamma = 0.0;
    }

}

static void
do_rotation(double gamma,
            int nbasis, int nocc_act,
            double *vec, int mo1, int mo2)
{
  double * RESTRICT tmp_ptr1 = &vec[mo1*nbasis];
  double * RESTRICT tmp_ptr2 = &vec[mo2*nbasis];
  double cgamma = cos(gamma);
  double sgamma = sin(gamma);
  for (int aoindex=0; aoindex<nbasis; aoindex++) {
      double tmp11 = *tmp_ptr1;
      double tmp22 = *tmp_ptr2;
      *tmp_ptr1 = cgamma*tmp11 + sgamma*tmp22;
      *tmp_ptr2 = -sgamma*tmp11 + cgamma*tmp22;
      tmp_ptr1 ++;
      tmp_ptr2 ++;
    }
}

RefSCMatrix 
pop_local_mo(const Ref<OneBodyWavefunction> &wfn, int nfzc, const sc::RefSymmSCMatrix &ao_overlap,
             const sc::Ref<sc::MessageGrp> &msg)
{
  Timer tim;

  // Returns the Pipek-Mezey localized ("population-localized")
  // occupied molecular orbitals
  msg->sync(); // imbn test
  
  tim.enter("pop_local_mo");
  tim.enter("Localization setup");

  Ref<Molecule> molecule = wfn->molecule();
  Ref<GaussianBasisSet> basis = wfn->basis();

  int nbasis = basis->nbasis();
  int natoms = molecule->natom();
  int nocc = wfn->nelectron()/2; // Number of doubly occ. orbitals (only works for closed shell systems)
  int nocc_act = nocc - nfzc;    
  
  // Get the SCF vector (AO -> MO transformation), including only active occupied MO's
//imbn  RefSCMatrix scf_vector = ao_to_mo(wfn);
  tim.enter("scf_vector"); // imbn debug print
  RefSCMatrix scf_vector = ao_to_occact_mo(wfn, nfzc);
  tim.exit("scf_vector"); // imbn debug print

  tim.enter("scf_vector_dat"); // imbn debug print
  double *scf_vector_dat = new double[nbasis*nocc_act];
  RefSCMatrix scf_vector_t = scf_vector.t();
  scf_vector_t.convert(scf_vector_dat);
  tim.exit("scf_vector_dat"); // imbn debug print
  
  double threshold = 1.0e-9; // threshold for performing rotation

  //**** Only works for closed shell systems for now ****//

#if 0 // imbn
  // Compute occ orb offset array (since symmetry causes occupied orbitals
  // to be spread out in vector).
  vector<int> mo_off(nocc);
  for (int iocc=0, ij=0, i=0; i<wfn->so_dimension()->blocks()->nblock(); i++) {
    for (int j=0; j<wfn->so_dimension()->blocks()->size(i); j++,ij++) {
      if (wfn->occupation(i,j) > 1.99999) {
        mo_off[iocc++] = ij;
        }
      }
    }
#endif // imbn

  ExEnv::out0() << "ao_overlap dim:" << std::endl;
  ao_overlap.dim().print();

  ExEnv::out0() << "scf_vector dims:" << std::endl;
  scf_vector.rowdim().print();
  scf_vector.coldim().print();

  tim.enter("shalf product");  // imbn debug print

  RefSCMatrix S_half_trans = ao_overlap*scf_vector;
  tim.exit("shalf product");  // imbn debug print
  tim.enter("shalf convert");  // imbn debug print
  double *S_half_trans_dat = new double[nbasis*nocc_act];
  S_half_trans.t().convert(S_half_trans_dat);
  tim.exit("shalf convert");  // imbn debug print
  msg->sync(); // imbn test
  
  tim.exit("Localization setup");

  tim.enter("Iterations");

  double delta_max = 2*threshold; // (need init. delta_max greater than threshold)
  double current_delta, current_gamma;
              
  int rot_counter = 0;
  int upd_counter = 0;
  while (delta_max > threshold) {
      delta_max = 0.0;
      for (int mo1=0; mo1<nocc_act; mo1++) {
          for (int mo2=0; mo2<mo1; mo2++) {
              
              tim.enter("compute_delta"); // imbn debug print
              compute_delta(natoms, basis, nocc_act, scf_vector_dat,
                            S_half_trans_dat, mo1, mo2, &current_delta, &current_gamma);
              tim.exit("compute_delta"); // imbn debug print
              if (fabs(current_delta) > delta_max) delta_max = fabs(current_delta);
              
              if (fabs(current_delta) > threshold) {
                  tim.enter("do_rotation 1"); // imbn debug print
                  do_rotation(current_gamma, nbasis, nocc_act, S_half_trans_dat, mo1, mo2);
                  tim.exit("do_rotation 1"); // imbn debug print
                  tim.enter("do_rotation 2"); // imbn debug print
                  do_rotation(current_gamma, nbasis, nocc_act, scf_vector_dat, mo1, mo2);
                  rot_counter++;
                  tim.exit("do_rotation 2"); // imbn debug print
                }
              upd_counter++;
            }
        }
    }
  tim.exit("Iterations");
  ExEnv::out0() << "finished iterations" << std::endl; // imbn debug print

  //std::cout << "rot_counter(old) = " << rot_counter << std::endl;
  //std::cout << "upd_counter(old) = " << upd_counter << std::endl;
  msg->sync(); // imbn test

  scf_vector_t->assign(scf_vector_dat);
  scf_vector->assign(scf_vector_t.t());
  delete[] scf_vector_dat;

  delete[] S_half_trans_dat;
  
  // scf_vector now contains the vector of localized occupied MO's
  // (and the virtual MO's which are left unchanged in the scf_vector)

# if WRITE_ORBITALS // NOT TESTED RECENTLY AFTER CHANGES TO CODE
  // Write file with values of orbitals in xz-plane
  // Write both localized and canonical orbitals for comparison
  // First get the bounding box of the molecule *** make this a separate function
  // p1 and p2 are vectors to opposite far corners of the box

  scf_vector->print("scf_vector");
  
  SCVector3 p1, p2;
  for (int i=0; i<3; i++) p1[i] = p2[i] = molecule()->r(0,i);
  for (int i=1; i<molecule()->natom(); i++) {
    for (int j=0; j<3; j++) {
      if (molecule()->r(i,j) < p1[i]) p1[i] = molecule()->r(i,j);
      if (molecule()->r(i,j) > p2[i]) p2[i] = molecule()->r(i,j);
      }
    }
  for (int i=0; i<3; i++) {
    p1[i] = p1[i] - 3.0;
    p2[i] = p2[i] + 3.0;
    }

  RefSCMatrix scf_canonical = ao_to_mo(wfn);
  int nx = 50;
  int nz = 50;
  double delta_x = (p2[0] - p1[0])/(nx-1);
  double delta_z = (p2[2] - p1[2])/(nz-1); 
  double y=0;
  FILE *outputfp1;
  FILE *outputfp2;
  outputfp1 = fopen("scf_canonical","w");
  outputfp2 = fopen("scf_localized","w");
  for (int iocc=0; iocc<nocc; iocc++) {
    fprintf(outputfp1, "Occupied MO number %d\n", iocc);
    fprintf(outputfp2, "Occupied MO number %d\n", iocc);
    for (int j=0; j<nx; j++) {
      double x = p1[0] + j*delta_x;
      for (int k=0; k<nz; k++) {
        double z = p1[2] + k*delta_z;
        SCVector3 r;
        r(0) = x;
        r(1) = y;
        r(2) = z;
        double orbval = Wavefunction::orbital(r,mo_off[iocc],scf_canonical);
        fprintf(outputfp1, "%20.10lf %20.10lf %20.10lf\n", x,z,orbval);
        orbval = Wavefunction::orbital(r,mo_off[iocc],scf_vector);
        fprintf(outputfp2, "%20.10lf %20.10lf %20.10lf\n", x,z,orbval);
        }
      }
    }
  fclose(outputfp1);
  fclose(outputfp2);
# endif

  tim.exit("pop_local_mo");

  return scf_vector;
  }

PipekMezeyLocalization::PipekMezeyLocalization(
    const Ref<OneBodyWavefunction> &wfn, int nfzc,
    const RefSymmSCMatrix &ao_overlap)
{
  wfn_ = wfn;
  nfzc_ = nfzc;
  ao_overlap_ = ao_overlap;
}
  

void
PipekMezeyLocalization::compute_rotation(OrbData &orbdata)
{
  int mo1 = orbdata.i;
  int mo2 = orbdata.j;

  // First compute the matrix elements <mo1|P_atom|mo2>
  // of the atomic population operator P_atom and the
  // intermediates A_mo1_mo2, B_mo1_mo2
  double A_mo1_mo2 = 0.0;
  double B_mo1_mo2 = 0.0;
  int offset = 0;
  int nbasis = wfn_->basis()->nbasis();
  for (int atom=0; atom<natom_; atom++) {
      double P_mo1_mo2 = 0.0;
      double P_mo1_mo1 = 0.0;
      double P_mo2_mo2 = 0.0;
      double * RESTRICT tmp_ptr1 = &scf_vector_dat_[offset + mo1*nbasis];
      double * RESTRICT tmp_ptr2 = &scf_vector_dat_[offset + mo2*nbasis];
      double * RESTRICT tmp_ptr3 = &S_half_trans_dat_[offset + mo1*nbasis];
      double * RESTRICT tmp_ptr4 = &S_half_trans_dat_[offset + mo2*nbasis];
      int index = basis_->nbasis_on_center(atom);
      for (int ao1=0; ao1<index; ao1++) {
          double tmp1 = *tmp_ptr1;
          double tmp2 = *tmp_ptr2;
          double tmp3 = *tmp_ptr3;
          double tmp4 = *tmp_ptr4;
          P_mo1_mo2 += tmp1*tmp4 + tmp2*tmp3;
          P_mo1_mo1 += tmp1*tmp3;
          P_mo2_mo2 += tmp2*tmp4;
          tmp_ptr1 ++;
          tmp_ptr2 ++;
          tmp_ptr3 ++;
          tmp_ptr4 ++;
        }
      offset += index;
      P_mo1_mo2 *= 0.5;
      // Compute intermediates A_mo1_mo2, B_mo1_mo2
      A_mo1_mo2 += P_mo1_mo2*P_mo1_mo2 - 0.25*(P_mo1_mo1 - P_mo2_mo2)*(P_mo1_mo1 - P_mo2_mo2);
      B_mo1_mo2 += P_mo1_mo2*(P_mo1_mo1 - P_mo2_mo2);
    }
        
  // Compute rotation angle gamma and the corresponding change, delta,
  // in the functional
  double tmp = sqrt(A_mo1_mo2*A_mo1_mo2 + B_mo1_mo2*B_mo1_mo2);
  if (tmp > DBL_EPSILON) {
      double cos_zeta = -A_mo1_mo2/tmp;
      double sin_zeta =  B_mo1_mo2/tmp;
      double zeta = acos(cos_zeta)*((sin_zeta < 0) ? -1:1);
      orbdata.gamma = 0.25*zeta;
      orbdata.delta = A_mo1_mo2 + tmp;
    }
  else {
      orbdata.delta = 0.0;
      orbdata.gamma = 0.0;
    }
}

void
PipekMezeyLocalization::rotate(const OrbData &orbdata,
                               double *vec)
{
  int mo1 = orbdata.i;
  int mo2 = orbdata.j;
  double gamma = orbdata.gamma;
  int nbasis = wfn_->basis()->nbasis();
  double * RESTRICT tmp_ptr1 = &vec[mo1*nbasis];
  double * RESTRICT tmp_ptr2 = &vec[mo2*nbasis];
  double cgamma = cos(gamma);
  double sgamma = sin(gamma);
  for (int aoindex=0; aoindex<nbasis_; aoindex++) {
      double tmp11 = *tmp_ptr1;
      double tmp22 = *tmp_ptr2;
      *tmp_ptr1 = cgamma*tmp11 + sgamma*tmp22;
      *tmp_ptr2 = -sgamma*tmp11 + cgamma*tmp22;
      tmp_ptr1 ++;
      tmp_ptr2 ++;
    }
}

void
PipekMezeyLocalization::rotate(const OrbData &orbdata)
{
  rotate(orbdata,&S_half_trans_dat_.front());
  rotate(orbdata,&scf_vector_dat_.front());
}

void
PipekMezeyLocalization::init_orbmap()
{
  orbmap_iters_.resize(nocc_);
  for (int i=0; i<nocc_act_; i++) {
      orbmap_iters_[i].resize(i);
      for (int j=0; j<i; j++) {
          OrbData data;
          data.i = i;
          data.j = j;
          compute_rotation(data);
          orbmap_t::iterator newiter
              = orbmap_.insert(std::make_pair(fabs(data.delta),data));
          orbmap_iters_[i][j] = newiter;
        }
    }
}

PipekMezeyLocalization::orbmap_t::iterator
PipekMezeyLocalization::largest_orbmap_entry()
{
  return orbmap_.begin();
}

void
PipekMezeyLocalization::zero_orbmap_entry(int I, int J)
{
  OrbData data;
  data.i = I;
  data.j = J;
  data.gamma = 0.0;
  data.delta = 0.0;
  data.obsolete = false;
  orbmap_t::iterator olditer = orbmap_iters_[I][J];
  orbmap_.erase(olditer);
  orbmap_t::iterator newiter = orbmap_.insert(std::make_pair(data.delta,data));
  orbmap_iters_[I][J] = newiter;
  remove_from_obsolete(I,J);
}

void
PipekMezeyLocalization::update_orbmap_entry(int I, int J)
{
  OrbData data;
  data.i = I;
  data.j = J;
  data.obsolete = false;
  compute_rotation(data);
  orbmap_t::iterator olditer = orbmap_iters_[I][J];
  orbmap_.erase(olditer);
  orbmap_t::iterator newiter = orbmap_.insert(std::make_pair(data.delta,data));
  orbmap_iters_[I][J] = newiter;
}

void
PipekMezeyLocalization::update_orbmap_entries(int orb)
{
  for (int i=0; i<nocc_act_; i++) {
      int I,J;
      if      (i == orb) continue;
      else if (i < orb)  { I = orb; J = i; }
      else               { J = orb; I = i; }
      OrbData data;
      data.i = I;
      data.j = J;
      compute_rotation(data);
      orbmap_t::iterator olditer = orbmap_iters_[I][J];
      orbmap_.erase(olditer);
      orbmap_t::iterator newiter = orbmap_.insert(std::make_pair(data.delta,data));
      orbmap_iters_[I][J] = newiter;
      remove_from_obsolete(I,J);
    }
}

void
PipekMezeyLocalization::obsolete_orbmap_entries(int orb)
{
  for (int i=0; i<nocc_act_; i++) {
      int I,J;
      if      (i == orb) continue;
      else if (i < orb)  { I = orb; J = i; }
      else               { J = orb; I = i; }
      orbmap_iters_[I][J]->second.obsolete = true;
      add_to_obsolete(I,J);
    }
}

long
PipekMezeyLocalization::update_obsolete_entries()
{
  int nupdate = 0;
  for (std::set<std::pair<int,int> >::iterator iter = obsolete_orbs_.begin();
       iter != obsolete_orbs_.end();
       iter++) {
      update_orbmap_entry(iter->first, iter->second);
      nupdate++;
    }
  obsolete_orbs_.clear();
  return nupdate;
}

void
PipekMezeyLocalization::add_to_obsolete(int i, int j)
{
  obsolete_orbs_.insert(std::make_pair(i,j));
}

void
PipekMezeyLocalization::remove_from_obsolete(int i, int j)
{
  obsolete_orbs_.erase(std::make_pair(i,j));
}

RefSCMatrix 
PipekMezeyLocalization::compute_orbitals()
{
  Timer overall_timer("PipekMezeyLocalization::orbitals");

  // Returns the Pipek-Mezey localized ("population-localized")
  // occupied molecular orbitals

  Timer section_timer("Localization setup");

  Ref<Molecule> molecule = wfn_->molecule();

  basis_ = wfn_->basis();
  nbasis_ = basis_->nbasis();
//imbn  noso_ = wfn_->oso_dimension();
  natom_ = molecule->natom();
  nocc_ = wfn_->nelectron()/2; // number of doubly occupied orbitals
  nocc_act_ = nocc_ - nfzc_;
  
  // Get the SCF vector (AO -> MO transformation), including active occupied MO's only
//imbn  scf_vector_ = ao_to_mo(wfn_);
  scf_vector_ = ao_to_occact_mo(wfn_, nfzc_);

  scf_vector_dat_.resize(nbasis_*nocc_act_);
  RefSCMatrix scf_vector_t = scf_vector_.t();
  scf_vector_t.convert(&scf_vector_dat_.front());
  
  double threshold = 1.0e-9; // threshold for performing rotation
  
  //**** Only works for closed shell systems for now ****//

  RefSCMatrix S_half_trans = ao_overlap_*scf_vector_;
  S_half_trans_dat_.resize(nbasis_*nocc_act_);
  S_half_trans.t().convert(&S_half_trans_dat_.front());

  section_timer.change("Iterations");

  init_orbmap();

  double last_complete_orbmap_delta = largest_orbmap_entry()->first;
  int rot_counter = 0;
  int upd_counter = 0;
  for (orbmap_t::iterator largest_delta = largest_orbmap_entry();
       (largest_delta->first > threshold
        || (update_obsolete_entries()
            && (largest_delta=largest_orbmap_entry())->first > threshold));
       largest_delta = largest_orbmap_entry()) {
      if (last_complete_orbmap_delta > 1.0e3 * largest_delta->first) {
          upd_counter += update_obsolete_entries();
          last_complete_orbmap_delta = largest_orbmap_entry()->first;
        }
      else {
          if (largest_delta->second.obsolete) {
              update_orbmap_entry(largest_delta->second.i,
                                  largest_delta->second.j);
              continue;
            }
          rotate(largest_delta->second);
          obsolete_orbmap_entries(largest_delta->second.i);
          obsolete_orbmap_entries(largest_delta->second.j);
          zero_orbmap_entry(largest_delta->second.i,
                            largest_delta->second.j);
          rot_counter++;
        }
    }

  std::cout << "rot_counter(new) = " << rot_counter << std::endl;
  std::cout << "upd_counter(new) = " << upd_counter << std::endl;

  section_timer.exit();

  scf_vector_t->assign(&scf_vector_dat_.front());
  scf_vector_.assign(scf_vector_t.t());

  scf_vector_dat_.clear();
  S_half_trans_dat_.clear();
  orbmap_.clear();
  orbmap_iters_.clear();

  // scf_vector now contains the vector of localized occupied MO's
  // (and the virtual MO's which are left unchanged in the scf_vector)

  return scf_vector_;
}

// NOT TESTED RECENTLY AFTER CHANGES TO CODE
// The compute_orbitals member must be called first.
// Write file with values of orbitals in xz-plane
// Write both localized and canonical orbitals for comparison
// First get the bounding box of the molecule *** make this a separate function
// p1 and p2 are vectors to opposite far corners of the box
void
PipekMezeyLocalization::write_orbitals()
{
  scf_vector_->print("scf_vector");
  
  SCVector3 p1, p2;
  for (int i=0; i<3; i++) p1[i] = p2[i] = wfn_->molecule()->r(0,i);
  for (int i=1; i<wfn_->molecule()->natom(); i++) {
    for (int j=0; j<3; j++) {
      if (wfn_->molecule()->r(i,j) < p1[i]) {
          p1[i] = wfn_->molecule()->r(i,j);
        }
      if (wfn_->molecule()->r(i,j) > p2[i]) {
          p2[i] = wfn_->molecule()->r(i,j);
        }
      }
    }
  for (int i=0; i<3; i++) {
    p1[i] = p1[i] - 3.0;
    p2[i] = p2[i] + 3.0;
    }

  RefSCMatrix scf_canonical = ao_to_mo(wfn_);

  int nx = 50;
  int nz = 50;
  double delta_x = (p2[0] - p1[0])/(nx-1);
  double delta_z = (p2[2] - p1[2])/(nz-1); 
  double y=0;
  FILE *outputfp1;
  FILE *outputfp2;
  outputfp1 = fopen("scf_canonical","w");
  outputfp2 = fopen("scf_localized","w");
  for (int iocc=0; iocc<nocc_; iocc++) {
    fprintf(outputfp1, "Occupied MO number %d\n", iocc);
    fprintf(outputfp2, "Occupied MO number %d\n", iocc);
    for (int j=0; j<nx; j++) {
      double x = p1[0] + j*delta_x;
      for (int k=0; k<nz; k++) {
        double z = p1[2] + k*delta_z;
        SCVector3 r;
        r(0) = x;
        r(1) = y;
        r(2) = z;
        }
      }
    }
  fclose(outputfp1);
  fclose(outputfp2);
}

}
