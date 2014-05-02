
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

#include <vector>
#include <algorithm>
#include <utility>
#include <util/misc/autovec.h>
#include <util/state/state_bin.h>
#include <util/group/mstate.h>
#include <math/scmat/matrix.h>
#include <math/scmat/vector3.h>
#include <chemistry/molecule/formula.h>
#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/molecule/energy.h>

#include <chemistry/qc/lmp2/pop_local.h>
#include <chemistry/qc/lmp2/domain.h>
#include <chemistry/qc/lmp2/sma.h>
#include <chemistry/qc/lmp2/parallel.h>

using namespace std;

namespace sc {

class GreaterAbsVal {
  public:
    bool operator()(double val1, double val2) const {
      return fabs(val1) > fabs(val2);
    }
};

static double
compute_residual(const RefSCMatrix &overlap_times_scf,
                 const RefSymmSCMatrix &S,
                 const vector<int> &domain, int *index_array,
                 const Ref<GaussianBasisSet> &basis, int mo_index)
{
  // Compute the residual for the domain of molecular orbital mo
  // This is computed from Eq. (9) in Hampel&Werner

  int domain_length = domain.size();
  int nbasis = basis->nbasis();

  Ref<SCMatrixKit> kit = overlap_times_scf->kit();

  // compute the number of basis functions in the domain
  int nbasis_in_domain = 0;
  for (int j=0; j<domain_length; j++) {
      int tmp = basis->nbasis_on_center(index_array[j]);
      nbasis_in_domain += tmp;
    }
  // Create an array of offsets for the basis functions on the
  // centers in index_array
  std::vector<int> center_offsets(domain_length);
  for (int j=0; j<domain_length; j++) {
    center_offsets[j] = 0;
    for (int k=0; k<index_array[j]; k++) {
      center_offsets[j] += basis->nbasis_on_center(k);
      }
    }

  // overlap_times_scf_domain is a part of overlap_times_scf:
  // row index in domain, col. index=mo.
  // NB: atoms in domain are ordered as listed in index_array, according
  // to their contribution to the current MO
  RefSCDimension nbasis_in_domain_dim(new SCDimension(nbasis_in_domain,1));
  nbasis_in_domain_dim->blocks()->set_subdim(0,new SCDimension(nbasis_in_domain,1));
  RefSCVector overlap_times_scf_domain(nbasis_in_domain_dim,kit);
  int index = 0;
  for (int j=0; j<domain_length; j++) {
    int k_offset = center_offsets[j];
    for (int k=0; k<basis->nbasis_on_center(index_array[j]); k++) {
      overlap_times_scf_domain(index) = overlap_times_scf(k+k_offset, mo_index);
      index++;
      }
    }
  // S_domain is a part of S; row and col. indices in domain
  RefSCMatrix S_domain(nbasis_in_domain_dim,nbasis_in_domain_dim,kit);
  int index1=0;
  for (int j=0; j<domain_length; j++) {
    int k_offset = center_offsets[j];
    for (int k=0; k<basis->nbasis_on_center(index_array[j]); k++) {
      int index2 = 0;
      for (int l=0; l<domain_length; l++) {
        int m_offset = center_offsets[l];
        for (int m=0; m<basis->nbasis_on_center(index_array[l]); m++) {
            S_domain(index1,index2) = S(k+k_offset,m+m_offset);
            index2++;
          }
        }
      index1++;
      }
    }

  RefSCVector R_prime = overlap_times_scf_domain.copy();
  S_domain.solve_lin(R_prime);

  double residual = 1.0 - R_prime.dot(overlap_times_scf_domain);

  return residual;
}

static double
domain_distance(const vector<int> &domain1, const vector<int> &domain2,
                const Ref<Molecule> &molecule)
{
  // Returns the distance between the domains domain1 and domain2.
  // The distance is defined as the minimum distance between an atom
  // in domain1 and an atom in domain2

  SCVector3 ri = molecule->r(domain1[0]);
  SCVector3 rj = molecule->r(domain2[0]);
  SCVector3 rij = ri - rj;
  double minimum_distance = rij.norm();

  for (int i=0; i<domain1.size(); i++) { // loop over atoms in domain1
    for (int j=0; j<domain2.size(); j++) { // loop over atoms in domain2
      ri = molecule->r(domain1[i]);
      rj = molecule->r(domain2[j]);
      rij = ri - rj;
      double distance = rij.norm();
      if (distance < minimum_distance) minimum_distance = distance;
      }
    }
  return minimum_distance;
}

static bool
very_distant(const vector<int> &domain1, const vector<int> &domain2,
             const Ref<Molecule> &molecule, double distance_threshold)
{
  // Return 1 if domain1 and domain2 are very distant, 0 otherwise.
  // The domains are very distant if the distance between all atom
  // pairs ij (i in domain1, j in domain2) is greater than or equal to
  // distance_threshold

  for (int i=0; i<domain1.size(); i++) { // loop over atoms in domain1
    for (int j=0; j<domain2.size(); j++) { // loop over atoms in domain2
      SCVector3 ri = molecule->r(domain1[i]);
      SCVector3 rj = molecule->r(domain2[j]);
      SCVector3 rij = ri - rj;
      double distance = rij.norm();
      if (distance < distance_threshold) return 0;
      }
    }
  return 1;
}

void
create_domains(const Ref<OneBodyWavefunction> &wfn, int nfzc,
               const RefSCMatrix &scf_local, vector<vector<int> > &domains,
               domainmap_t &domainmap, double distance_threshold,
               double completeness_threshold,
               bool all_nondist_pairs, std::vector<double>& interdomain_distances,
               sma2::Array<2> &S_ao,
               sma2::Array<2> &L, double bound, const sc::Ref<sc::MessageGrp> &msg)
{
  domainmap_t weak_pair_domainmap;
  create_domains(wfn, nfzc, scf_local, domains,
                 domainmap, distance_threshold,
                 completeness_threshold,
                 weak_pair_domainmap, -1.0,
                 all_nondist_pairs, interdomain_distances,
                 S_ao, L, bound, msg);
}

void
create_domains(const Ref<OneBodyWavefunction> &wfn, int nfzc,
               const RefSCMatrix &scf_local, vector<vector<int> > &domains,
               domainmap_t &domainmap, double distance_threshold,
               double completeness_threshold,
               domainmap_t &weak_pair_domainmap, double weak_pair_distance_threshold,
               bool all_nondist_pairs, std::vector<double>& interdomain_distances,
               sma2::Array<2> &S_ao, sma2::Array<2> &L,
               double bound, const sc::Ref<sc::MessageGrp> &msg)
{
  const bool store_distances = (interdomain_distances.empty() == false);
  // Creates the vector "domains" containing the orbital domain
  // associated with each occupied localized MO (the domain of an MO
  // is here defined as the set of atoms whose AO's will constitute
  // the excitation domain for that MO);
  // also creates "domainmap" containing the coresponding pair domains.
  // NB: if "all_nondist_pairs" is true, all non-distant ij pairs are included
  //     in the domainmap; otherwise, the restriction j<=i is used.
  //     Currently, all_nondist_pairs must be true for running the lccsd code
  //     and false for running the lmp2 code.
  // The employed procedure is outlined in Hampel&Werner, JCP v.104,
  // p.6286 (1996).
  // "domains" is a vector with nocc_acc elements; each element is an
  // integer vector: the elements of the vector are the numbers of the
  // atoms in the domain for this mo

  Timer tim("various inits");

  Ref<GaussianBasisSet> basis = wfn->basis();

  int natoms = wfn->molecule()->natom();
  double occ_threshold = 0.9;           // **** may need to adjust this - or read in
  double mo_pop;

  int nocc = wfn->nelectron()/2; // number of doubly occupied orbitals
  int nocc_act = nocc - nfzc;

  //**** Only works for closed shell systems for now ****//

  int nbasis = basis->nbasis();

  Ref<PetiteList> pl = wfn->integral()->petite_list();

  tim.change("overlap matrix");
  // Compute the AO overlap matrix S
  RefSymmSCMatrix blockedS = pl->to_AO_basis(wfn->overlap());
  RefSymmSCMatrix S = blockedS;
  blockedS = 0;

  tim.change("overlap_times_scf");
  // Compute product
  sma2::Range ao_range(S_ao.index(0));
  sma2::Range occ_act_range(L.index(1));
  sma2::Array<2> overlap_times_scf_array(ao_range, occ_act_range, "overlap_times_scf_array", bound);
  overlap_times_scf_array("r","i") |= S_ao("r","p") * L("p","i");
  overlap_times_scf_array.zero();
  int nproc = msg->n();
  int me = msg->me();
  for (int index=0; index<nbasis; index++) {
      if (index%nproc != me) continue;
      sma2::Index fixed_index("index",index);
      sma2::Index i("i"), p("p");
      overlap_times_scf_array(fixed_index,i).skip_bounds_update()
          += S_ao(fixed_index,p) * L(p,i);
    }
  overlap_times_scf_array.parallel_accumulate(msg);
  // Pack into RefSCMatrix
  Ref<SCMatrixKit> kit = basis->matrixkit();
  RefSCDimension ao_dim(new SCDimension(basis->nbasis()));
  RefSCDimension occ_act_dim(new SCDimension(nocc_act));
  RefSCMatrix overlap_times_scf(ao_dim, occ_act_dim, kit);
  pack_array_into_matrix(overlap_times_scf_array, overlap_times_scf);

  tim.change("create domain for each MO");
  // Create the domain associated with each occupied MO
  for (int i=nfzc; i<nocc; i++) { // loop over active occupied MO's

    // ordered_pop contains Mulliken populations for the current MO
    // and the corresponding atoms
    // (populations ordered according to decreasing absolute value)
    std::multimap<double, int, GreaterAbsVal> ordered_pop;

    int k_offset = 0;
    for (int j=0; j<natoms; j++) { // compute contrib. to population of i
      // from AO's on atom j
      double pop = 0.0;
      for (int k=0; k<basis->nbasis_on_center(j); k++) {
        pop += scf_local(k+k_offset, i-nfzc) * overlap_times_scf(k+k_offset, i-nfzc);
      }
      ordered_pop.insert(std::make_pair(pop,j));
      k_offset += basis->nbasis_on_center(j);
    }

    // Add atoms to the domain
    mo_pop = 0.0;
    for (std::multimap<double, int, GreaterAbsVal>::iterator iter = ordered_pop.begin();
         iter != ordered_pop.end(); iter++) {
      mo_pop += iter->first;
      domains[i-nfzc].push_back(iter->second);
      if (mo_pop >= occ_threshold) break;
    }

    // Create index_array which is an array of atoms, ordered as they
    // are in ordered_pop
    int *index_array = new int[natoms];
    int index = 0;
    for (std::multimap<double, int, GreaterAbsVal>::iterator iter=ordered_pop.begin();
         iter != ordered_pop.end(); iter++, index++) {
      index_array[index] = iter->second;
    }

    // Test for completeness of the domain and expand domain, if necessary
    double residual = compute_residual(overlap_times_scf, S, domains[i-nfzc],
                                       index_array, basis, i-nfzc);
    while (residual > completeness_threshold) {
      int domain_length = domains[i-nfzc].size();
      if (domain_length == natoms) {
        std::cout << indent << "All atoms included in domain, but residual is "
          "still too large - quitting" << endl ;
        abort();
      }
      domains[i-nfzc].push_back(index_array[domain_length]); // add atom to domain
      residual = compute_residual(overlap_times_scf, S, domains[i-nfzc], index_array,
                                  basis, i-nfzc);
    }

    delete[] index_array;

#if 1
    // imbn debug print
    ExEnv::out0() << indent << scprintf("  MO   Atoms (and pop) in domain") << std::endl;
    int tmpindex=0;
    for (std::multimap<double, int, GreaterAbsVal>::iterator iter=ordered_pop.begin();
         iter != ordered_pop.end(); iter++) {
      ExEnv::out0() << scprintf("%4d (%6.3f) ", domains[i-nfzc][tmpindex], iter->first);
      tmpindex++;
      if (tmpindex == domains[i-nfzc].size()) break;
    }
    ExEnv::out0() << std::endl;
    // imbn end of debug print
#endif

  } // End i loop

  // Print out the domains
  ExEnv::out0() << indent << scprintf("Orbital Domains:") << std::endl;
  ExEnv::out0() << indent << scprintf("  MO            Atoms in domain"
                                      " (in order of decreasing contrib. "
                                      "to the MO pop.)") << std::endl;
  ExEnv::out0() << indent << scprintf("------         -----------------") << std::endl;
  for (int i=nfzc; i<nocc; i++) {
      ExEnv::out0() << indent << scprintf("%4d          ", i);
    for (int j=0; j<domains[i-nfzc].size(); j++) {
        ExEnv::out0() << scprintf("%4d ", domains[i-nfzc][j]);
      }
    ExEnv::out0() << std::endl;
    }
  ExEnv::out0() << std::endl;

  // Create pair domains
  // A pair domain for orbitals i and j is the union of domains of i and j,
  // Redundant orbitals (due to the redundancy of the projected AO's) will
  // be removed at a later stage

  tim.change("create pair domains");

  // First find out which pairs to include
  vector<pair<int,int> > included_pairs; // occ-occ pairs that are not "very distant"
  vector<pair<int,int> > weak_pairs; // "weak" occ-occ pairs
  int jmax;
  for (int i=nfzc; i<nocc; i++) {
      jmax = all_nondist_pairs ? nocc-1:i;
      for (int j=nfzc; j<=jmax; j++) {
          double dist
              = domain_distance(domains[i-nfzc], domains[j-nfzc], wfn->molecule());
          if(store_distances)
            interdomain_distances[(i-nfzc)*nocc_act + (j-nfzc)] = dist;
          if (dist < distance_threshold) {
              pair<int,int> current_pair(i,j); // the current pair
              included_pairs.push_back(current_pair);
            }
          else if (dist < weak_pair_distance_threshold) {
              pair<int,int> current_pair(i,j); // the current pair
              weak_pairs.push_back(current_pair);
            }
        }
    }
  int npairs = included_pairs.size();

  // sort the domains in order of decreasing elements
  for (int i=0; i<domains.size(); i++) {
      sort(domains[i].begin(), domains[i].end());
    }

  vector<vector<int> > pair_domains;
  pair_domains.resize(npairs);
  for (int index=0; index<npairs; index++) {
    int i = included_pairs[index].first;
    int j = included_pairs[index].second;
    // create the union of domains of i and j
    set_union(domains[i-nfzc].begin(), domains[i-nfzc].end(),
              domains[j-nfzc].begin(), domains[j-nfzc].end(),
              back_inserter(pair_domains[index]));
    domainmap[std::make_pair(i,j)].insert(pair_domains[index].begin(),
                                          pair_domains[index].end());
    }

  for (int index=0; index<weak_pairs.size(); index++) {
    int i = weak_pairs[index].first;
    int j = weak_pairs[index].second;
    // create the union of domains of i and j
    vector<int> tmp_pair_domains;
    set_union(domains[i-nfzc].begin(), domains[i-nfzc].end(),
              domains[j-nfzc].begin(), domains[j-nfzc].end(),
              back_inserter(tmp_pair_domains));
    weak_pair_domainmap[std::make_pair(i,j)].insert(tmp_pair_domains.begin(),
                                                    tmp_pair_domains.end());
    }

//    // Print out the pair domains
//    ExEnv::out0() << indent << scprintf("Pair Domains:") << endl;
//    ExEnv::out0() << indent << scprintf("  MO Pair       Atoms in domain") << std::endl;
//    ExEnv::out0() << indent << scprintf("---------      -----------------") << std::endl;
//    for (int index=0; index<npairs; index++) {
//      int i = included_pairs[index].first;
//      int j = included_pairs[index].second;
//      ExEnv::out0() << indent << scprintf("%4d %4d     ", i, j);
//      for (int k=0; k<pair_domains[index].size(); k++) {
//        ExEnv::out0() << scprintf("%4d ", pair_domains[index][k]);
//        }
//      std::cout << endl;
//      }
//    std::cout << endl;

  // Compute the number of double substitution amplitudes
  double ndoubles = 0;
  for (int i=0; i<npairs; i++) {
    int nbasis_in_pair_domain = 0;
    for (int j=0; j<pair_domains[i].size(); j++) {
      nbasis_in_pair_domain += basis->nbasis_on_center(pair_domains[i][j]);
      }
    ndoubles += nbasis_in_pair_domain*nbasis_in_pair_domain;
    }


  // Print out the number of pairs and double substitution amplitudes
  int npairs_total;
  if (all_nondist_pairs) {
      npairs_total = (nocc-nfzc)*(nocc-nfzc);
    }
  else {
      npairs_total = ((nocc-nfzc)*(nocc-nfzc+1))/2;
    }
  ExEnv::out0() << indent << scprintf("Total number of ij pairs:    %10d", npairs_total) << std::endl;
  ExEnv::out0() << indent << scprintf("Number of included ij pairs: %10d", npairs) << std::endl;
  ExEnv::out0() << indent << scprintf("Number of local double substitution "
                                      "amplitudes:    %15.0f", ndoubles) << std::endl;
  int nvir = basis->nbasis() - nocc;
  ndoubles = (((double)nocc*nvir)*(nocc*nvir+1))/2;
  ExEnv::out0() << indent << scprintf("Number of canonical double substitution "
                                      "amplitudes: %15.0f", ndoubles) << std::endl;

  }

}
