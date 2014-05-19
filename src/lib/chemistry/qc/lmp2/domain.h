
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

#ifndef _chemistry_qc_lmp2_domain_h
#define _chemistry_qc_lmp2_domain_h

#include <chemistry/qc/wfn/obwfn.h>
#include <vector>
#include <map>
#include <set>
#include <chemistry/qc/lmp2/sma.h>

namespace sc {

typedef std::set<int> domainmapvirbs_t;
typedef std::map<std::pair<int,int>, domainmapvirbs_t > domainmap_t;


/// \brief Create maps of occupied orbital pairs to atoms in their domain.
void
create_domains(const sc::Ref<sc::OneBodyWavefunction> &wfn, int nfzc,
               const sc::RefSCMatrix &scf_local,
               std::vector<std::vector<int> > &domains, domainmap_t &domainmap,
               double distance_threshold, double completeness_threshold,
               bool all_nondist_pairs,
               std::vector<double>& interdomain_distances,
               sma2::Array<2> &S_ao, sma2::Array<2> &L,
               double bound, const sc::Ref<sc::MessageGrp> &msg);


/// \brief Create maps of occupied orbital pairs to atoms in their
/// domain for both strong and weak pairs.
void
create_domains(const sc::Ref<sc::OneBodyWavefunction> &wfn, int nfzc,
               const sc::RefSCMatrix &scf_local,
               std::vector<std::vector<int> > &domains,
               domainmap_t &domainmap, double distance_threshold,
               double completeness_threshold,
               domainmap_t &weak_pair_domainmap, double weak_pair_distance_threshold,
               bool all_nondist_pairs, std::vector<double>& interdomain_distances,
               sma2::Array<2> &S_ao, sma2::Array<2> &L,
               double bound, const sc::Ref<sc::MessageGrp> &msg);

}

#endif
