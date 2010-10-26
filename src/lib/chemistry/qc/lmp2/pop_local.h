
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

#ifndef _chemistry_qc_lmp2_pop_local_h
#define _chemistry_qc_lmp2_pop_local_h

#include <chemistry/qc/wfn/obwfn.h>

namespace sc {

sc::RefSCMatrix
convert_complete_to_occupied_vector(
    const sc::Ref<sc::OneBodyWavefunction> &wfn,
    const sc::RefSCMatrix &vec);

sc::RefSCMatrix
convert_complete_to_occupied_vector_nosymm(
    const sc::Ref<sc::OneBodyWavefunction> &wfn,
    int nfzc, const sc::RefSCMatrix &vec);


sc::RefSCMatrix pop_local_mo(const sc::Ref<sc::OneBodyWavefunction> &wfn, int nfzc, const sc::RefSymmSCMatrix &ao_overlap,
                             const sc::Ref<sc::MessageGrp> &msg);

/// \brief Performs a Pipek-Mezey orbital localization.
class PipekMezeyLocalization {
    sc::Ref<sc::OneBodyWavefunction> wfn_;
    int nfzc_;
    sc::RefSymmSCMatrix ao_overlap_;

    // filled in by compute_orbitals
    sc::RefSCMatrix scf_vector_;

    class OrbData {
      public:
        double delta,gamma;
        int i,j;
        mutable bool obsolete;
        OrbData(): obsolete(true) {};
    };

    // temporary data and routines used by compute_orbitals
    typedef std::multimap<double, OrbData, std::greater<double> > orbmap_t;
    orbmap_t orbmap_;
    std::set<std::pair<int,int> > obsolete_orbs_;
    std::vector<std::vector<orbmap_t::iterator> > orbmap_iters_;
    std::vector<double> scf_vector_dat_, S_half_trans_dat_;
    int nocc_;
    int nocc_act_;
    int natom_;
    int noso_;
    int nbasis_;
    sc::Ref<sc::GaussianBasisSet> basis_;
    void compute_rotation(OrbData &orbdata);
    void rotate(const OrbData &orbdata, double *vec);
    void rotate(const OrbData &orbdata);
    void init_orbmap();
    orbmap_t::iterator largest_orbmap_entry();
    void update_orbmap_entries(int orb);
    void obsolete_orbmap_entries(int orb);
    void update_orbmap_entry(int I, int J);
    void zero_orbmap_entry(int I, int J);
    void remove_from_obsolete(int i, int j);
    void add_to_obsolete(int i, int j);
    long update_obsolete_entries();

  public:
    PipekMezeyLocalization(const sc::Ref<sc::OneBodyWavefunction> &wfn,
                           int nfzc,
                           const sc::RefSymmSCMatrix &ao_overlap);
    sc::RefSCMatrix compute_orbitals();
    sc::RefSCMatrix write_orbitals();
};

}

#endif
