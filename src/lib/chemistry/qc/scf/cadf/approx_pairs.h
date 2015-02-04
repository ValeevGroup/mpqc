//
// approx_pairs.h
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Jun 6, 2014
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

#ifndef _chemistry_qc_scf_cadf_approx_pairs_h
#define _chemistry_qc_scf_cadf_approx_pairs_h

#include "cadfclhf.h"

namespace sc {

class ApproximatePairWriter : public DescribedClass {
  private:

    Ref<CADFCLHF> wfn_;

    std::string filename_;

    Ref<GaussianBasisSet> exbs_;

    // Pairs of shell indices to print
    std::vector<std::pair<int, int>> pairs_;

    Ref<TwoBodyTwoCenterInt> ints2c_ex_;
    Ref<TwoBodyThreeCenterInt> ints3c_ex_;

    std::unordered_map<int, int> shell_map_gbs_;
    std::unordered_map<int, int> shell_map_dfbs_;
    std::unordered_map<int, int> shell_map_exbs_;
    std::vector<std::array<int,4>> pairs_tmp_;

    void initialize();

    bool initialized_ = false;

    void compute_pairs_ex();

    void write_atoms_section(std::ostream &out);
    void write_gto_section(std::ostream &out);
    void write_mo_section(std::ostream &out);

    bool normalize_pairs_ = true;

    std::unordered_map<std::pair<int,int>, Eigen::MatrixXd, sc::hash<std::pair<int, int>>> excoefs_;

    enum {
      GBS_orbs,
      DFBS_pair,
      EX_pair
    };

    void write_mo(std::ostream& out, const Eigen::VectorXd& coefs, int orb_type, const std::string& name);

  public:

    ApproximatePairWriter(const Ref<KeyVal>& keyval);

    void write_pairs();


    friend class CADFCLHF;
    static ClassDesc cd_;

};

} // end namespace sc

#endif /* _chemistry_qc_scf_cadf_approx_pairs_h */
