//
// tiledbasisset.cpp
//
// Copyright (C) 2013 Drew Lewis
//
// Authors: Drew Lewis
// Maintainer: Drew Lewis and Edward Valeev
//
// This file is part of the MPQC Toolkit.
//
// The MPQC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The MPQC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the MPQC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <chemistry/qc/basis/tiledbasisset.hpp>
#include <Eigen/Dense>
#include <functional>
#include <numeric>
#include <vector>
#include <string>

using namespace mpqc;
using namespace mpqc::TA;

static sc::ClassDesc TiledBasisSet_cd( typeid(TiledBasisSet), "TiledBasisSet",
                1, "public GaussianBasisSet",
                0,  sc::create<TiledBasisSet>, sc::create<TiledBasisSet>);

TiledBasisSet::TiledBasisSet(const sc::Ref<sc::KeyVal> &keyval):
    ntiles_(keyval->intvalue("ntiles", sc::KeyValValueint(1))),
    SRange_()
{
    Basis basis;
    basis << keyval->describedclassvalue("basis");
    if(basis.null()){
        basis = new sc::GaussianBasisSet(keyval);
        if(basis.null()){
            throw sc::InputError("Could not construct a GaussianBasisSet",
                                 __FILE__, __LINE__,
                                 "basis", 0, class_desc());
        }
    }
    // make sure that molecule_ is set so that can construct Shell objects
    molecule_ = basis->molecule();

    ShellOrder ordering(basis);
    std::vector<Shell> shells = ordering.ordered_shells(ntiles_, this);
    SRange_ = ordering.shell_ranges();

    init(converted_name(basis->name()),
         converted_name(basis->label()),
         basis->molecule(),
         shells);
}

TiledBasisSet::TiledBasisSet(const sc::Ref<sc::GaussianBasisSet>& bs,
                             size_t ntiles) :
                ntiles_(ntiles),
                SRange_()
{
  // make sure that molecule_ is set so that can construct Shell objects
  molecule_ = bs->molecule();
  ShellOrder ordering(bs);
  std::vector<Shell> shells = ordering.ordered_shells(ntiles_, this);
  SRange_ = ordering.shell_ranges();

  init(converted_name(bs->name()),
       converted_name(bs->label()),
       bs->molecule(),
       shells);
}

std::string TiledBasisSet::converted_name(const std::string& name) {
    if(name.empty()) return name;
    std::string new_name = "Tiled(";
    new_name += name;
    new_name += ")";
    return new_name;
}

TiledArray::TiledRange1 TiledBasisSet::trange1() const{
    std::vector<std::size_t> tilesizes;
    tilesizes.reserve(ntiles_);

    // Loop over clusters
    for(auto i = 0; i < ntiles_; ++i){
       // Get the first function in the shell
       tilesizes.push_back(shell_to_function(SRange_[i]));
    }

    // Get the last function in the shell since our loop doesn't cover it.
    tilesizes.push_back(nbasis());

    return TiledArray::TiledRange1(tilesizes.begin(), tilesizes.end());

}

void
TiledBasisSet::print(std::ostream &os_in ) const {
  // Print base information
  sc::GaussianBasisSet::print(os_in);

  // Create vector of clusters use ref_wrapper so as to not copy atoms
  std::vector< unsigned int > clusters;
  clusters.reserve(ntiles_);

  os_in << sc::incindent;

  // Loop over Tiles to get clusters
  for(auto tile = 0; tile < ntiles_; ++tile){
    // Temp to hold the number of functions on a given cluster
    unsigned int funcs_on_tile = 0;
    // Loop over the shells on the tile
    for(auto shell_i = SRange_[tile]; shell_i < SRange_[tile+1]; ++shell_i){
      funcs_on_tile += shell(shell_i).nfunction();
    }
    clusters.push_back(funcs_on_tile);
  }

  auto begin = std::begin(clusters);
  auto end = std::end(clusters);
  auto size = clusters.size();

  double cluster_avg = static_cast<double>(std::accumulate(begin, end, 0ul))/
                                                   static_cast<double>(size);
  double sq_sum = std::inner_product(begin, end, begin, 0);
  double std_dev = std::sqrt(sq_sum / size - cluster_avg * cluster_avg);
  auto cluster_minmax = std::minmax_element(begin, end);

  os_in << "Average Cluster Size = " << cluster_avg
     << "\n\tStandard Deviation =   " << std_dev
     << "\n\tMinimum cluster size = " << *cluster_minmax.first
     << "\n\tMaximum cluster size = " << *cluster_minmax.second
     << std::endl;
}

void TiledBasisSet::save_data_state(sc::StateOut& s){
    sc::GaussianBasisSet::save_data_state(s);
}

TiledBasisSet::TiledBasisSet(sc::StateIn& s):
    sc::SavableState(s), sc::GaussianBasisSet(s)
{}



