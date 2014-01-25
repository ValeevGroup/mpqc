//
// tiledbasisset.hpp
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

#ifndef MPQC_BASIS_TILEDBASISSET_HPP
#define MPQC_BASIS_TILEDBASISSET_HPP

#include <tiled_array.h>
#include <vector>

#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/molecule/atom.h>
#include <chemistry/molecule/molecule.h>

#include <Eigen/Dense>
#include <mpqc/utility/foreach.hpp>
#include "shellorder.hpp"
#include "kcluster.hpp"
#include <string>

namespace mpqc{

    /**
     * Will create a basis using k-means clustering.
     */
    class TiledBasisSet : public sc::GaussianBasisSet {

    public:
        typedef basis::ShellOrder::Shell Shell;
        typedef sc::Ref<sc::GaussianBasisSet> Basis;
        typedef basis::ShellOrder::ShellRange ShellRange;

        /**
         * Constructs a TiledBasisSet from a sc::Keyval object.
         * @param[in] keyval is a sc::Ref<sc::Keyval> object which contains
         *     information about the molecule and basis set desired.
         * @param[in] nclusters is the number of KClusters the molecule should
         *     be divided into.
         */
        TiledBasisSet(const sc::Ref<sc::KeyVal> &keyval):
            SRange_(), nclusters_(keyval->intvalue("cluster_size",
                                                   sc::KeyValValueint(2)))
        {
            Basis basis;
            basis << keyval->describedclassvalue("basis");
            if(basis == 0){
                basis = new sc::GaussianBasisSet(keyval);
                if(basis == 0){
                    throw sc::InputError("Could not construct a GaussianBasisSet",
                                     __FILE__, __LINE__,
                                     "basis", 0, class_desc());
                }
            }

            basis::ShellOrder ordering(basis);
            std::vector<Shell> shells = ordering.ordered_shells(nclusters_);
            SRange_ = ordering.shell_ranges();

            init(name_conv_TBS(basis->name()),
                 name_conv_TBS(basis->label()),
                 basis->molecule(),
                 shells);
        }

        static
        std::string
        name_conv_TBS(const std::string& name){
            if(name.empty()) return name;
            std::string new_name = "TiledArray(";
            new_name += name;
            new_name += ")";
            return new_name;
        }

        /**
         * Returns a TiledArray::TiledRange1 for tiles made by k-means clustering.
         */
        TiledArray::TiledRange1
        trange1() const {
            std::vector<std::size_t> tilesizes;
            tilesizes.reserve(nclusters_);

            // Loop over clusters
            for(auto i = 0; i < nclusters_; ++i){
               // Get the first function in the shell
               tilesizes.push_back(shell_to_function(SRange_[i]));
            }

            // Get the last function in the shell since our loop doesn't cover it.
            tilesizes.push_back(nbasis());

            return TiledArray::TiledRange1(tilesizes.begin(), tilesizes.end());

        }

        TiledBasisSet(sc::StateIn& s):
            sc::SavableState(s), sc::GaussianBasisSet(s)
        {}

        void save_data_state(sc::StateOut& s){
            sc::GaussianBasisSet::save_data_state(s);
        }

    private:
        std::size_t nclusters_;
        ShellRange SRange_;
    };

    static sc::ClassDesc TiledBasisSet_cd(
                    typeid(TiledBasisSet), "TiledBasisSet", 1,
                    "public GaussianBasisSet",
                    0,  sc::create<TiledBasisSet>, sc::create<TiledBasisSet>);
} // namespace mpqc


#endif /* MPQC_BASIS_TILEDBASISSET_HPP */
