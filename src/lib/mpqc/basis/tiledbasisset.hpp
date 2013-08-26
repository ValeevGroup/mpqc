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
#include <chemistry/molecule/atom.h>
#include <chemistry/molecule/molecule.h>

#include <Eigen/Dense>
#include <mpqc/utility/foreach.hpp>
#include "kcluster.hpp"
#include <string>

namespace mpqc{
    static sc::ClassDesc TiledBasisSet_cd(
                    typeid(TiledBasisSet), "TiledBasisSet", 1,
                    "public GaussianBasisSet",
                    0,  create<TiledBasisSet>, create<TiledBasisSet>);

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
        TiledBasisSet(const sc::Ref<sc::Keyval> &keyval, std::size_t nclusters):
            nclusters_(nclusters), SRange()
        {
            Basis basis;
            basis << keyval->describedclassvalue("basis");
            if(basis.null()){
                basis = new GaussianBasisSet(keyval);
                if(basis.null()){
                    throw InputError("Could not construct a GaussianBasisSet",
                                     __FILE__, __LINE__,
                                     "basis", 0, class_desc());
                }
            }

            basis::ShellOrder ordering(basis);
            std::vector<Shell> shells = ordering.ordered_shells(nclusters_);
            SRange = ordering.shell_ranges();

            init(name_conv_TBS(basis->name()),
                 name_conv_TBS(basis_-label()),
                 basis_->molecule(),
                 shells);
        }

        static
        std::string
        name_conv_TBS(const std::string& name){
            if(name.empty()) return name;
            std::string new_name = "TiledArray(";
            newname += name;
            newname += ")";
            return newname;
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
               tilesizes.push_back(function_to_shell(SRange[i]));
            }

            // Get the last function in the shell since our loop doesn't cover it.
            tilesizes.push_back(nbasis());

            return TiledArray::TiledRange1(tilesizes.begin(), tilesizes.end());

        }

    private:
        std::size_t nclusters_;
        ShellRange SRange_;
    };
} // namespace mpqc


#endif /* MPQC_BASIS_TILEDBASISSET_HPP */
