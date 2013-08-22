//
// k_means_basis.hpp
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

#ifndef MPQC_BASIS_KMEANSBASIS_HPP
#define MPQC_BASIS_KMEANSBASIS_HPP

#include <tiled_array.h>
#include <vector>

#include <chemistry/qc/basis/basis.h>
#include <chemistry/molecule/atom.h>
#include <chemistry/molecule/molecule.h>

#include <Eigen/Dense>
#include "kcluster.hpp"
#include <string>

namespace mpqc{
namespace basis{

    /**
     * Will create a basis using k-means clustering.
     */
    class TiledBasisSet : public sc::GaussianBasisSet {

    public:
        typedef KCluster::Atom Atom;
        typedef KCluster::Vector3 Vector3;
        typedef sc::Ref<sc::Molecule> Mol;
        typedef sc::Ref<sc::GaussianBasisSet> Basis;
        typedef std::string string;

        TiledBasisSet(const sc::Ref<sc::Keyval> &keyval){
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
            // Do stuff
        }


    private:
        RefMol mol_;
        std::size_t ncenters_;
        Vector3 mol_com_;
        std::vector<KCluster> kclusters_;
    };

} // namespace basis
} // namespace mpqc


#endif /* MPQC_BASIS_KMEANSBASIS_HPP */
