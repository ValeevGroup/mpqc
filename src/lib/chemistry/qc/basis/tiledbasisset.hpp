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

#ifndef CHEMISTRY_QC_BASIS_TILEDBASISSET_HPP
#define CHEMISTRY_QC_BASIS_TILEDBASISSET_HPP

#include <tiled_array.h>

#include "basis.h"
#include "gaussbas.h"
#include <chemistry/molecule/atom.h>
#include <chemistry/molecule/molecule.h>

#include <mpqc/utility/foreach.hpp>
#include "shellorder.hpp"
#include "kcluster.hpp"

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
        TiledBasisSet(const sc::Ref<sc::KeyVal> &keyval);

        static std::string name_conv_TBS(const std::string& name);

        /**
         * Returns a TiledArray::TiledRange1 for tiles made by k-means clustering.
         */
        TiledArray::TiledRange1 trange1() const;

        TiledBasisSet(sc::StateIn& s);

        void save_data_state(sc::StateOut& s);

    private:
        std::size_t nclusters_;
        ShellRange SRange_;
    };

} // namespace mpqc


#endif /* CHEMISTRY_QC_BASIS_TILEDBASISSET_HPP */
