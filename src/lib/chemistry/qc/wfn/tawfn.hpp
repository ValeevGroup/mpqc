//
// tawfn.hpp
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

#ifndef CHEMISTRY_WFN_TAWFN_HPP
#define CHEMISTRY_WFN_TAWFN_HPP

#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/tiledbasisset.hpp>
#include <chemistry/molecule/energy.h>
#include <tiled_array.h>

namespace mpqc {
    /// @add to group TAWFN
    /// @{
    
    /** The TiledArrayWavefunction class inherits from the MolecularEnergy class.
     * It will compute energies using TiledArrayBasisSet which is a specialization
     * of GaussianBasisSet
     */
    class TiledArrayWavefunction: public sc::MolecularEnergy {
    private:
        sc::RefSCDimension aodim_;
        sc::RefSCDimension sodim_;
        sc::Ref<TiledBasisSet> tbs_;
        sc::Ref<sc::Integral> integral_;
        sc::ResultRefSymmSCMatrix overlap_;
        sc::ResultRefSymmSCMatrix hcore_;

    public:
        // Short hand for TiledArray Matrix
        // Fix made double for short term
        using TAMat = TiledArray::Array<double, 2>; // TiledArray::Array<double,2>;
        //TiledArrayWavefunction(sc::StateIn &s);
        /** The KeyVal constructor.
         *
         */
        TiledArrayWavefunction(const sc::Ref<sc::KeyVal> &kval);
        virtual ~TiledArrayWavefunction() {};

        //void save_data_state(sc::StateOut &s);

        /// total charge of system.
        double total_charge() const;

        /// number of electrons in the system.
        virtual int nelectron() = 0;
        /// Returns the AO density.
        virtual TAMat ao_density();
        /// Returns the AO overlap.
        virtual TAMat ao_overlap();
        /// Return basis set.
        sc::Ref<TiledBasisSet> basis() const;
        
    };

/// @}

}// namespace mpqc

#endif /* CHEMISTRY_WFN_TAWFN_HPP */
