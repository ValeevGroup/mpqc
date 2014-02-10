//
// tascf.hpp
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

#ifndef CHEMISTRY_QC_SCF_TASCF_HPP
#define CHEMISTRY_QC_SCF_TASCF_HPP

#include <chemistry/qc/wfn/tawfn.hpp>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/tbint.h>

namespace mpqc{
    class TiledArrayScf : public TiledArrayWavefunction {
    public:
        typedef TiledArrayWavefunction::TAMat TAMat;

    protected:
        // Number of iterations to use
        unsigned int maxiter_ = 1000;
        unsigned int miniter_ = 0;

        // Integral objects
        sc::Ref<sc::TwoBodyInt> *tbints_;

    public:
        TiledArrayScf(const sc::Ref<sc::KeyVal> &kval);
        virtual ~TiledArrayScf();
        virtual void compute() override;
        virtual int nelectron() override;

        virtual TAMat ao_fock();
        virtual TAMat ao_density() override;
        virtual TAMat ao_overlap() override;
        virtual double scf_energy();

    };

} //namespace mpqc


#endif /* CHEMISTRY_QC_SCF_TASCF_HPP */
