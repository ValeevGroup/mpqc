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

#ifndef _MPQC_CHEMISTRY_QC_SCF_TASCF_HPP_
#define _MPQC_CHEMISTRY_QC_SCF_TASCF_HPP_

#include <chemistry/qc/wfn/tawfn.hpp>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/tbint.h>
#include <TiledArray/algebra/diis.h>

namespace mpqc{
  namespace TA {

    class SCF : public Wavefunction {
    public:
        typedef Wavefunction::Matrix Matrix;

        SCF(const sc::Ref<sc::KeyVal> &kval);
        virtual ~SCF();
        virtual void compute() override;

        virtual const Matrix& fock(); // Final Fock matrix
        virtual const Matrix& rdm1(sc::SpinCase1) override;
        virtual const Matrix& hcore();
        virtual double scf_energy() = 0;

        /// @return the number of electrons in the system
        virtual size_t nelectron() const override;

        virtual size_t occupation(){return occupation_;}

    protected:
        ::TiledArray::DIIS<Matrix> diis;
        virtual Matrix& scf_fock(){return fock_.result_noupdate();}
        size_t & occ(){return occupation_;}
        virtual double iter_energy() = 0;

    private:
        // Number of iterations to use
        unsigned int maxiter_;
        unsigned int miniter_;

        // Number of electrons
        size_t occupation_;

        // Fock Matrix
        ResultMatrix fock_;

        // Hcore
        Matrix hcore_;

        static sc::ClassDesc class_desc_;

    };
  } // namespace mpqc::TA
} //namespace mpqc


#endif /* _MPQC_CHEMISTRY_QC_SCF_TASCF_HPP_ */
