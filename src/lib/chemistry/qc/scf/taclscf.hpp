//
// taclscf.hpp
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

#ifndef _MPQC_CHEMISTRY_QC_SCF_TACLSCF_HPP_
#define _MPQC_CHEMISTRY_QC_SCF_TACLSCF_HPP_

#include <chemsitry/qc/scf/tascf.hpp>

namespace mpqc{
namespace TA {

    /** The taclscf class is the base class for implementing self-consistent
     * proceedure for closed-shell molecules in MPQC3.
     */
    class CLSCF : public SCF {
    public :
        typedef SCF::Matrix Matrix;

        /** Key Value constructor . . . */
        CLSCF(const sc::Ref<sc::KeyVal> &kval);
        virtual ~CLSCF();

        virtual void compute() override;
        virtual const Matrix& fock() override;

        virtual const Matrix& ao_density() override;

    private:
        static sc::ClassDesc class_desc_;

    }; // Class CLSCF

} // namespace TA
} // namespace mpqc

#endif /* _MPQC_CHEMISTRY_QC_SCF_TACLSCF_HPP_ */
