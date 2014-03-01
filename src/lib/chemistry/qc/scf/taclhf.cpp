//
// tahf.cpp
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

#ifndef _MPQC_CHEMISTRY_QC_SCF_TACLHF_HPP_
#define _MPQC_CHEMISTRY_QC_SCF_TACLHF_HPP_

#include <chemistry/qc/scf/taclscf.hpp>

namespace mpqc{
namespace TA{
class CLHF: public CLSCF {
public:
    CLHF(const Ref<KeyVal>&);
    ~CLHF();

protected:
    void ao_fock(double accuracy);

}; //class CLHF

} // namespace TA
} // namespace mpqc



#endif /* _MPQC_CHEMISTRY_QC_SCF_TACLHF_HPP_ */
