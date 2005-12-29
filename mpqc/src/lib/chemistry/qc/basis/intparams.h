//
// intparams.h
//
// Copyright (C) 2005 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
// Maintainer: EV
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUG__
#pragma interface
#endif

#include <vector>
#include <utility>
#include <util/ref/ref.h>

#ifndef _chemistry_qc_basis_intparams_h
#define _chemistry_qc_basis_intparams_h

namespace sc {
  
  /** This class passes optional operator parameters. These parameters can be
      passed to
      1) Factory methods which initialize XXXBodyInt objects
      2) compute_shell methods of XXXBodyInt objects
  */
  class IntParams : public RefCount {
    public:
      IntParams(unsigned int nparams=0);
      virtual ~IntParams();
      
      unsigned int nparams() const;
      
    private:
      unsigned int nparams_;
  };

  /** Passes params to Integral::g12() */
  class IntParamsG12 : public IntParams {
    public:
      typedef std::pair<double,double> PrimitiveGeminal;
      typedef std::vector<PrimitiveGeminal> ContractedGeminal;
      static ContractedGeminal zero_exponent_geminal;
      
      IntParamsG12(const ContractedGeminal& bra,
                   const ContractedGeminal& ket);
      ~IntParamsG12();

      const ContractedGeminal& bra() const;
      const ContractedGeminal& ket() const;
      
    private:
      ContractedGeminal bra_;
      ContractedGeminal ket_;
  };
  

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
