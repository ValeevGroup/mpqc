//
// linearr12.timpl.h
//
// Copyright (C) 2007 Edward Valeev
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

#ifndef _chemistry_qc_mbptr12_linearr12timpl_h
#define _chemistry_qc_mbptr12_linearr12timpl_h

#include <cmath>
#include <algorithm>
#include <chemistry/qc/mbptr12/gaussianfit.h>
#include <chemistry/qc/mbptr12/linearr12.h>

namespace sc {

    namespace LinearR12 {

      template <class CorrFactor, class Fitter>
        Ref<CorrelationFactor> stg_to_g12(const Fitter& fitter, double gamma) {

	  using sc::mbptr12::Slater1D;
	  typedef typename Fitter::Gaussians Gaussians;
	  Slater1D stg(gamma);
	  Gaussians gtgs = fitter(stg);

	  // feed to the constructor of CorrFactor
	  typedef IntParamsG12::PrimitiveGeminal PrimitiveGeminal;
	  typedef IntParamsG12::ContractedGeminal ContractedGeminal;
	  ContractedGeminal geminal;
	  typedef typename Gaussians::const_iterator citer;
	  for(citer g=gtgs.begin(); g!=gtgs.end(); ++g) {
	      geminal.push_back(*g);
	  }
	  std::vector<ContractedGeminal> geminals(1,geminal);

	  Ref<CorrelationFactor> cf = new CorrFactor(geminals);
	  return cf;
      }

      template <class Fitter>
      Ref<CorrelationFactor> angstg_to_geng12(const Fitter& fitter, double alpha, double gamma) {

	  const double halfalpha = alpha/2.0;

	  using sc::mbptr12::Slater1D;
	  typedef typename Fitter::Gaussians Gaussians;
	  Slater1D stg(gamma);
	  Gaussians gtgs = fitter(stg);

	  // feed to the constructor of CorrFactor
	  typedef IntParamsGenG12::PrimitiveGeminal PrimitiveGeminal;
	  typedef IntParamsGenG12::ContractedGeminal ContractedGeminal;
	  ContractedGeminal geminal;
	  typedef typename Gaussians::const_iterator citer;
	  for(citer g=gtgs.begin(); g!=gtgs.end(); ++g) {
	      const double alpha_i = halfalpha;
	      const double gamma_i = (*g).first - halfalpha;
	      const double C_i = (*g).second;
	      // see basis/intparams.h
	      PrimitiveGeminal i = std::make_pair(std::make_pair(alpha_i,gamma_i),C_i);
	      geminal.push_back(i);
	  }
	  std::vector<ContractedGeminal> geminals(1,geminal);

	  Ref<CorrelationFactor> cf = new GenG12CorrelationFactor(geminals);
	  return cf;
      }

      template <class IntParam>
      bool
      CorrParamCompare<IntParam>::equiv(const ContractedGeminals& A, const ContractedGeminals& B)
      {
	  unsigned int nf = A.size();
	  if (nf != B.size()) return false;

	  for(unsigned int f=0; f<nf; ++f) {
	      const ContractedGeminal& Af = A[f];
	      const ContractedGeminal& Bf = B[f];
	      unsigned int np = Af.size();
	      if (np != Bf.size()) return false;
	      for(unsigned int p=0; p<np; ++p) {
		  if (!equiv(Af[p],Bf[p])) return false;
	      }
	  }

	  return true;
      }
	
    };

};
    
#endif // include guards
    
/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
