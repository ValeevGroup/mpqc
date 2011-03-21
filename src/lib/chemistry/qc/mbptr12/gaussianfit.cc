//
// gaussianfit.cc
//
// Copyright (C) 2007 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
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
#pragma implementation
#endif

#include <chemistry/qc/mbptr12/gaussianfit.h>
#include <chemistry/qc/mbptr12/gaussianfit.timpl.h>

using namespace std;
using namespace sc;

extern "C" {

    void
    __eval_slater(double* params, double* f, int nparam, int np, void *extraparams)
    {
	eval_f<sc::mbptr12::Slater1D,sc::mbptr12::Gaussian1D>(params,f,nparam,np,extraparams);
    }

    void
    __eval_slater_dfdp(double* params, double* f, int nparam, int np, void *extraparams)
    {
	eval_dfdp<sc::mbptr12::Slater1D,sc::mbptr12::Gaussian1D>(params,f,nparam,np,extraparams);
    }

    void
    __eval_slater_pgauss(double* params, double* f, int nparam, int np, void *extraparams)
    {
	eval_f<sc::mbptr12::Slater1D,sc::mbptr12::PowerGaussian1D>(params,f,nparam,np,extraparams);
    }

    void
    __eval_slater_dfdp_pgauss(double* params, double* f, int nparam, int np, void *extraparams)
    {
	eval_dfdp<sc::mbptr12::Slater1D,sc::mbptr12::PowerGaussian1D>(params,f,nparam,np,extraparams);
    }

};

namespace sc {
  namespace mbptr12 {
    template<>
    eval_f_ptr __to_extern_C_eval<sc::mbptr12::Slater1D,sc::mbptr12::Gaussian1D>::f_ptr(__eval_slater);
    template<>
    eval_dfdp_ptr __to_extern_C_eval<sc::mbptr12::Slater1D,sc::mbptr12::Gaussian1D>::dfdp_ptr(__eval_slater_dfdp);
    template<>
    eval_f_ptr __to_extern_C_eval<sc::mbptr12::Slater1D,sc::mbptr12::PowerGaussian1D>::f_ptr(__eval_slater_pgauss);
    template<>
    eval_dfdp_ptr __to_extern_C_eval<sc::mbptr12::Slater1D,sc::mbptr12::PowerGaussian1D>::dfdp_ptr(__eval_slater_dfdp_pgauss);
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
