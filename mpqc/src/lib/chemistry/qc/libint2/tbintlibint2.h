//
// tbintlibint2.h
//
// Copyright (C) 2001 Edward Valeev
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

#ifndef _chemistry_qc_libint2_tbint_h
#define _chemistry_qc_libint2_tbint_h

#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/libint2/int2e.h>

namespace sc {

// Types of available 2-electron integral evaluators
typedef enum {erieval = 0, grteval = 1, g12eval = 2} tbinteval;

/** This implements electron repulsion integrals in the IntLibint2 library. */
class TwoBodyIntLibint2 : public TwoBodyInt {

    unsigned int num_tbint_types_;

  protected:
    Ref<Int2eLibint2> int2elibint2_;
  
  public:
    typedef IntParamsG12::ContractedGeminal ContractedGeminal;
    TwoBodyIntLibint2(Integral*integral,
                 const Ref<GaussianBasisSet>&b1,
                 const Ref<GaussianBasisSet>&b2,
                 const Ref<GaussianBasisSet>&b3,
                 const Ref<GaussianBasisSet>&b4,
                 size_t storage, tbinteval int2etype,
                 const ContractedGeminal& gamma_bra,
                 const ContractedGeminal& gamma_ket = IntParamsG12::null_geminal);
    ~TwoBodyIntLibint2();

    unsigned int num_tbint_types() const {
      return num_tbint_types_;
    }

    int log2_shell_bound(int,int,int,int);
    void compute_shell(int,int,int,int);

    size_t used_storage() const { return int2elibint2_->storage_used(); }
    void set_integral_storage(size_t storage);

    const double *buffer(tbint_type te_type) const {
      return int2elibint2_->buffer(te_type);
    }
};

/** This implements electron repulsion derivative integrals in the IntV3
    library. */
class TwoBodyDerivIntLibint2 : public TwoBodyDerivInt {
  protected:
    Ref<Int2eLibint2> int2elibint2_;

  public:
    TwoBodyDerivIntLibint2(Integral*integral,
                      const Ref<GaussianBasisSet>&b1,
                      const Ref<GaussianBasisSet>&b2,
                      const Ref<GaussianBasisSet>&b3,
                      const Ref<GaussianBasisSet>&b4,
                      size_t storage, tbinteval int2etype);
    ~TwoBodyDerivIntLibint2();

    int log2_shell_bound(int,int,int,int);
    void compute_shell(int,int,int,int,DerivCenters&);

    size_t used_storage() const { return int2elibint2_->storage_used(); }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
