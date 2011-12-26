//
// gto.cc
//
// Copyright (C) 2011 Edward Valeev
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

#include <cmath>
#include <algorithm>
#include <scconfig.h>
#include <util/ref/ref.h>
#include <chemistry/qc/libint2/gto.h>
#include <chemistry/qc/libint2/macros.h>

using namespace sc;

Ref<GTOInfo> GTOInfo::instance_ = new GTOInfo;

const Ref<GTOInfo>&
GTOInfo::instance() {
  return instance_;
}

GTOInfo::GTOInfo() : fp1_(1.0, lmax_*lmax_*lmax_)
{
  df_.resize( std::max(lmax_*2u,3u) );
  df_[0] = 1.0;
  df_[1] = 1.0;
  df_[2] = 1.0;
  for(unsigned int i=3u; i<df_.size(); i++){
    df_[i] = (i-1)*df_[i-2];
  }

  norm_.resize(lmax_+1);
  for(unsigned int l=0; l<=lmax_; l++) {
    std::vector<double>& norm_l = norm_[l];

#if INTEGRALLIBINT2_NORMCONV == INTEGRALLIBINT2_NORMCONV_CCA
    const unsigned int nbf_l = INT_NCART_NN(l);
    norm_l.resize(nbf_l);
    std::fill(norm_l.begin(), norm_l.end(), 1.0);
#else
    int i, j, k;
    FOR_CART(i,j,k,l)
      norm_l.push_back( std::sqrt(df_[2*l]/(df_[2*i]*df_[2*j]*df_[2*k])) );
    END_FOR_CART
#endif
  }
}

GTOInfo::~GTOInfo() {
}

const double*
GTOInfo::fp1() {
  return const_cast<const double*>(&(fp1_[0]));
}

const double*
GTOInfo::norm(unsigned int l) {
  return const_cast<const double*>(&(norm_.at(l)[0]));
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
