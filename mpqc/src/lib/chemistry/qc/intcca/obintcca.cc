//
// obintcca.cc
//
// Copyright (C) 2004 Sandia National Laboratories
//
// Author: Joe Kenny <jpkenny@sandia.gov>
// Maintainer: Joe Kenny
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

#include <chemistry/qc/intcca/obintcca.h>

using namespace Chemistry::QC::GaussianBasis;
using namespace sc;

////////////////////////////////////////////////////////////////////////////
// OneBodyIntCCA

OneBodyIntCCA::OneBodyIntCCA(Integral* integral,
                           const Ref<GaussianBasisSet>&bs1,
                           const Ref<GaussianBasisSet>&bs2,
                           IntegralFunction ifunc,
			   IntegralEvaluatorFactory eval_factory,
                           bool use_opaque ):
  OneBodyInt(integral,bs1,bs2), eval_factory_(eval_factory), 
  use_opaque_(use_opaque) 
{
  int1ecca_ = new Int1eCCA(integral,bs1,bs2,0,eval_factory,use_opaque);
  intfunc_ = ifunc;
  buffer_ = int1ecca_->buffer();
}

OneBodyIntCCA::~OneBodyIntCCA()
{
}

void
OneBodyIntCCA::compute_shell(int i, int j)
{
  //ExEnv::out0() << "OneBodyIntCCA::compute_shell() called\n";
  (int1ecca_.pointer()->*intfunc_)(i, j);
}

bool
OneBodyIntCCA::cloneable()
{
  return true;
}

Ref<OneBodyInt>
OneBodyIntCCA::clone()
{
  return new OneBodyIntCCA(integral_, bs1_, bs2_, intfunc_, 
                           eval_factory_, use_opaque_ );
}

// ////////////////////////////////////////////////////////////////////////////
// // PointChargeIntV3

// PointChargeIntV3::PointChargeIntV3(
//     Integral *integral,
//     const Ref<GaussianBasisSet>&bs1,
//     const Ref<GaussianBasisSet>&bs2,
//     const Ref<PointChargeData>&dat):
//   OneBodyInt(integral,bs1,bs2),
//   data_(dat)
// {
//   int1ev3_ = new Int1eV3(integral,bs1,bs2,0);
//   buffer_ = int1ev3_->buffer();
// }

// PointChargeIntV3::~PointChargeIntV3()
// {
// }

// void
// PointChargeIntV3::compute_shell(int i,int j)
// {
//   int1ev3_->point_charge(i,j,
//                          data_->ncharges(),
//                          data_->charges(),
//                          data_->positions());
// }

// ////////////////////////////////////////////////////////////////////////////
// // EfieldDotVectorIntV3

// EfieldDotVectorIntV3::EfieldDotVectorIntV3(
//     Integral *integral,
//     const Ref<GaussianBasisSet>&bs1,
//     const Ref<GaussianBasisSet>&bs2,
//     const Ref<EfieldDotVectorData>&dat) :
//   OneBodyInt(integral,bs1,bs2),
//   data_(dat)
// {
//   int1ev3_ = new Int1eV3(integral,bs1,bs2,0);
//   buffer_ = int1ev3_->buffer();
// }

// EfieldDotVectorIntV3::~EfieldDotVectorIntV3()
// {
// }

// void
// EfieldDotVectorIntV3::compute_shell(int i,int j)
// {
//   int nbfi = basis1()->shell(i).nfunction();
//   int nbfj = basis2()->shell(j).nfunction();
//   int nint = nbfi*nbfj;
//   double *tmp;
//   int ii,jj;

//   int1ev3_->efield(i,j,data_->position);

//   tmp = int1ev3_->buffer();
//   for (ii=0; ii<nint; ii++) {
//       double tmpval = 0.0;
//       for (jj=0; jj<3; jj++) {
//           tmpval += *tmp++ * data_->vector[jj];
//         }
//       buffer_[ii] = tmpval;
//     }
// }

// ////////////////////////////////////////////////////////////////////////////
// // DipoleIntV3

// DipoleIntV3::DipoleIntV3(Integral *integral,
//                          const Ref<GaussianBasisSet>&bs1,
//                          const Ref<GaussianBasisSet>&bs2,
//                          const Ref<DipoleData>&dat) :
//   OneBodyInt(integral,bs1,bs2),
//   data_(dat)
// {
//   int1ev3_ = new Int1eV3(integral,bs1,bs2,0);
//   buffer_ = int1ev3_->buffer();
//   if (data_.null()) {
//       data_ = new DipoleData;
//     }
// }

// DipoleIntV3::~DipoleIntV3()
// {
// }

// void
// DipoleIntV3::compute_shell(int i,int j)
// {
//   int1ev3_->dipole(i,j,data_->origin);
// }

// ////////////////////////////////////////////////////////////////////////////
// // OneBodyDerivIntV3

// OneBodyDerivIntV3::OneBodyDerivIntV3(Integral *integral,
//                                      const Ref<GaussianBasisSet>&bs1,
//                                      const Ref<GaussianBasisSet>&bs2,
//                                      IntegralFunction ifunc):
//   OneBodyDerivInt(integral,bs1,bs2)
// {
//   int1ev3_ = new Int1eV3(integral,bs1,bs2,1);
//   intfunc_ = ifunc;
//   buffer_ = int1ev3_->buffer();
// }

// OneBodyDerivIntV3::~OneBodyDerivIntV3()
// {
// }

// void
// OneBodyDerivIntV3::compute_shell(int i, int j, DerivCenters& c)
// {
//   (int1ev3_.pointer()->*intfunc_)(i,j,0,basis1()->shell_to_center(i));
//   c.clear();
//   c.add_center(0,basis1(),i);
//   c.add_omitted(1,basis2(),j);
// }

// void
// OneBodyDerivIntV3::compute_shell(int i, int j, int c)
// {
//   (int1ev3_.pointer()->*intfunc_)(i,j,0,c);
// }

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
