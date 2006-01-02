//
// creator.cc
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
#pragma implementation
#endif

#include <chemistry/qc/mbptr12/creator.h>

using namespace sc;

NamedTransformCreator::NamedTransformCreator(
  Ref<R12IntEval>& r12eval,
  const Ref<MOIndexSpace>& space1,
  const Ref<MOIndexSpace>& space2,
  const Ref<MOIndexSpace>& space3,
  const Ref<MOIndexSpace>& space4,
  bool CorrFunctionInBra,
  bool CorrFunctionInKet
  ) : RangeCreator<ObjT>((CorrFunctionInBra ? (CorrFunctionInKet ? r12eval->corrfactor()->nfunctions() * r12eval->corrfactor()->nfunctions()
                                                                 : r12eval->corrfactor()->nfunctions())
                                            : (CorrFunctionInKet ? r12eval->corrfactor()->nfunctions()
                                                                 : 1))),
      r12eval_(r12eval), space1_(space1), space2_(space2),
      space3_(space3), space4_(space4),
      CorrFunctionInBraKet_(CorrFunctionInBra && CorrFunctionInKet),
      nf12bra_((CorrFunctionInBra ? r12eval->corrfactor()->nfunctions() : 1)),
      nf12ket_((CorrFunctionInKet ? r12eval->corrfactor()->nfunctions() : 1)),
      braindex_(0),
      ketindex_(0)
  {
  }

NamedTransformCreator::ObjT
NamedTransformCreator::operator()()
{
  if (!can_create())
    return 0;
  
  ObjT result;
  if (CorrFunctionInBraKet_)
    result = r12eval_->get_tform_(r12eval_->transform_label(space1_,space2_,space3_,space4_,braindex_,ketindex_));
  else
    result = r12eval_->get_tform_(r12eval_->transform_label(space1_,space2_,space3_,space4_,braindex_));
  
  increment_indices();
  return result;
}

void
NamedTransformCreator::increment_indices()
{
  if (CorrFunctionInBraKet_) {
    ++ketindex_;
    if (ketindex_ == nf12ket_) {
      ++braindex_;
      ketindex_ = 0;
    }
  }
  else {
    ++braindex_;
  }
  next();
}

////

TwoBodyIntDescrCreator::TwoBodyIntDescrCreator(
  const Ref<CorrelationFactor>& corrfactor,
  const Ref<Integral>& integral,
  bool CorrFunctionInBra,
  bool CorrFunctionInKet
  ) : RangeCreator<ObjT>((CorrFunctionInBra ? (CorrFunctionInKet ? corrfactor->nfunctions() * corrfactor->nfunctions()
                                                                 : corrfactor->nfunctions())
                                            : (CorrFunctionInKet ? corrfactor->nfunctions()
                                                                 : 1))),
      corrfactor_(corrfactor), integral_(integral),
      CorrFunctionInBraKet_(CorrFunctionInBra && CorrFunctionInKet),
      nf12bra_((CorrFunctionInBra ? corrfactor->nfunctions() : 1)),
      nf12ket_((CorrFunctionInKet ? corrfactor->nfunctions() : 1)),
      braindex_(0),
      ketindex_(0)
  {
  }

TwoBodyIntDescrCreator::ObjT
TwoBodyIntDescrCreator::operator()()
{
  if (!can_create())
    return 0;
  
  ObjT result;
  if (CorrFunctionInBraKet_)
    result = corrfactor_->tbintdescr(integral_,braindex_,ketindex_);
  else
    result = corrfactor_->tbintdescr(integral_,braindex_);
  
  increment_indices();
  return result;
}

void
TwoBodyIntDescrCreator::increment_indices()
{
  if (CorrFunctionInBraKet_) {
    ++ketindex_;
    if (ketindex_ == nf12ket_) {
      ++braindex_;
      ketindex_ = 0;
    }
  }
  else {
    ++braindex_;
  }
  next();
}

