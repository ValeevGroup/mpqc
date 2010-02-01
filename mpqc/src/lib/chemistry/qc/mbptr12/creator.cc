//
// creator.cc
//
// Copyright (C) 2005 Edward Valeev
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

#include <chemistry/qc/mbptr12/creator.h>

using namespace sc;

DistArray4Creator::DistArray4Creator(const Ref<TwoBodyFourCenterMOIntsRuntime>& moints_rtime,
                                     const std::vector<std::string>& tform_keys) :
  RangeCreator<ObjT> (tform_keys.size()), moints_rtime_(moints_rtime),
      tform_keys_(tform_keys) {
}

DistArray4Creator::ObjT DistArray4Creator::operator()() {
  if (!can_create()) return null();

  const std::string& tform_key = tform_keys_.at(ncreated());
  const ObjT& result = moints_rtime_->get(tform_key)->ints_acc();

  return result;
}

////

TwoBodyIntDescrCreator::TwoBodyIntDescrCreator(
                                               const Ref<R12Technology::R12Technology::CorrelationFactor>& corrfactor,
                                               const Ref<Integral>& integral,
                                               bool CorrFunctionInBra,
                                               bool CorrFunctionInKet) :
      RangeCreator<ObjT> (
                          (CorrFunctionInBra ? (CorrFunctionInKet ? corrfactor->nfunctions()
                                                                     * corrfactor->nfunctions()
                                                                  : corrfactor->nfunctions())
                                             : (CorrFunctionInKet ? corrfactor->nfunctions()
                                                                  : 1))),
      corrfactor_(corrfactor), integral_(integral),
      CorrFunctionInBraKet_(CorrFunctionInBra && CorrFunctionInKet),
      nf12bra_((CorrFunctionInBra ? corrfactor->nfunctions() : 1)),
      nf12ket_((CorrFunctionInKet ? corrfactor->nfunctions() : 1)),
      braindex_(0), ketindex_(0) {
}

TwoBodyIntDescrCreator::ObjT TwoBodyIntDescrCreator::operator()() {
  if (!can_create()) return null();

  ObjT result;
  if (CorrFunctionInBraKet_)
    result = corrfactor_->tbintdescr(integral_, braindex_, ketindex_);
  else
    result = corrfactor_->tbintdescr(integral_, braindex_);

  increment_indices();
  return result;
}

void TwoBodyIntDescrCreator::increment_indices() {
  if (CorrFunctionInBraKet_) {
    ++ketindex_;
    if (ketindex_ == nf12ket_) {
      ++braindex_;
      ketindex_ = 0;
    }
  } else {
    ++braindex_;
  }
  next();
}

////

R12TwoBodyIntKeyCreator::R12TwoBodyIntKeyCreator(const Ref<TwoBodyFourCenterMOIntsRuntime>& moints_rtime,
                                                 const Ref<OrbitalSpace>& bra1,
                                                 const Ref<OrbitalSpace>& ket1,
                                                 const Ref<OrbitalSpace>& bra2,
                                                 const Ref<OrbitalSpace>& ket2,
                                                 const Ref<R12Technology::CorrelationFactor>& corrfactor,
                                                 bool CorrFunctionInBra,
                                                 bool CorrFunctionInKet,
                                                 std::string layout) :
      RangeCreator<ObjT> (
                          (CorrFunctionInBra ? (CorrFunctionInKet ? corrfactor->nfunctions()
                                                                     * corrfactor->nfunctions()
                                                                  : corrfactor->nfunctions())
                                             : (CorrFunctionInKet ? corrfactor->nfunctions()
                                                                  : 1))),
      corrfactor_(corrfactor), layout_key_(layout), moints_rtime_(moints_rtime),
      bra1_(bra1), bra2_(bra2), ket1_(ket1), ket2_(ket2),
      CorrFunctionInBra_(CorrFunctionInBra), CorrFunctionInKet_(CorrFunctionInKet),
      nf12bra_((CorrFunctionInBra ? corrfactor->nfunctions() : 1)),
      nf12ket_((CorrFunctionInKet ? corrfactor->nfunctions() : 1)),
      braindex_(0), ketindex_(0)
{
  if (! (CorrFunctionInBra_ || CorrFunctionInKet_) ) {
    throw ProgrammingError("R12TwoBodyIntKeyCreator used but no correlation factor is present on bra or ket",__FILE__,__LINE__);
  }
}

R12TwoBodyIntKeyCreator::ObjT R12TwoBodyIntKeyCreator::null() const { return std::string(""); }

R12TwoBodyIntKeyCreator::ObjT R12TwoBodyIntKeyCreator::operator()()
{
  if (!can_create()) return null();

  const Ref<Integral>& integral = moints_rtime_->factory()->integral();
  Ref<TwoBodyIntDescr> descr;
  if (CorrFunctionInBra_ && CorrFunctionInKet_)
    descr = corrfactor_->tbintdescr(integral,braindex_,ketindex_);
  else if (CorrFunctionInBra_)
    descr = corrfactor_->tbintdescr(integral,braindex_);
  else // if (CorrFunctionInKet_)
    descr = corrfactor_->tbintdescr(integral,ketindex_);
  const std::string descr_key = moints_rtime_->descr_key(descr);
  ObjT result = ParsedTwoBodyFourCenterIntKey::key(bra1_->id(),bra2_->id(),ket1_->id(),ket2_->id(),descr_key,layout_key_);

  increment_indices();
  return result;
}

void R12TwoBodyIntKeyCreator::increment_indices() {
  if (CorrFunctionInBra_ && CorrFunctionInKet_) {
    ++ketindex_;
    if (ketindex_ == nf12ket_) {
      ++braindex_;
      ketindex_ = 0;
    }
  }
  else if (CorrFunctionInBra_) {
    ++braindex_;
  }
  else if (CorrFunctionInKet_) {
    ++ketindex_;
  }
  next();
}

