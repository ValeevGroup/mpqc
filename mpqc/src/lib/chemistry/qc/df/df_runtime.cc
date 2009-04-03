//
// df_runtime.cc
//
// Copyright (C) 2009 Edward Valeev
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

#include <chemistry/qc/df/df_runtime.h>

using namespace sc;

/////////////////////////////////////////////////////////////////////////////

namespace {
  // pop off str from beginning up to token.
  std::string
  pop_till_token(std::string& str,
                 char token) {
    const size_t next_token_pos = str.find_first_of(token);
    std::string result;
    if (next_token_pos != std::string::npos) {
      result = str.substr(0,next_token_pos);
      str.erase(0,next_token_pos+1);
    }
    else {
      result = str;
      str.clear();
    }
    return result;
  }
}

ParsedDensityFittingKey::ParsedDensityFittingKey(const std::string& key) :
  key_(key)
{
  std::string keycopy(key);

  // pop off the leading '('
  assert(keycopy[0] == '(');
  keycopy.erase(keycopy.begin());
  // get space1
  space1_ = pop_till_token(keycopy,' ');
  // get space2
  space2_ = pop_till_token(keycopy,'|');
  // get rid of the "DF" token
  std::string crap = pop_till_token(keycopy,'|');
  // get fspace
  fspace_ = pop_till_token(keycopy,')');

#if 0
  ExEnv::out0() << indent << "ParsedDensityFittingKey::ParsedDensityFittingKey():" << std::endl << incindent;
  ExEnv::out0() << indent << "key = " << key_ << std::endl;
  ExEnv::out0() << indent << "space1 = " << space1_ << std::endl;
  ExEnv::out0() << indent << "space2 = " << space2_ << std::endl;
  ExEnv::out0() << indent << "fspace = " << fspace_ << std::endl;
#endif
}

std::string
ParsedDensityFittingKey::key(const std::string& space1,
                             const std::string& space2,
                             const std::string& fspace)
{
  std::ostringstream oss;
  oss << "(" << space1 << " " << space2 << "|DF|" << fspace << ")";
  return oss.str();
}

/////////////////////////////////////////////////////////////////////////////

ClassDesc
DensityFittingRuntime::class_desc_(typeid(this_type),
                                   "DensityFittingRuntime",
                                   1,
                                   "virtual public SavableState", 0, 0,
                                   create<this_type> );

DensityFittingRuntime::DensityFittingRuntime(const Ref<MOIntsRuntime>& r) :
  moints_runtime_(r),
  results_(ResultRegistry::instance())
{
}

DensityFittingRuntime::DensityFittingRuntime(StateIn& si)
{
  moints_runtime_ << SavableState::restore_state(si);
  results_ = ResultRegistry::restore_instance(si);
}

void
DensityFittingRuntime::save_data_state(StateOut& so)
{
  SavableState::save_state(moints_runtime_.pointer(),so);
  ResultRegistry::save_instance(results_,so);
}

bool
DensityFittingRuntime::exists(const std::string& key) const
{
  return results_->key_exists(key);
}


DensityFittingRuntime::ResultRef
DensityFittingRuntime::get(const std::string& key)
{
  if (results_->key_exists(key)) {
    return results_->value(key);
  }
  else {  // if not found
    try { ParsedResultKey parsedkey(key); }
    catch (...) {
      std::ostringstream oss;
      oss << "DensityFittingRuntime::get() -- key " << key << " does not match the format";
      throw ProgrammingError(oss.str().c_str(),__FILE__,__LINE__);
    }
    // then create evaluated tform
    const ResultRef& result = create_result(key);
    return result;
  }
  assert(false); // unreachable
}

const DensityFittingRuntime::ResultRef&
DensityFittingRuntime::create_result(const std::string& key)
{
  // parse the key
  ParsedResultKey pkey(key);
  const std::string& space1_str = pkey.space1();
  const std::string& space2_str = pkey.space2();
  const std::string& fspace_str = pkey.fspace();

  // get the spaces and construct the descriptor
  Ref<OrbitalSpaceRegistry> idxreg = OrbitalSpaceRegistry::instance();
  Ref<OrbitalSpace> space1 = idxreg->value(space1_str);
  Ref<OrbitalSpace> space2 = idxreg->value(space2_str);
  Ref<OrbitalSpace> fspace = idxreg->value(fspace_str);

  // create the result
  Ref<DensityFitting> df = new DensityFitting(moints_runtime_, "1/r_{12}",
                                              space1, space2, fspace->basis());
  df->compute();
  ResultRef result = df->C();

  // add to the map
  results_->add(key,result);

  return results_->value(key);
}

/////////////////////////////////////////////////////////////////////////////

ClassDesc
DensityFittingInfo::class_desc_(typeid(this_type),
                                "DensityFittingInfo",
                                1,
                                "virtual public SavableState", 0, 0,
                                create<this_type> );

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
