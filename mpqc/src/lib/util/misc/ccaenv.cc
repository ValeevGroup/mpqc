//
// ccaenv.cc
//
// Copyright (C) 1997 Limit Point Systems, Inc.
//
// Author: Joseph Kenny <jpkenny@sandia.gov>
// Maintainer: JK
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/misc/ccaenv.h>

using namespace sc;

Ref<CCAFramework> CCAEnv::ccafw_;

void 
CCAEnv::init(const Ref<CCAFramework> &ccafw)
{ 
  ccafw_ = ccafw;
}

int
CCAEnv::initialized() 
{ return ccafw_.nonnull(); }

ccaffeine::AbstractFramework* 
CCAEnv::get_framework()
{ return ccafw_->get_framework(); }

gov::cca::Services* 
CCAEnv::get_services()
{ return ccafw_->get_services(); }

gov::cca::ports::BuilderService* 
CCAEnv::get_builder_service()
{ return ccafw_->get_builder_service(); }

gov::cca::TypeMap* 
CCAEnv::get_type_map()
{ return ccafw_->get_type_map(); }

gov::cca::ComponentID* 
CCAEnv::get_component_id()
{ return ccafw_->get_component_id(); }

CCAFramework::~CCAFramework()
{
}

