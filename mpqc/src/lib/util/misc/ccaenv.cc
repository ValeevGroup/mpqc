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

#include <util/misc/ccaenv.h>

using namespace sc;

int CCAEnv::initialized_ = 0;
//ccaffeine::AbstractFramework CCAEnv::fw_ = 
//  ccaffeine::AbstractFramework::AbstractFramework();
ccaffeine::AbstractFramework CCAEnv::fw_;

void 
CCAEnv::init(std::string &args)
{ 
  fw_ = ccaffeine::AbstractFramework::_create();
  fw_.initialize(args); 
  initialized_=1; 
}
int
CCAEnv::initialized() 
{ 
  return initialized_; 
}

gov::cca::AbstractFramework 
CCAEnv::get_framework() 
{ 
  return fw_; 
}

