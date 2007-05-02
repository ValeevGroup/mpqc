//
// mpqc/cca.cc
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

#include <scdirlist.h>
#include <cca.h>

using namespace sc;

MPQC_CCAFramework::MPQC_CCAFramework(const std::string &args)
{ 
  fw_ = ccaffeine::AbstractFramework::_create();
  fw_.initialize(args,0,0); 
  type_map_ = fw_.createTypeMap();
  services_ = fw_.getServices("uber","UberComponent",type_map_);
  my_id_    = services_.getComponentID();
  services_.registerUsesPort("bs","gov.cca.BuilderService",type_map_);
  bs_ = babel_cast<gov::cca::ports::BuilderService>( services_.getPort("bs") );

  component_factory_ = MPQC::ComponentFactory::_create();
  component_factory_.addDescription("Chemistry.IntegralSuperFactory",
                                    "Chemistry.IntegralSuperFactory");
  component_factory_.addDescription("MPQC.IntV3EvaluatorFactory",
                                    "MPQC.IntV3EvaluatorFactory");
#ifdef HAVE_CINTS
  component_factory_.addDescription("MPQC.CintsEvaluatorFactory",
                                    "MPQC.CintsEvaluatorFactory");
#endif
  services_.addProvidesPort(component_factory_, 
                            "ccaffeine.ports.ComponentFactory",
                            "ccaffeine.ports.ComponentFactory",
                            type_map_);
}

ccaffeine::AbstractFramework* 
MPQC_CCAFramework::get_framework() 
{ return &fw_; }

gov::cca::Services* 
MPQC_CCAFramework::get_services()
{ return &services_; }

gov::cca::ports::BuilderService* 
MPQC_CCAFramework::get_builder_service()
{ return &bs_; }

gov::cca::TypeMap* 
MPQC_CCAFramework::get_type_map()
{ return &type_map_; }

gov::cca::ComponentID* 
MPQC_CCAFramework::get_component_id()
{ return &my_id_; }
