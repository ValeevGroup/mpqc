//
// ccaenv.h
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

#ifdef __GNUG__
#pragma interface
#endif

#ifndef _util_misc_ccaenv_h
#define _util_misc_ccaenv_h

#include <ccaffeine_AbstractFramework.hh>
#include <gov_cca.hh>
#include <MPQC_ComponentFactory.hh>

namespace sc {

/** The CCAEnv class handles embedded CCA frameworks. */
class CCAEnv {
  private:
    static int initialized_;
    static ccaffeine::AbstractFramework fw_;
    static gov::cca::Services services_;
    static gov::cca::ports::BuilderService bs_;
    static gov::cca::TypeMap type_map_;
    static gov::cca::ComponentID my_id_;
    static MPQC::ComponentFactory component_factory_;

  public:
    /// Initialize the framework
    static void init(std::string &args);
    /// Return nonzero if CCAEnv has been initialized.
    static int initialized();
    /// Returns pointer to framework
    static ccaffeine::AbstractFramework* get_framework();
    /// Returns pointer to Services object
    static gov::cca::Services* get_services();
    /// Returns pointer to BuilderService object
    static gov::cca::ports::BuilderService* get_builder_service();
    /// Returns pointer to type map
    static gov::cca::TypeMap* get_type_map(); 
    /// Returns pointer to "uber" component's ComponentID
    static gov::cca::ComponentID* get_component_id();
};

}

#endif
