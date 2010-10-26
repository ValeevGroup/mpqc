//
// mpqc/cca.h
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

#ifndef _mpqc_cca_h
#define _mpqc_cca_h

#include <util/misc/ccaenv.h>
#include <MPQC_ComponentFactory.hxx>
#include <ccaffeine_AbstractFramework.hxx>

namespace sc {

/** The MPQC_CCAFramework class handles embedded CCA frameworks. */
class MPQC_CCAFramework: public AbstractCCAFramework {

    ccaffeine::AbstractFramework fw_;
    gov::cca::Services services_;
    gov::cca::ports::BuilderService bs_;
    gov::cca::TypeMap type_map_;
    gov::cca::ComponentID my_id_;
    MPQC::ComponentFactory component_factory_;

  public:
    /// Initialize the framework.
    MPQC_CCAFramework(const std::string &);
    /// Returns pointer to framework
    virtual ccaffeine::AbstractFramework* get_framework();
    /// Returns pointer to Services object
    virtual gov::cca::Services* get_services();
    /// Returns pointer to BuilderService object
    virtual gov::cca::ports::BuilderService* get_builder_service();
    /// Returns pointer to type map
    virtual gov::cca::TypeMap* get_type_map();
    /// Returns pointer to "uber" component's ComponentID
    virtual gov::cca::ComponentID* get_component_id();
};

}

#endif
