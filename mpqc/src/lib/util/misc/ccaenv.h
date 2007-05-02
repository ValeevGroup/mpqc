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

#include <util/ref/ref.h>
//#include <gov_cca_AbstractFramework.hxx>
#include <ccaffeine_AbstractFramework.hxx>
#include <gov_cca.hxx>

namespace sc {

/** The CCAFramework class defines an interface to CCA frameworks. */
class CCAFramework: public RefCount {
  public:
    virtual ~CCAFramework();
    /// Returns pointer to Services object
    virtual gov::cca::Services* get_services() = 0;
    /// Returns pointer to BuilderService object
    virtual gov::cca::ports::BuilderService* get_builder_service() = 0;
    /// Returns pointer to type map
    virtual gov::cca::TypeMap* get_type_map() = 0;
    /// Returns pointer to framework's ComponentID
    virtual gov::cca::ComponentID* get_component_id() = 0;
};

/** The AbstractCCAFramework class defines an interface to
    abstract CCA frameworks. */
class AbstractCCAFramework: public CCAFramework {
  public:
    virtual ~AbstractCCAFramework(){ }
    /// Returns pointer to abstract framework
    virtual ccaffeine::AbstractFramework* get_framework() = 0;
};

/** The Ext_CCAFramework class handles externally initialized 
    CCA frameworks. */
class Ext_CCAFramework: public CCAFramework {

    gov::cca::Services services_;
    gov::cca::ports::BuilderService bs_;
    gov::cca::TypeMap type_map_;
    gov::cca::ComponentID my_id_;

  public:
    /// Initialize the framework.
    Ext_CCAFramework(gov::cca::Services &services);
    /// Returns pointer to Services object
    virtual gov::cca::Services* get_services();
    /// Returns pointer to BuilderService object
    virtual gov::cca::ports::BuilderService* get_builder_service();
    /// Returns pointer to type map
    virtual gov::cca::TypeMap* get_type_map();
    /// Returns pointer to framework's ComponentID
    virtual gov::cca::ComponentID* get_component_id();
};

/** The CCAEnv class provides a CCA environment. */
class CCAEnv {
  private:
    static Ref<CCAFramework> ccafw_;

  public:
    /// Set the global CCA framework.
    static void init(const Ref<CCAFramework> &);
    /// Return nonzero if CCAEnv has been initialized.
    static int initialized();
    /// Returns the CCAFramework object.
    static Ref<CCAFramework> ccaframework() { return ccafw_; }
    /// Returns pointer to Services object
    static gov::cca::Services* get_services();
    /// Returns pointer to BuilderService object
    static gov::cca::ports::BuilderService* get_builder_service();
    /// Returns pointer to type map
    static gov::cca::TypeMap* get_type_map(); 
    /// Returns pointer to framework's ComponentID
    static gov::cca::ComponentID* get_component_id();
};

}

#endif
