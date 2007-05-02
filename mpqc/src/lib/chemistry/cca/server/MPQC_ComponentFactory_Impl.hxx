// 
// File:          MPQC_ComponentFactory_Impl.hxx
// Symbol:        MPQC.ComponentFactory-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.ComponentFactory
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 

#ifndef included_MPQC_ComponentFactory_Impl_hxx
#define included_MPQC_ComponentFactory_Impl_hxx

#ifndef included_sidl_cxx_hxx
#include "sidl_cxx.hxx"
#endif
#ifndef included_MPQC_ComponentFactory_IOR_h
#include "MPQC_ComponentFactory_IOR.h"
#endif
#ifndef included_MPQC_ComponentFactory_hxx
#include "MPQC_ComponentFactory.hxx"
#endif
#ifndef included_ccaffeine_ports_ComponentFactory_hxx
#include "ccaffeine_ports_ComponentFactory.hxx"
#endif
#ifndef included_gov_cca_CCAException_hxx
#include "gov_cca_CCAException.hxx"
#endif
#ifndef included_gov_cca_Component_hxx
#include "gov_cca_Component.hxx"
#endif
#ifndef included_gov_cca_ComponentClassDescription_hxx
#include "gov_cca_ComponentClassDescription.hxx"
#endif
#ifndef included_sidl_BaseClass_hxx
#include "sidl_BaseClass.hxx"
#endif
#ifndef included_sidl_BaseInterface_hxx
#include "sidl_BaseInterface.hxx"
#endif
#ifndef included_sidl_ClassInfo_hxx
#include "sidl_ClassInfo.hxx"
#endif
#ifndef included_sidl_RuntimeException_hxx
#include "sidl_RuntimeException.hxx"
#endif


// DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory._includes)
// Put additional includes or other arbitrary code here...
// DO-NOT-DELETE splicer.end(MPQC.ComponentFactory._includes)

namespace MPQC { 

  /**
   * Symbol "MPQC.ComponentFactory" (version 0.2)
   */
  class ComponentFactory_impl : public virtual ::MPQC::ComponentFactory 
  // DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory._inherits)

  /** ComponentFactory implements a CCA standard component
      interface for component factories.  This class is used to
      inform the embedded framework of available components in a
      statically linked executable.

      This is an implementation of a SIDL interface.
      The stub code is generated by the Babel tool.  Do not make
      modifications outside of splicer blocks, as these will be lost.
      This is a server implementation for a Babel class, the Babel
      client code is provided by the cca-spec-babel package.
   */

  // DO-NOT-DELETE splicer.end(MPQC.ComponentFactory._inherits)
  {

  // All data marked protected will be accessable by 
  // descendant Impl classes
  protected:

    bool _wrapped;

    // DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory._implementation)
    std::vector< gov::cca::ComponentClassDescription > descriptions;
    // DO-NOT-DELETE splicer.end(MPQC.ComponentFactory._implementation)

  public:
    // default constructor, used for data wrapping(required)
    ComponentFactory_impl();
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    ComponentFactory_impl( struct MPQC_ComponentFactory__object * s ) : 
      StubBase(s,true), _wrapped(false) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~ComponentFactory_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // true if this object was created by a user newing the impl
    inline bool _isWrapped() {return _wrapped;}

    // static class initializer
    static void _load();

  public:

    /**
     * user defined non-static method.
     */
    void
    addDescription_impl (
      /* in */const ::std::string& className,
      /* in */const ::std::string& classAlias
    )
    ;


    /**
     *  
     * Collect the currently obtainable class name strings from
     * factories known to the builder and the from the
     * already instantiated components.
     * @return The list of class description, which may be empty, that are
     * known a priori to contain valid values for the className
     * argument of createInstance. 
     * @throws CCAException in the event of error.
     */
    ::sidl::array< ::gov::cca::ComponentClassDescription>
    getAvailableComponentClasses_impl() // throws:
    //     ::gov::cca::CCAException
    //     ::sidl::RuntimeException
    ;

    /**
     *  the component instance returned is nil if the name is unknown
     * to the factory. The component is raw: it has been constructed
     * but not initialized via setServices.
     */
    ::gov::cca::Component
    createComponentInstance_impl (
      /* in */const ::std::string& className
    )
    ;


    /**
     *  reclaim any resources the factory may have associated with
     * the port it is using. This will occur after the
     * normal component shutdown  (ala componentrelease) is finished. 
     */
    void
    destroyComponentInstance_impl (
      /* in */const ::std::string& className,
      /* in */::gov::cca::Component c
    )
    ;

  };  // end class ComponentFactory_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory._misc)
// Put miscellaneous things here...
// DO-NOT-DELETE splicer.end(MPQC.ComponentFactory._misc)

#endif
