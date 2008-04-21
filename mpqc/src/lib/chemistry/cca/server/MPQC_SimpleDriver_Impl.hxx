// 
// File:          MPQC_SimpleDriver_Impl.hxx
// Symbol:        MPQC.SimpleDriver-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.SimpleDriver
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 

#ifndef included_MPQC_SimpleDriver_Impl_hxx
#define included_MPQC_SimpleDriver_Impl_hxx

#ifndef included_sidl_cxx_hxx
#include "sidl_cxx.hxx"
#endif
#ifndef included_MPQC_SimpleDriver_IOR_h
#include "MPQC_SimpleDriver_IOR.h"
#endif
#ifndef included_MPQC_SimpleDriver_hxx
#include "MPQC_SimpleDriver.hxx"
#endif
#ifndef included_gov_cca_CCAException_hxx
#include "gov_cca_CCAException.hxx"
#endif
#ifndef included_gov_cca_Component_hxx
#include "gov_cca_Component.hxx"
#endif
#ifndef included_gov_cca_Services_hxx
#include "gov_cca_Services.hxx"
#endif
#ifndef included_gov_cca_ports_GoPort_hxx
#include "gov_cca_ports_GoPort.hxx"
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


// DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver._includes)
#include "Chemistry_QC_ModelFactoryInterface.hxx"
#include "Chemistry_QC_ModelInterface.hxx"
#include "Chemistry_MoleculeInterface.hxx"
#include "dc/babel.new/babel-cca/server/ccaffeine_TypeMap.hxx"
#include "parameters/parametersStar.h"
#include "gov_cca_ports_ParameterPortFactory.hxx"
#include "gov_cca_ports_ParameterPort.hxx"
// DO-NOT-DELETE splicer.end(MPQC.SimpleDriver._includes)

namespace MPQC { 

  /**
   * Symbol "MPQC.SimpleDriver" (version 0.2)
   */
  class SimpleDriver_impl : public virtual ::MPQC::SimpleDriver 
  // DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver._inherits)
  // Insert-Code-Here {MPQC.SimpleDriver._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(MPQC.SimpleDriver._inherits)
  {

  // All data marked protected will be accessable by 
  // descendant Impl classes
  protected:

    bool _wrapped;

    // DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver._implementation)

    gov::cca::Services services_;
    Chemistry::QC::ModelFactoryInterface factory_;
    Chemistry::QC::ModelInterface model_;
    Chemistry::MoleculeInterface molecule_;
    BoolParameter *grad_param_;
    bool do_grad_;

    gov::cca::TypeMap tm_;
    gov::cca::ports::ParameterPortFactory ppf_;
    gov::cca::ports::ParameterPort pp_;

    // DO-NOT-DELETE splicer.end(MPQC.SimpleDriver._implementation)

  public:
    // default constructor, used for data wrapping(required)
    SimpleDriver_impl();
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    SimpleDriver_impl( struct MPQC_SimpleDriver__object * s ) : StubBase(s,
      true), _wrapped(false) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~SimpleDriver_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // true if this object was created by a user newing the impl
    inline bool _isWrapped() {return _wrapped;}

    // static class initializer
    static void _load();

  public:


    /**
     *  Starts up a component presence in the calling framework.
     * @param services the component instance's handle on the framework world.
     * Contracts concerning services and setServices:
     * 
     * The component interaction with the CCA framework
     * and Ports begins on the call to setServices by the framework.
     * 
     * This function is called exactly once for each instance created
     * by the framework.
     * 
     * The argument services will never be nil/null.
     * 
     * Those uses ports which are automatically connected by the framework
     * (so-called service-ports) may be obtained via getPort during
     * setServices.
     */
    void
    setServices_impl (
      /* in */::gov::cca::Services services
    )
    // throws:
    //     ::gov::cca::CCAException
    //     ::sidl::RuntimeException
    ;


    /**
     *  
     * Execute some encapsulated functionality on the component. 
     * Return 0 if ok, -1 if internal error but component may be 
     * used further, and -2 if error so severe that component cannot
     * be further used safely.
     */
    int32_t
    go_impl() ;
  };  // end class SimpleDriver_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver._misc)
// Insert-Code-Here {MPQC.SimpleDriver._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(MPQC.SimpleDriver._misc)

#endif
