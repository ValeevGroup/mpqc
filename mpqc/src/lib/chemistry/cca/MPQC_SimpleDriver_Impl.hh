// 
// File:          MPQC_SimpleDriver_Impl.hh
// Symbol:        MPQC.SimpleDriver-v0.2
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for MPQC.SimpleDriver
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// xml-url       = /home/jpkenny/src/mpqc-libint2.build-shared/src/lib/chemistry/cca/../../../../lib/cca/repo/MPQC.SimpleDriver-v0.2.xml
// 

#ifndef included_MPQC_SimpleDriver_Impl_hh
#define included_MPQC_SimpleDriver_Impl_hh

#ifndef included_sidl_cxx_hh
#include "sidl_cxx.hh"
#endif
#ifndef included_MPQC_SimpleDriver_IOR_h
#include "MPQC_SimpleDriver_IOR.h"
#endif
// 
// Includes for all method dependencies.
// 
#ifndef included_MPQC_SimpleDriver_hh
#include "MPQC_SimpleDriver.hh"
#endif
#ifndef included_gov_cca_CCAException_hh
#include "gov_cca_CCAException.hh"
#endif
#ifndef included_gov_cca_Services_hh
#include "gov_cca_Services.hh"
#endif
#ifndef included_sidl_BaseInterface_hh
#include "sidl_BaseInterface.hh"
#endif
#ifndef included_sidl_ClassInfo_hh
#include "sidl_ClassInfo.hh"
#endif


// DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver._includes)
#include "Chemistry_QC_ModelFactory.hh"
#include "Chemistry_QC_Model.hh"
#include "Chemistry_Molecule.hh"
#include "dc/babel/babel-cca/server/ccaffeine_TypeMap.hh"
#include "dc/babel/babel-cca/server/ccaffeine_ports_PortTranslator.hh"
#include "cca.h"
#include "util/IO.h"
#include "jc++/jc++.h"
#include "jc++/util/jc++util.h"
#include "parameters/parametersStar.h"
#include "port/portInterfaces.h"
#include "port/supportInterfaces.h"
// DO-NOT-DELETE splicer.end(MPQC.SimpleDriver._includes)

namespace MPQC { 

  /**
   * Symbol "MPQC.SimpleDriver" (version 0.2)
   */
  class SimpleDriver_impl
  // DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver._inherits)
  // Insert-Code-Here {MPQC.SimpleDriver._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(MPQC.SimpleDriver._inherits)
  {

  private:
    // Pointer back to IOR.
    // Use this to dispatch back through IOR vtable.
    SimpleDriver self;

    // DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver._implementation)

    gov::cca::Services services_;
    Chemistry::QC::ModelFactory factory_;
    Chemistry::QC::Model model_;
    Chemistry::Molecule molecule_;
    BoolParameter *grad_param_;
    bool do_grad_;

    ConfigurableParameterPort*
    setup_parameters(ConfigurableParameterFactory *);

    // DO-NOT-DELETE splicer.end(MPQC.SimpleDriver._implementation)

  private:
    // private default constructor (required)
    SimpleDriver_impl() 
    {} 

  public:
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    SimpleDriver_impl( struct MPQC_SimpleDriver__object * s ) : self(s,
      true) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~SimpleDriver_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // static class initializer
    static void _load();

  public:


    /**
     * Starts up a component presence in the calling framework.
     * @param services the component instance's handle on the framework world.
     * Contracts concerning Svc and setServices:
     * 
     * The component interaction with the CCA framework
     * and Ports begins on the call to setServices by the framework.
     * 
     * This function is called exactly once for each instance created
     * by the framework.
     * 
     * The argument Svc will never be nil/null.
     * 
     * Those uses ports which are automatically connected by the framework
     * (so-called service-ports) may be obtained via getPort during
     * setServices.
     */
    void
    setServices (
      /* in */ ::gov::cca::Services services
    )
    throw ( 
      ::gov::cca::CCAException
    );


    /**
     * Execute some encapsulated functionality on the component. 
     * Return 0 if ok, -1 if internal error but component may be 
     * used further, and -2 if error so severe that component cannot
     * be further used safely.
     */
    int32_t
    go() throw () 
    ;
  };  // end class SimpleDriver_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver._misc)
// Insert-Code-Here {MPQC.SimpleDriver._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(MPQC.SimpleDriver._misc)

#endif
