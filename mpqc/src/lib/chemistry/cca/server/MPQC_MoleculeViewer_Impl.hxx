// 
// File:          MPQC_MoleculeViewer_Impl.hxx
// Symbol:        MPQC.MoleculeViewer-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.MoleculeViewer
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 

#ifndef included_MPQC_MoleculeViewer_Impl_hxx
#define included_MPQC_MoleculeViewer_Impl_hxx

#ifndef included_sidl_cxx_hxx
#include "sidl_cxx.hxx"
#endif
#ifndef included_MPQC_MoleculeViewer_IOR_h
#include "MPQC_MoleculeViewer_IOR.h"
#endif
#ifndef included_Chemistry_MoleculeInterface_hxx
#include "Chemistry_MoleculeInterface.hxx"
#endif
#ifndef included_Chemistry_MoleculeViewerInterface_hxx
#include "Chemistry_MoleculeViewerInterface.hxx"
#endif
#ifndef included_MPQC_MoleculeViewer_hxx
#include "MPQC_MoleculeViewer.hxx"
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


// DO-NOT-DELETE splicer.begin(MPQC.MoleculeViewer._includes)

#define USE_SOCKET 1
#if USE_SOCKET
#include "socket.h"
#endif // USE_SOCKET

// DO-NOT-DELETE splicer.end(MPQC.MoleculeViewer._includes)

namespace MPQC { 

  /**
   * Symbol "MPQC.MoleculeViewer" (version 0.2)
   */
  class MoleculeViewer_impl : public virtual ::MPQC::MoleculeViewer 
  // DO-NOT-DELETE splicer.begin(MPQC.MoleculeViewer._inherits)
  // Insert-Code-Here {MPQC.MoleculeViewer._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(MPQC.MoleculeViewer._inherits)
  {

  // All data marked protected will be accessable by 
  // descendant Impl classes
  protected:

    bool _wrapped;

    // DO-NOT-DELETE splicer.begin(MPQC.MoleculeViewer._implementation)

      gov::cca::Services services_;
      Chemistry::MoleculeInterface molecule_;
      int is_updated;
#if USE_SOCKET
      TCPClientConnection socket_;
#endif // USE_SOCKET

    // DO-NOT-DELETE splicer.end(MPQC.MoleculeViewer._implementation)

  public:
    // default constructor, used for data wrapping(required)
    MoleculeViewer_impl();
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    MoleculeViewer_impl( struct MPQC_MoleculeViewer__object * s ) : StubBase(s,
      true), _wrapped(false) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~MoleculeViewer_impl() { _dtor(); }

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
    set_molecule_impl (
      /* in */::Chemistry::MoleculeInterface molecule
    )
    ;

    /**
     * user defined non-static method.
     */
    void
    set_coor_impl (
      /* in */const ::std::string& coords
    )
    ;

    /**
     * user defined non-static method.
     */
    void
    run_gui_impl() ;
    /**
     * user defined non-static method.
     */
    void
    draw_impl() ;

    /**
     *  This should be called when the object is no longer needed.
     * No other members may be called after finalize. 
     */
    int32_t
    finalize_impl() ;

    /**
     *  Starts up a component presence in the calling framework.
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
    setServices_impl (
      /* in */::gov::cca::Services services
    )
    // throws:
    //     ::gov::cca::CCAException
    //     ::sidl::RuntimeException
    ;

  };  // end class MoleculeViewer_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.MoleculeViewer._misc)
// Insert-Code-Here {MPQC.MoleculeViewer._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(MPQC.MoleculeViewer._misc)

#endif
