// 
// File:          MPQC_Chemistry_QC_GaussianShell_Impl.hh
// Symbol:        MPQC.Chemistry_QC_GaussianShell-v0.2
// Symbol Type:   class
// Babel Version: 0.9.8
// Description:   Server-side implementation for MPQC.Chemistry_QC_GaussianShell
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.9.8
// 

#ifndef included_MPQC_Chemistry_QC_GaussianShell_Impl_hh
#define included_MPQC_Chemistry_QC_GaussianShell_Impl_hh

#ifndef included_sidl_cxx_hh
#include "sidl_cxx.hh"
#endif
#ifndef included_MPQC_Chemistry_QC_GaussianShell_IOR_h
#include "MPQC_Chemistry_QC_GaussianShell_IOR.h"
#endif
// 
// Includes for all method dependencies.
// 
#ifndef included_Chemistry_QC_GaussianBasis_AngularType_hh
#include "Chemistry_QC_GaussianBasis_AngularType.hh"
#endif
#ifndef included_MPQC_Chemistry_QC_GaussianShell_hh
#include "MPQC_Chemistry_QC_GaussianShell.hh"
#endif
#ifndef included_sidl_BaseInterface_hh
#include "sidl_BaseInterface.hh"
#endif
#ifndef included_sidl_ClassInfo_hh
#include "sidl_ClassInfo.hh"
#endif


// DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_GaussianShell._includes)
// Put additional includes or other arbitrary code here...
// DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_GaussianShell._includes)

namespace MPQC { 

  /**
   * Symbol "MPQC.Chemistry_QC_GaussianShell" (version 0.2)
   */
  class Chemistry_QC_GaussianShell_impl
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_GaussianShell._inherits)
  // Put additional inheritance here...
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_GaussianShell._inherits)
  {

  private:
    // Pointer back to IOR.
    // Use this to dispatch back through IOR vtable.
    Chemistry_QC_GaussianShell self;

    // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_GaussianShell._implementation)
    // Put additional implementation details here...
    // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_GaussianShell._implementation)

  private:
    // private default constructor (required)
    Chemistry_QC_GaussianShell_impl() {} 

  public:
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    Chemistry_QC_GaussianShell_impl( struct 
      MPQC_Chemistry_QC_GaussianShell__object * s ) : self(s,true) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~Chemistry_QC_GaussianShell_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

  public:

    /**
     * user defined non-static method.
     */
    int64_t
    get_n_contraction() throw () 
    ;
    /**
     * user defined non-static method.
     */
    int64_t
    get_n_primitive() throw () 
    ;
    /**
     * user defined non-static method.
     */
    double
    get_contraction_coef (
      /*in*/ int64_t connum,
      /*in*/ int64_t expnum
    )
    throw () 
    ;

    /**
     * user defined non-static method.
     */
    double
    get_exponent (
      /*in*/ int64_t expnum
    )
    throw () 
    ;

    /**
     * user defined non-static method.
     */
    int64_t
    get_angular_momentum (
      /*in*/ int64_t connum
    )
    throw () 
    ;

    /**
     * user defined non-static method.
     */
    ::Chemistry::QC::GaussianBasis::AngularType
    get_angular_type() throw () 
    ;
  };  // end class Chemistry_QC_GaussianShell_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_GaussianShell._misc)
// Put miscellaneous things here...
// DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_GaussianShell._misc)

#endif
