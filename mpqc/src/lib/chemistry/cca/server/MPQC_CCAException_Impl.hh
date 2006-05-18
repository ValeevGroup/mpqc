// 
// File:          MPQC_CCAException_Impl.hh
// Symbol:        MPQC.CCAException-v0.2
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for MPQC.CCAException
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// 

#ifndef included_MPQC_CCAException_Impl_hh
#define included_MPQC_CCAException_Impl_hh

#ifndef included_sidl_cxx_hh
#include "sidl_cxx.hh"
#endif
#ifndef included_MPQC_CCAException_IOR_h
#include "MPQC_CCAException_IOR.h"
#endif
// 
// Includes for all method dependencies.
// 
#ifndef included_MPQC_CCAException_hh
#include "MPQC_CCAException.hh"
#endif
#ifndef included_gov_cca_CCAExceptionType_hh
#include "gov_cca_CCAExceptionType.hh"
#endif
#ifndef included_sidl_BaseInterface_hh
#include "sidl_BaseInterface.hh"
#endif
#ifndef included_sidl_ClassInfo_hh
#include "sidl_ClassInfo.hh"
#endif


// DO-NOT-DELETE splicer.begin(MPQC.CCAException._includes)
#include <vector>
#include <string>
// DO-NOT-DELETE splicer.end(MPQC.CCAException._includes)

namespace MPQC { 

  /**
   * Symbol "MPQC.CCAException" (version 0.2)
   */
  class CCAException_impl
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException._inherits)
  // Insert-Code-Here {MPQC.CCAException._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(MPQC.CCAException._inherits)
  {

  private:
    // Pointer back to IOR.
    // Use this to dispatch back through IOR vtable.
    CCAException self;

    // DO-NOT-DELETE splicer.begin(MPQC.CCAException._implementation)
      std::string note_;
      std::vector<std::string> trace_;
      gov::cca::CCAExceptionType type_;
    // DO-NOT-DELETE splicer.end(MPQC.CCAException._implementation)

  private:
    // private default constructor (required)
    CCAException_impl() 
    {} 

  public:
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    CCAException_impl( struct MPQC_CCAException__object * s ) : self(s,
      true) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~CCAException_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // static class initializer
    static void _load();

  public:

    /**
     * user defined non-static method.
     */
    void
    setCCAExceptionType (
      /* in */ ::gov::cca::CCAExceptionType type
    )
    throw () 
    ;


    /**
     * Return the message associated with the exception.
     */
    ::std::string
    getNote() throw () 
    ;

    /**
     * Set the message associated with the exception.
     */
    void
    setNote (
      /* in */ const ::std::string& message
    )
    throw () 
    ;


    /**
     * Returns formatted string containing the concatenation of all 
     * tracelines.
     */
    ::std::string
    getTrace() throw () 
    ;

    /**
     * Adds a stringified entry/line to the stack trace.
     */
    void
    add (
      /* in */ const ::std::string& traceline
    )
    throw () 
    ;


    /**
     * Formats and adds an entry to the stack trace based on the 
     * file name, line number, and method name.
     */
    void
    add (
      /* in */ const ::std::string& filename,
      /* in */ int32_t lineno,
      /* in */ const ::std::string& methodname
    )
    throw () 
    ;

    /**
     * user defined non-static method.
     */
    ::gov::cca::CCAExceptionType
    getCCAExceptionType() throw () 
    ;
  };  // end class CCAException_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.CCAException._misc)
// Insert-Code-Here {MPQC.CCAException._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(MPQC.CCAException._misc)

#endif
