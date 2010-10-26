// 
// File:          MPQC_CCAException_Impl.hxx
// Symbol:        MPQC.CCAException-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.CCAException
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 

#ifndef included_MPQC_CCAException_Impl_hxx
#define included_MPQC_CCAException_Impl_hxx

#ifndef included_sidl_cxx_hxx
#include "sidl_cxx.hxx"
#endif
#ifndef included_MPQC_CCAException_IOR_h
#include "MPQC_CCAException_IOR.h"
#endif
#ifndef included_MPQC_CCAException_hxx
#include "MPQC_CCAException.hxx"
#endif
#ifndef included_gov_cca_CCAException_hxx
#include "gov_cca_CCAException.hxx"
#endif
#ifndef included_gov_cca_CCAExceptionType_hxx
#include "gov_cca_CCAExceptionType.hxx"
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
#ifndef included_sidl_io_Deserializer_hxx
#include "sidl_io_Deserializer.hxx"
#endif
#ifndef included_sidl_io_Serializer_hxx
#include "sidl_io_Serializer.hxx"
#endif


// DO-NOT-DELETE splicer.begin(MPQC.CCAException._hincludes)
#include <vector>
#include <string>
// DO-NOT-DELETE splicer.end(MPQC.CCAException._hincludes)

namespace MPQC { 

  /**
   * Symbol "MPQC.CCAException" (version 0.2)
   */
  class CCAException_impl : public virtual ::MPQC::CCAException 
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException._inherits)
  // Insert-Code-Here {MPQC.CCAException._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(MPQC.CCAException._inherits)

  {

  // All data marked protected will be accessable by 
  // descendant Impl classes
  protected:

    bool _wrapped;

    // DO-NOT-DELETE splicer.begin(MPQC.CCAException._implementation)
      std::string note_;
      std::vector<std::string> trace_;
      gov::cca::CCAExceptionType type_;
    // DO-NOT-DELETE splicer.end(MPQC.CCAException._implementation)

  public:
    // default constructor, used for data wrapping(required)
    CCAException_impl();
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
      CCAException_impl( struct MPQC_CCAException__object * ior ) : StubBase(
        ior,true), 
      ::sidl::io::Serializable((ior==NULL) ? NULL : &((
        *ior).d_sidl_io_serializable)),
      ::sidl::BaseException((ior==NULL) ? NULL : &((
        *ior).d_sidl_baseexception)),
    ::gov::cca::CCAException((ior==NULL) ? NULL : &((
      *ior).d_gov_cca_ccaexception)) , _wrapped(false) {
      ior->d_data = this;
      _ctor();
    }


    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~CCAException_impl() { _dtor(); }

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
    setCCAExceptionType_impl (
      /* in */::gov::cca::CCAExceptionType type
    )
    ;

    /**
     * user defined non-static method.
     */
    void
    packObj_impl (
      /* in */::sidl::io::Serializer& ser
    )
    ;

    /**
     * user defined non-static method.
     */
    void
    unpackObj_impl (
      /* in */::sidl::io::Deserializer& des
    )
    ;


    /**
     * Return the message associated with the exception.
     */
    ::std::string
    getNote_impl() ;

    /**
     * Set the message associated with the exception.
     */
    void
    setNote_impl (
      /* in */const ::std::string& message
    )
    ;


    /**
     * Returns formatted string containing the concatenation of all 
     * tracelines.
     */
    ::std::string
    getTrace_impl() ;

    /**
     * Adds a stringified entry/line to the stack trace.
     */
    void
    add_impl (
      /* in */const ::std::string& traceline
    )
    ;


    /**
     * Formats and adds an entry to the stack trace based on the 
     * file name, line number, and method name.
     */
    void
    add_impl (
      /* in */const ::std::string& filename,
      /* in */int32_t lineno,
      /* in */const ::std::string& methodname
    )
    ;

    /**
     * user defined non-static method.
     */
    ::gov::cca::CCAExceptionType
    getCCAExceptionType_impl() ;
  };  // end class CCAException_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.CCAException._hmisc)
// Insert-Code-Here {MPQC.CCAException._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(MPQC.CCAException._hmisc)

#endif
