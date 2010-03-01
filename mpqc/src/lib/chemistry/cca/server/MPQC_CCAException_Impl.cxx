// 
// File:          MPQC_CCAException_Impl.cxx
// Symbol:        MPQC.CCAException-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.CCAException
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 
#include "MPQC_CCAException_Impl.hxx"

// 
// Includes for all method dependencies.
// 
#ifndef included_gov_cca_CCAExceptionType_hxx
#include "gov_cca_CCAExceptionType.hxx"
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
#ifndef included_sidl_NotImplementedException_hxx
#include "sidl_NotImplementedException.hxx"
#endif
// DO-NOT-DELETE splicer.begin(MPQC.CCAException._includes)
#include <sstream>
#include <iomanip>
// DO-NOT-DELETE splicer.end(MPQC.CCAException._includes)

// special constructor, used for data wrapping(required).  Do not put code here unless you really know what you're doing!
MPQC::CCAException_impl::CCAException_impl() : StubBase(reinterpret_cast< 
  void*>(::MPQC::CCAException::_wrapObj(reinterpret_cast< void*>(this))),false) 
  , _wrapped(true){ 
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException._ctor2)
  // Insert-Code-Here {MPQC.CCAException._ctor2} (ctor2)
  // DO-NOT-DELETE splicer.end(MPQC.CCAException._ctor2)
}

// user defined constructor
void MPQC::CCAException_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException._ctor)
  type_ = gov::cca::CCAExceptionType_Nonstandard;
  // DO-NOT-DELETE splicer.end(MPQC.CCAException._ctor)
}

// user defined destructor
void MPQC::CCAException_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException._dtor)
  // Insert-Code-Here {MPQC.CCAException._dtor} (destructor)
  // DO-NOT-DELETE splicer.end(MPQC.CCAException._dtor)
}

// static class initializer
void MPQC::CCAException_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException._load)
  // Insert-Code-Here {MPQC.CCAException._load} (class initialization)
  // DO-NOT-DELETE splicer.end(MPQC.CCAException._load)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 * Method:  setCCAExceptionType[]
 */
void
MPQC::CCAException_impl::setCCAExceptionType_impl (
  /* in */::gov::cca::CCAExceptionType type ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException.setCCAExceptionType)
  type_ = type;
  // DO-NOT-DELETE splicer.end(MPQC.CCAException.setCCAExceptionType)
}

/**
 * Method:  packObj[]
 */
void
MPQC::CCAException_impl::packObj_impl (
  /* in */::sidl::io::Serializer& ser ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException.packObj)
  // Insert-Code-Here {MPQC.CCAException.packObj} (packObj method)
  // 
  // This method has not been implemented
  // 
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "packObj");
    throw ex;
  // DO-NOT-DELETE splicer.end(MPQC.CCAException.packObj)
}

/**
 * Method:  unpackObj[]
 */
void
MPQC::CCAException_impl::unpackObj_impl (
  /* in */::sidl::io::Deserializer& des ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException.unpackObj)
  // Insert-Code-Here {MPQC.CCAException.unpackObj} (unpackObj method)
  // 
  // This method has not been implemented
  // 
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "unpackObj");
    throw ex;
  // DO-NOT-DELETE splicer.end(MPQC.CCAException.unpackObj)
}

/**
 * Return the message associated with the exception.
 */
::std::string
MPQC::CCAException_impl::getNote_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException.getNote)
  return note_;
  // DO-NOT-DELETE splicer.end(MPQC.CCAException.getNote)
}

/**
 * Set the message associated with the exception.
 */
void
MPQC::CCAException_impl::setNote_impl (
  /* in */const ::std::string& message ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException.setNote)
  note_ = message;
  // DO-NOT-DELETE splicer.end(MPQC.CCAException.setNote)
}

/**
 * Returns formatted string containing the concatenation of all 
 * tracelines.
 */
::std::string
MPQC::CCAException_impl::getTrace_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException.getTrace)
  std::ostringstream oss;
  oss << "MPQC::CCAException: trace:" << std::endl;
  for (int i=0; i<trace_.size(); i++) {
      oss << "  " << std::setw(2) << i << trace_[i] << std::endl;
    }
  return oss.str();
  // DO-NOT-DELETE splicer.end(MPQC.CCAException.getTrace)
}

/**
 * Adds a stringified entry/line to the stack trace.
 */
void
MPQC::CCAException_impl::add_impl (
  /* in */const ::std::string& traceline ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException.addLine)
  trace_.push_back(traceline);
  // DO-NOT-DELETE splicer.end(MPQC.CCAException.addLine)
}

/**
 * Formats and adds an entry to the stack trace based on the 
 * file name, line number, and method name.
 */
void
MPQC::CCAException_impl::add_impl (
  /* in */const ::std::string& filename,
  /* in */int32_t lineno,
  /* in */const ::std::string& methodname ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException.add)
  std::ostringstream oss;
  oss << filename << ":" << lineno << " " << methodname;
  add(oss.str());
  // DO-NOT-DELETE splicer.end(MPQC.CCAException.add)
}

/**
 * Method:  getCCAExceptionType[]
 */
::gov::cca::CCAExceptionType
MPQC::CCAException_impl::getCCAExceptionType_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException.getCCAExceptionType)
  return type_;
  // DO-NOT-DELETE splicer.end(MPQC.CCAException.getCCAExceptionType)
}


// DO-NOT-DELETE splicer.begin(MPQC.CCAException._misc)
// Insert-Code-Here {MPQC.CCAException._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(MPQC.CCAException._misc)

