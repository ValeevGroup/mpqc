// 
// File:          MPQC_CCAException_Impl.cc
// Symbol:        MPQC.CCAException-v0.2
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for MPQC.CCAException
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// 
#include "MPQC_CCAException_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.CCAException._includes)
#include <sstream>
#include <iomanip>
// DO-NOT-DELETE splicer.end(MPQC.CCAException._includes)

// user-defined constructor.
void MPQC::CCAException_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException._ctor)
  type_ = gov::cca::CCAExceptionType_Nonstandard;
  // DO-NOT-DELETE splicer.end(MPQC.CCAException._ctor)
}

// user-defined destructor.
void MPQC::CCAException_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException._dtor)
  // Insert-Code-Here {MPQC.CCAException._dtor} (destructor)
  // DO-NOT-DELETE splicer.end(MPQC.CCAException._dtor)
}

// static class initializer.
void MPQC::CCAException_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException._load)
  // Insert-Code-Here {MPQC.CCAException._load} (class initialization)
  // DO-NOT-DELETE splicer.end(MPQC.CCAException._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Method:  setCCAExceptionType[]
 */
void
MPQC::CCAException_impl::setCCAExceptionType (
  /* in */ ::gov::cca::CCAExceptionType type ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException.setCCAExceptionType)
  type_ = type;
  // DO-NOT-DELETE splicer.end(MPQC.CCAException.setCCAExceptionType)
}

/**
 * Return the message associated with the exception.
 */
::std::string
MPQC::CCAException_impl::getNote ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException.getNote)
  return note_;
  // DO-NOT-DELETE splicer.end(MPQC.CCAException.getNote)
}

/**
 * Set the message associated with the exception.
 */
void
MPQC::CCAException_impl::setNote (
  /* in */ const ::std::string& message ) 
throw () 
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
MPQC::CCAException_impl::getTrace ()
throw () 

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
MPQC::CCAException_impl::add (
  /* in */ const ::std::string& traceline ) 
throw () 
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
MPQC::CCAException_impl::add (
  /* in */ const ::std::string& filename,
  /* in */ int32_t lineno,
  /* in */ const ::std::string& methodname ) 
throw () 
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
MPQC::CCAException_impl::getCCAExceptionType ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.CCAException.getCCAExceptionType)
  return type_;
  // DO-NOT-DELETE splicer.end(MPQC.CCAException.getCCAExceptionType)
}


// DO-NOT-DELETE splicer.begin(MPQC.CCAException._misc)
// Insert-Code-Here {MPQC.CCAException._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(MPQC.CCAException._misc)

