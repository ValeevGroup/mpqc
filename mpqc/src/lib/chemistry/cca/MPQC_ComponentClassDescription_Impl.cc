// 
// File:          MPQC_ComponentClassDescription_Impl.cc
// Symbol:        MPQC.ComponentClassDescription-v0.2
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for MPQC.ComponentClassDescription
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// xml-url       = /home/jpkenny/src/mpqc-libint2.build-shared/src/lib/chemistry/cca/../../../../lib/cca/repo/MPQC.ComponentClassDescription-v0.2.xml
// 
#include "MPQC_ComponentClassDescription_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.ComponentClassDescription._includes)
// Put additional includes or other arbitrary code here...
// DO-NOT-DELETE splicer.end(MPQC.ComponentClassDescription._includes)

// user-defined constructor.
void MPQC::ComponentClassDescription_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.ComponentClassDescription._ctor)
  // add construction details here
  // DO-NOT-DELETE splicer.end(MPQC.ComponentClassDescription._ctor)
}

// user-defined destructor.
void MPQC::ComponentClassDescription_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.ComponentClassDescription._dtor)
  // add destruction details here
  // DO-NOT-DELETE splicer.end(MPQC.ComponentClassDescription._dtor)
}

// static class initializer.
void MPQC::ComponentClassDescription_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.ComponentClassDescription._load)
  // guaranteed to be called at most once before any other method in this class
  // DO-NOT-DELETE splicer.end(MPQC.ComponentClassDescription._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Method:  initialize[]
 */
void
MPQC::ComponentClassDescription_impl::initialize (
  /* in */ const ::std::string& className,
  /* in */ const ::std::string& classAlias ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.ComponentClassDescription.initialize)
  cName = className;
  cAlias = classAlias;
  // DO-NOT-DELETE splicer.end(MPQC.ComponentClassDescription.initialize)
}

/**
 *  Returns the class name provided in 
 *   <code>BuilderService.createInstance()</code>
 *   or in
 *   <code>AbstractFramework.getServices()</code>.
 *  <p>
 *  Throws <code>CCAException</code> if <code>ComponentClassDescription</code> is invalid.
 */
::std::string
MPQC::ComponentClassDescription_impl::getComponentClassName ()
throw ( 
  ::gov::cca::CCAException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.ComponentClassDescription.getComponentClassName)
  return cName;
  // DO-NOT-DELETE splicer.end(MPQC.ComponentClassDescription.getComponentClassName)
}


// DO-NOT-DELETE splicer.begin(MPQC.ComponentClassDescription._misc)
// Put miscellaneous code here
// DO-NOT-DELETE splicer.end(MPQC.ComponentClassDescription._misc)

