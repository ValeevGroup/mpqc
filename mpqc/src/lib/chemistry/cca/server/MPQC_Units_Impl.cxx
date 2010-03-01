// 
// File:          MPQC_Units_Impl.cxx
// Symbol:        MPQC.Units-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.Units
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 
#include "MPQC_Units_Impl.hxx"

// 
// Includes for all method dependencies.
// 
#ifndef included_sidl_BaseInterface_hxx
#include "sidl_BaseInterface.hxx"
#endif
#ifndef included_sidl_ClassInfo_hxx
#include "sidl_ClassInfo.hxx"
#endif
#ifndef included_sidl_NotImplementedException_hxx
#include "sidl_NotImplementedException.hxx"
#endif
// DO-NOT-DELETE splicer.begin(MPQC.Units._includes)
// Insert-Code-Here {MPQC.Units._includes} (additional includes or code)
// DO-NOT-DELETE splicer.end(MPQC.Units._includes)

// special constructor, used for data wrapping(required).  Do not put code here unless you really know what you're doing!
MPQC::Units_impl::Units_impl() : StubBase(reinterpret_cast< void*>(
  ::MPQC::Units::_wrapObj(reinterpret_cast< void*>(this))),false) , _wrapped(
  true){ 
  // DO-NOT-DELETE splicer.begin(MPQC.Units._ctor2)
  // Insert-Code-Here {MPQC.Units._ctor2} (ctor2)
  // DO-NOT-DELETE splicer.end(MPQC.Units._ctor2)
}

// user defined constructor
void MPQC::Units_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.Units._ctor)
  // Insert-Code-Here {MPQC.Units._ctor} (constructor)
  // DO-NOT-DELETE splicer.end(MPQC.Units._ctor)
}

// user defined destructor
void MPQC::Units_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.Units._dtor)
  // Insert-Code-Here {MPQC.Units._dtor} (destructor)
  // DO-NOT-DELETE splicer.end(MPQC.Units._dtor)
}

// static class initializer
void MPQC::Units_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.Units._load)
  // Insert-Code-Here {MPQC.Units._load} (class initialization)
  // DO-NOT-DELETE splicer.end(MPQC.Units._load)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 *  Initializes the units as a human readable string
 * options are "angstroms" or "bohr" 
 */
void
MPQC::Units_impl::initialize_impl (
  /* in */const ::std::string& unitname ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Units.initialize)
  // Insert-Code-Here {MPQC.Units.initialize} (initialize method)
  // 
  // This method has not been implemented
  // 
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "initialize");
    throw ex;
  // DO-NOT-DELETE splicer.end(MPQC.Units.initialize)
}

/**
 *  Returns the units as a human readable string. 
 */
::std::string
MPQC::Units_impl::get_unit_name_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Units.get_unit_name)

  return units->string_rep();

  // DO-NOT-DELETE splicer.end(MPQC.Units.get_unit_name)
}

/**
 *  Converts from self's units to the given unit name. 
 */
double
MPQC::Units_impl::convert_to_impl (
  /* in */const ::std::string& unitname ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Units.convert_to)

  if (units.null()) return 0;
  sc::Ref<sc::Units> u = new sc::Units(unitname.c_str());
  return units->to(u);

  // DO-NOT-DELETE splicer.end(MPQC.Units.convert_to)
}

/**
 *  Converts to self's units from the given unit name. 
 */
double
MPQC::Units_impl::convert_from_impl (
  /* in */const ::std::string& unitname ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Units.convert_from)

  if (units.null()) return 0;
  sc::Ref<sc::Units> u = new sc::Units(unitname.c_str());
  return units->from(u);

  // DO-NOT-DELETE splicer.end(MPQC.Units.convert_from)
}


// DO-NOT-DELETE splicer.begin(MPQC.Units._misc)

void
MPQC::Units_impl::set_units(const sc::Ref<sc::Units> &u)
{
  units = u;
}

// DO-NOT-DELETE splicer.end(MPQC.Units._misc)

