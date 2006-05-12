// 
// File:          MPQC_Physics_Units_Impl.cc
// Symbol:        MPQC.Physics_Units-v0.2
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for MPQC.Physics_Units
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// xml-url       = /home/jpkenny/src/mpqc-libint2.build-shared/src/lib/chemistry/cca/server/../../../../../lib/cca/repo/MPQC.Physics_Units-v0.2.xml
// 
#include "MPQC_Physics_Units_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.Physics_Units._includes)
// Put additional includes or other arbitrary code here...
// DO-NOT-DELETE splicer.end(MPQC.Physics_Units._includes)

// user-defined constructor.
void MPQC::Physics_Units_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.Physics_Units._ctor)
  // add construction details here
  // DO-NOT-DELETE splicer.end(MPQC.Physics_Units._ctor)
}

// user-defined destructor.
void MPQC::Physics_Units_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.Physics_Units._dtor)
  // add destruction details here
  // DO-NOT-DELETE splicer.end(MPQC.Physics_Units._dtor)
}

// static class initializer.
void MPQC::Physics_Units_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.Physics_Units._load)
  // guaranteed to be called at most once before any other method in this class
  // DO-NOT-DELETE splicer.end(MPQC.Physics_Units._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Initializes the units as a human readable string
 * options are "angstroms" or "bohr" 
 */
void
MPQC::Physics_Units_impl::initialize (
  /* in */ const ::std::string& unitname ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Physics_Units.initialize)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.Physics_Units.initialize)
}

/**
 * Returns the units as a human readable string. 
 */
::std::string
MPQC::Physics_Units_impl::get_unit_name ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Physics_Units.get_unit_name)
  return units->string_rep();
  // DO-NOT-DELETE splicer.end(MPQC.Physics_Units.get_unit_name)
}

/**
 * Converts from self's units to the given unit name. 
 */
double
MPQC::Physics_Units_impl::convert_to (
  /* in */ const ::std::string& unitname ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Physics_Units.convert_to)
  if (units.null()) return 0;
  sc::Ref<sc::Units> u = new sc::Units(unitname.c_str());
  return units->to(u);
  // DO-NOT-DELETE splicer.end(MPQC.Physics_Units.convert_to)
}

/**
 * Converts to self's units from the given unit name. 
 */
double
MPQC::Physics_Units_impl::convert_from (
  /* in */ const ::std::string& unitname ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Physics_Units.convert_from)
  if (units.null()) return 0;
  sc::Ref<sc::Units> u = new sc::Units(unitname.c_str());
  return units->from(u);
  // DO-NOT-DELETE splicer.end(MPQC.Physics_Units.convert_from)
}


// DO-NOT-DELETE splicer.begin(MPQC.Physics_Units._misc)
void
MPQC::Physics_Units_impl::set_units(const sc::Ref<sc::Units> &u)
{
  units = u;
}
// DO-NOT-DELETE splicer.end(MPQC.Physics_Units._misc)

