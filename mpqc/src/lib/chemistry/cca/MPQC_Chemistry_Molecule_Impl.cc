// 
// File:          MPQC_Chemistry_Molecule_Impl.cc
// Symbol:        MPQC.Chemistry_Molecule-v0.2
// Symbol Type:   class
// Babel Version: 0.8.6
// Description:   Server-side implementation for MPQC.Chemistry_Molecule
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.8.6
// 
#include "MPQC_Chemistry_Molecule_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.Chemistry_Molecule._includes)
#include "MPQC_Physics_Units_Impl.hh"
// DO-NOT-DELETE splicer.end(MPQC.Chemistry_Molecule._includes)

// user defined constructor
void MPQC::Chemistry_Molecule_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_Molecule._ctor)
  net_charge = 0;
  mol = new sc::Molecule;
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_Molecule._ctor)
}

// user defined destructor
void MPQC::Chemistry_Molecule_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_Molecule._dtor)
  // add destruction details here
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_Molecule._dtor)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 * Method:  initialize_pointer[]
 */
void
MPQC::Chemistry_Molecule_impl::initialize_pointer (
  /*in*/ void* ptr ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_Molecule.initialize_pointer)

  mol = static_cast<sc::Molecule*>(ptr);

  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_Molecule.initialize_pointer)
}

/**
 * Obtain Services handle, through which the 
 * component communicates with the framework. 
 * This is the one method that every CCA Component
 * must implement. The component will be called
 * with a nil/null Services pointer when it is
 * to shut itself down.
 */
void
MPQC::Chemistry_Molecule_impl::setServices (
  /*in*/ ::gov::cca::Services services ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_Molecule.setServices)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_Molecule.setServices)
}

/**
 * Method:  initialize[]
 */
void
MPQC::Chemistry_Molecule_impl::initialize (
  /*in*/ int32_t natom ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_Molecule.initialize)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_Molecule.initialize)
}

/**
 * Method:  get_units[]
 */
::Physics::Units
MPQC::Chemistry_Molecule_impl::get_units () 
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_Molecule.get_units)
  MPQC::Physics_Units units = MPQC::Physics_Units::_create();
  MPQC::Physics_Units_impl *unitsi
      = reinterpret_cast<MPQC::Physics_Units_impl*>(units._get_ior()->d_data);
  unitsi->set_units(new sc::Units("bohr"));
  return units;
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_Molecule.get_units)
}

/**
 * Method:  get_n_atom[]
 */
int64_t
MPQC::Chemistry_Molecule_impl::get_n_atom () 
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_Molecule.get_n_atom)
  return mol->natom();
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_Molecule.get_n_atom)
}

/**
 * Method:  get_atomic_number[]
 */
int64_t
MPQC::Chemistry_Molecule_impl::get_atomic_number (
  /*in*/ int64_t atomnum ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_Molecule.get_atomic_number)
  return mol->Z(atomnum);
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_Molecule.get_atomic_number)
}

/**
 * Method:  set_atomic_number[]
 */
void
MPQC::Chemistry_Molecule_impl::set_atomic_number (
  /*in*/ int64_t atomnum,
  /*in*/ int64_t atomic_number ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_Molecule.set_atomic_number)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_Molecule.set_atomic_number)
}

/**
 * Method:  get_net_charge[]
 */
double
MPQC::Chemistry_Molecule_impl::get_net_charge () 
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_Molecule.get_net_charge)
  return net_charge;
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_Molecule.get_net_charge)
}

/**
 * Method:  set_net_charge[]
 */
void
MPQC::Chemistry_Molecule_impl::set_net_charge (
  /*in*/ double charge ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_Molecule.set_net_charge)
  net_charge = charge;
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_Molecule.set_net_charge)
}

/**
 * Method:  get_cart_coor[]
 */
double
MPQC::Chemistry_Molecule_impl::get_cart_coor (
  /*in*/ int64_t atomnum,
  /*in*/ int32_t xyz ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_Molecule.get_cart_coor)
  return mol->r(atomnum)[xyz];
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_Molecule.get_cart_coor)
}

/**
 * Method:  set_cart_coor[]
 */
void
MPQC::Chemistry_Molecule_impl::set_cart_coor (
  /*in*/ int64_t atomnum,
  /*in*/ int32_t xyz,
  /*in*/ double val ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_Molecule.set_cart_coor)
  mol->r(atomnum)[xyz] = val;
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_Molecule.set_cart_coor)
}

/**
 * Method:  get_atomic_label[]
 */
::std::string
MPQC::Chemistry_Molecule_impl::get_atomic_label (
  /*in*/ int64_t atomnum ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_Molecule.get_atomic_label)
  return mol->label(atomnum);
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_Molecule.get_atomic_label)
}

/**
 * Method:  set_atomic_label[]
 */
void
MPQC::Chemistry_Molecule_impl::set_atomic_label (
  /*in*/ int64_t atomnum,
  /*in*/ const ::std::string& label ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_Molecule.set_atomic_label)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_Molecule.set_atomic_label)
}

/**
 * Method:  get_symmetry[]
 */
::Physics::PointGroup
MPQC::Chemistry_Molecule_impl::get_symmetry () 
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_Molecule.get_symmetry)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_Molecule.get_symmetry)
}

/**
 * Method:  get_coor[]
 */
::sidl::array<double>
MPQC::Chemistry_Molecule_impl::get_coor () 
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_Molecule.get_coor)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_Molecule.get_coor)
}

/**
 * Method:  set_coor[]
 */
void
MPQC::Chemistry_Molecule_impl::set_coor (
  /*in*/ ::sidl::array<double> x ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_Molecule.set_coor)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_Molecule.set_coor)
}


// DO-NOT-DELETE splicer.begin(MPQC.Chemistry_Molecule._misc)
// Put miscellaneous code here
// DO-NOT-DELETE splicer.end(MPQC.Chemistry_Molecule._misc)

