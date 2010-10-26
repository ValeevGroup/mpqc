// 
// File:          MPQC_GaussianBasisShell_Impl.cxx
// Symbol:        MPQC.GaussianBasisShell-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.GaussianBasisShell
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 
#include "MPQC_GaussianBasisShell_Impl.hxx"

// 
// Includes for all method dependencies.
// 
#ifndef included_Chemistry_QC_GaussianBasis_AngularType_hxx
#include "Chemistry_QC_GaussianBasis_AngularType.hxx"
#endif
#ifndef included_sidl_BaseInterface_hxx
#include "sidl_BaseInterface.hxx"
#endif
#ifndef included_sidl_ClassInfo_hxx
#include "sidl_ClassInfo.hxx"
#endif
#ifndef included_sidl_NotImplementedException_hxx
#include "sidl_NotImplementedException.hxx"
#endif
// DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisShell._includes)
// Insert-Code-Here {MPQC.GaussianBasisShell._includes} (additional includes or code)
// DO-NOT-DELETE splicer.end(MPQC.GaussianBasisShell._includes)

// special constructor, used for data wrapping(required).  Do not put code here unless you really know what you're doing!
MPQC::GaussianBasisShell_impl::GaussianBasisShell_impl() : StubBase(
  reinterpret_cast< void*>(::MPQC::GaussianBasisShell::_wrapObj(
  reinterpret_cast< void*>(this))),false) , _wrapped(true){ 
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisShell._ctor2)
  // Insert-Code-Here {MPQC.GaussianBasisShell._ctor2} (ctor2)
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisShell._ctor2)
}

// user defined constructor
void MPQC::GaussianBasisShell_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisShell._ctor)
  // Insert-Code-Here {MPQC.GaussianBasisShell._ctor} (constructor)
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisShell._ctor)
}

// user defined destructor
void MPQC::GaussianBasisShell_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisShell._dtor)
  // Insert-Code-Here {MPQC.GaussianBasisShell._dtor} (destructor)
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisShell._dtor)
}

// static class initializer
void MPQC::GaussianBasisShell_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisShell._load)
  // Insert-Code-Here {MPQC.GaussianBasisShell._load} (class initialization)
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisShell._load)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 * Method:  initialize[]
 */
void
MPQC::GaussianBasisShell_impl::initialize_impl (
  /* in */void* scshell ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisShell.initialize)

  shell_ptr_ = static_cast<GaussianShell*>(scshell);
  sc_shell_.assign_pointer( shell_ptr_ );  
  if(sc_shell_.null())
    cout << "Shell: sc::GaussianShell is null" << endl;

  max_am_ = sc_shell_->max_angular_momentum();

  // determine angular type
  int has_pure = 0;
  int has_cartesian = 0;
  for(int i=0; i<sc_shell_->ncontraction(); ++i) {
    if( sc_shell_->is_cartesian(i) )
      ++has_cartesian;
    else
      ++has_pure;
  }
  if(has_pure && has_cartesian)
    angular_type_ = AngularType_MIXED;
  else if(has_pure)
    angular_type_ = AngularType_SPHERICAL;
  else if(has_cartesian)
    angular_type_ = AngularType_CARTESIAN;

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisShell.initialize)
}

/**
 *  Get the number of contractions in the shell. 
 * @return number of contractions 
 */
int32_t
MPQC::GaussianBasisShell_impl::get_n_contraction_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisShell.get_n_contraction)

  return sc_shell_->ncontraction();

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisShell.get_n_contraction)
}

/**
 *  Get the number of primitives in the shell.
 * @return number of primitives 
 */
int32_t
MPQC::GaussianBasisShell_impl::get_n_primitive_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisShell.get_n_primitive)

  return sc_shell_->nprimitive();

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisShell.get_n_primitive)
}

/**
 *  Get the coefficient for an unnormalized primitive 
 * in a contraction.
 * @param connum contraction number
 * @param expnum primitive number
 * @return contraction coefficient 
 */
double
MPQC::GaussianBasisShell_impl::get_contraction_coef_impl (
  /* in */int32_t connum,
  /* in */int32_t expnum ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisShell.get_contraction_coef)

  return sc_shell_->coefficient_unnorm(connum,expnum);

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisShell.get_contraction_coef)
}

/**
 *  Get the exponent for a primitive.
 * @param expnum primitive id number
 * @return exponent 
 */
double
MPQC::GaussianBasisShell_impl::get_exponent_impl (
  /* in */int32_t expnum ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisShell.get_exponent)

  return sc_shell_->exponent(expnum);

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisShell.get_exponent)
}

/**
 *  Get the angular momentum for a single contraction.
 * @param connum contraction id number
 * @return angular momentum value 
 */
int32_t
MPQC::GaussianBasisShell_impl::get_angular_momentum_impl (
  /* in */int32_t connum ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisShell.get_angular_momentum)

  return sc_shell_->am(connum);

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisShell.get_angular_momentum)
}

/**
 *  Get the max angular momentum, considering all contractions 
 * in the shell.
 * @return maximum angular momentum value 
 */
int32_t
MPQC::GaussianBasisShell_impl::get_max_angular_momentum_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisShell.get_max_angular_momentum)

  return max_am_;

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisShell.get_max_angular_momentum)
}

/**
 *  Get the angular type for a single contraction.
 * @param connum contraction number
 * @return enum AngularType 
 */
::Chemistry::QC::GaussianBasis::AngularType
MPQC::GaussianBasisShell_impl::get_contraction_angular_type_impl (
  /* in */int32_t connum ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisShell.get_contraction_angular_type)

  AngularType angular;
  
  if(sc_shell_->is_cartesian(connum) )
    angular = AngularType_CARTESIAN;
  else 
    angular = AngularType_SPHERICAL;
  
  return angular;

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisShell.get_contraction_angular_type)
}

/**
 *  Get the angular type.
 * @return enum AngularType 
 */
::Chemistry::QC::GaussianBasis::AngularType
MPQC::GaussianBasisShell_impl::get_angular_type_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisShell.get_angular_type)

  return angular_type_;

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisShell.get_angular_type)
}

/**
 *  Print the shell data. 
 */
void
MPQC::GaussianBasisShell_impl::print_shell_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisShell.print_shell)

 std::cout << "\n    shell:";
    std::cout << "\n      type: [";
    for(int icon=0; icon<get_n_contraction(); ++icon) 
      std::cout << " am = " << get_angular_momentum(icon);
    if( max_am_ > 1 ) {
      if( angular_type_ == AngularType_CARTESIAN )
         std::cout << " puream = 0";
      else if( angular_type_ == AngularType_SPHERICAL )
         std::cout << " puream = 1";
      else if( angular_type_ == AngularType_MIXED )
         std::cerr << " mixed angular types?";
    }
    std::cout << "]\n";
    // {exp coef:<am> ...} = {
    std::cout << "      exp";
    for(int icon=0; icon<get_n_contraction(); ++icon)
      std::cout << " coef:" << icon;
    std::cout << "\n";
    // <exp> <coef> ...
    for(int iprim=0; iprim<get_n_primitive(); ++iprim) {
      std::cout << "\t" << get_exponent(iprim);
      for(int icon=0; icon<get_n_contraction(); ++icon)
        std::cout << "\t" << get_contraction_coef(icon, iprim);
      std::cout << endl;
    }

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisShell.print_shell)
}


// DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisShell._misc)
// Insert-Code-Here {MPQC.GaussianBasisShell._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(MPQC.GaussianBasisShell._misc)

