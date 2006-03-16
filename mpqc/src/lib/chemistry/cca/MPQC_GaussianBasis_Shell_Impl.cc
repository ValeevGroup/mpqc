// 
// File:          MPQC_GaussianBasis_Shell_Impl.cc
// Symbol:        MPQC.GaussianBasis_Shell-v0.2
// Symbol Type:   class
// Babel Version: 0.10.2
// Description:   Server-side implementation for MPQC.GaussianBasis_Shell
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.2
// 
#include "MPQC_GaussianBasis_Shell_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Shell._includes)
// Put additional includes or other arbitrary code here...
// DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Shell._includes)

// user-defined constructor.
void MPQC::GaussianBasis_Shell_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Shell._ctor)
  // add construction details here
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Shell._ctor)
}

// user-defined destructor.
void MPQC::GaussianBasis_Shell_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Shell._dtor)
  // add destruction details here
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Shell._dtor)
}

// static class initializer.
void MPQC::GaussianBasis_Shell_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Shell._load)
  // guaranteed to be called at most once before any other method in this class
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Shell._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Method:  initialize[]
 */
void
MPQC::GaussianBasis_Shell_impl::initialize (
  /* in */ void* scshell ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Shell.initialize)

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
  
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Shell.initialize)
}

/**
 * Get the number of contractions in the shell. 
 * @return number of contractions 
 */
int32_t
MPQC::GaussianBasis_Shell_impl::get_n_contraction ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Shell.get_n_contraction)
  return sc_shell_->ncontraction();
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Shell.get_n_contraction)
}

/**
 * Get the number of primitives in the shell.
 * @return number of primitives 
 */
int32_t
MPQC::GaussianBasis_Shell_impl::get_n_primitive ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Shell.get_n_primitive)
  return sc_shell_->nprimitive();
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Shell.get_n_primitive)
}

/**
 * Get the coefficient for an unnormalized primitive 
 * in a contraction.
 * @param connum contraction number
 * @param expnum primitive number
 * @return contraction coefficient 
 */
double
MPQC::GaussianBasis_Shell_impl::get_contraction_coef (
  /* in */ int32_t connum,
  /* in */ int32_t expnum ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Shell.get_contraction_coef)
  return sc_shell_->coefficient_unnorm(connum,expnum);
  //return sc_shell_->coefficient_norm(connum,expnum);
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Shell.get_contraction_coef)
}

/**
 * Get the exponent for a primitive.
 * @param expnum primitive id number
 * @return exponent 
 */
double
MPQC::GaussianBasis_Shell_impl::get_exponent (
  /* in */ int32_t expnum ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Shell.get_exponent)
  return sc_shell_->exponent(expnum);
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Shell.get_exponent)
}

/**
 * Get the angular momentum for a single contraction.
 * @param connum contraction id number
 * @return angular momentum value 
 */
int32_t
MPQC::GaussianBasis_Shell_impl::get_angular_momentum (
  /* in */ int32_t connum ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Shell.get_angular_momentum)
  return sc_shell_->am(connum);
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Shell.get_angular_momentum)
}

/**
 * Get the max angular momentum, considering all contractions 
 * in the shell.
 * @return maximum angular momentum value 
 */
int32_t
MPQC::GaussianBasis_Shell_impl::get_max_angular_momentum ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Shell.get_max_angular_momentum)
  return max_am_;
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Shell.get_max_angular_momentum)
}

/**
 * Get the angular type for a single contraction.
 * @param connum contraction number
 * @return enum AngularType 
 */
::Chemistry::QC::GaussianBasis::AngularType
MPQC::GaussianBasis_Shell_impl::get_contraction_angular_type (
  /* in */ int32_t connum ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Shell.get_contraction_angular_type)

  AngularType angular;
  
  if(sc_shell_->is_cartesian(connum) )
    angular = AngularType_CARTESIAN;
  else 
    angular = AngularType_SPHERICAL;
  
  return angular;

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Shell.get_contraction_angular_type)
}

/**
 * Get the angular type.
 * @return enum AngularType 
 */
::Chemistry::QC::GaussianBasis::AngularType
MPQC::GaussianBasis_Shell_impl::get_angular_type ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Shell.get_angular_type)
  return angular_type_;
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Shell.get_angular_type)
}

/**
 * Print the shell data. 
 */
void
MPQC::GaussianBasis_Shell_impl::print_shell ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Shell.print_shell)
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
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Shell.print_shell)
}


// DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Shell._misc)
// Put miscellaneous code here
// DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Shell._misc)

