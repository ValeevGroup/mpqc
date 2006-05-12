// 
// File:          MPQC_GaussianBasis_Atomic_Impl.cc
// Symbol:        MPQC.GaussianBasis_Atomic-v0.2
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for MPQC.GaussianBasis_Atomic
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// xml-url       = /home/jpkenny/src/mpqc-libint2.build-shared/src/lib/chemistry/cca/server/../../../../../lib/cca/repo/MPQC.GaussianBasis_Atomic-v0.2.xml
// 
#include "MPQC_GaussianBasis_Atomic_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Atomic._includes)
// Put additional includes or other arbitrary code here...
// DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Atomic._includes)

// user-defined constructor.
void MPQC::GaussianBasis_Atomic_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Atomic._ctor)
  // add construction details here
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Atomic._ctor)
}

// user-defined destructor.
void MPQC::GaussianBasis_Atomic_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Atomic._dtor)

  // JPK: problems here
  //for(int i=0; i<nshell_; ++i)
  //  delete &shell_array_[i];
  //delete shell_array_;

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Atomic._dtor)
}

// static class initializer.
void MPQC::GaussianBasis_Atomic_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Atomic._load)
  // guaranteed to be called at most once before any other method in this class
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Atomic._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Method:  initialize[]
 */
void
MPQC::GaussianBasis_Atomic_impl::initialize (
  /* in */ void* scbasis,
  /* in */ int32_t atomnum ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Atomic.initialize)
  
  atomnum_ = atomnum;
  
  gbs_ptr_ = static_cast<GaussianBasisSet*>(scbasis);
  sc_gbs_.assign_pointer( gbs_ptr_ );  
  if(sc_gbs_.null())
    cout << "Atomic: sc::GaussianBasisSet is null" << endl;

  // create shell array
  nshell_ = sc_gbs_->nshell_on_center(atomnum_);
  shell_array_ = new GaussianBasis_Shell[nshell_];
  for(int i=0; i<nshell_; ++i) {
    shell_array_[i] = GaussianBasis_Shell::_create();
    GaussianShell &shell_ref = sc_gbs_->shell(atomnum_,i);
    shell_array_[i].initialize( &shell_ref );
  }

  // determine max am
  max_am_ = 0;
  int temp_am;
  for(int i=0; i<nshell_; ++i) {
    for(int j=0; j<shell_array_[i].get_n_contraction(); ++j) {
      temp_am = shell_array_[i].get_angular_momentum(j);
      if( temp_am > max_am_ )
	max_am_ = temp_am;
    }
  } 

  // determine angular type
  int has_pure = 0;
  int has_cartesian = 0;
  for(int i=0; i<sc_gbs_->nshell_on_center(atomnum_); ++i) {
    for(int j=0; j<sc_gbs_->shell(atomnum_,i).ncontraction(); ++j) {
      if( sc_gbs_->shell(atomnum_,i).is_cartesian(j) )
	++has_cartesian;
      if( sc_gbs_->shell(atomnum_,i).is_pure(j) )
	++has_pure;
    }
  }
  if(has_pure && has_cartesian)
    angular_type_ = AngularType_MIXED;
  else if(has_pure)
    angular_type_ = AngularType_SPHERICAL;
  else if(has_cartesian)
    angular_type_ = AngularType_CARTESIAN;

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Atomic.initialize)
}

/**
 * Get the canonical basis set name. 
 * @return canonical name 
 */
::std::string
MPQC::GaussianBasis_Atomic_impl::get_name ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Atomic.get_name)
  return sc_gbs_->name();
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Atomic.get_name)
}

/**
 * Get the number of basis functions.
 * @return number of functions 
 */
int32_t
MPQC::GaussianBasis_Atomic_impl::get_n_basis ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Atomic.get_n_basis)
  return sc_gbs_->nbasis_on_center(atomnum_);
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Atomic.get_n_basis)
}

/**
 * Get the number of shells.
 * @return number of shells 
 */
int32_t
MPQC::GaussianBasis_Atomic_impl::get_n_shell ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Atomic.get_n_shell)
  return nshell_;
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Atomic.get_n_shell)
}

/**
 * Get the max angular momentum for any shell on the atom.
 * @return max angular momentum 
 */
int32_t
MPQC::GaussianBasis_Atomic_impl::get_max_angular_momentum ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Atomic.get_max_angular_momentum)
  return max_am_;
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Atomic.get_max_angular_momentum)
}

/**
 * Get the angular type for the atom.
 * @return enum AngularType 
 */
::Chemistry::QC::GaussianBasis::AngularType
MPQC::GaussianBasis_Atomic_impl::get_angular_type ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Atomic.get_angular_type)
  return angular_type_;
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Atomic.get_angular_type)
}

/**
 * Get a gaussian shell. 
 * @param shellnum shell number
 * @return Shell 
 */
::Chemistry::QC::GaussianBasis::Shell
MPQC::GaussianBasis_Atomic_impl::get_shell (
  /* in */ int32_t shellnum ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Atomic.get_shell)
  return shell_array_[shellnum];
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Atomic.get_shell)
}

/**
 * Print the atomic basis data. 
 */
void
MPQC::GaussianBasis_Atomic_impl::print_atomic ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Atomic.print_atomic)
  std::cout << "\n  Atomic basis set:";
  for( int i=0; i<nshell_; ++i ) 
    shell_array_[i].print_shell();
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Atomic.print_atomic)
}


// DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Atomic._misc)
// Put miscellaneous code here
// DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Atomic._misc)

