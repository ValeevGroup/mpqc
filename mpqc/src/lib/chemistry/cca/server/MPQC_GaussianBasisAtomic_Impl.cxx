// 
// File:          MPQC_GaussianBasisAtomic_Impl.cxx
// Symbol:        MPQC.GaussianBasisAtomic-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.GaussianBasisAtomic
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 
#include "MPQC_GaussianBasisAtomic_Impl.hxx"

// 
// Includes for all method dependencies.
// 
#ifndef included_Chemistry_QC_GaussianBasis_AngularType_hxx
#include "Chemistry_QC_GaussianBasis_AngularType.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_ShellInterface_hxx
#include "Chemistry_QC_GaussianBasis_ShellInterface.hxx"
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
// DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisAtomic._includes)
// Insert-Code-Here {MPQC.GaussianBasisAtomic._includes} (additional includes or code)
// DO-NOT-DELETE splicer.end(MPQC.GaussianBasisAtomic._includes)

// special constructor, used for data wrapping(required).  Do not put code here unless you really know what you're doing!
MPQC::GaussianBasisAtomic_impl::GaussianBasisAtomic_impl() : StubBase(
  reinterpret_cast< void*>(::MPQC::GaussianBasisAtomic::_wrapObj(
  reinterpret_cast< void*>(this))),false) , _wrapped(true){ 
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisAtomic._ctor2)
  // Insert-Code-Here {MPQC.GaussianBasisAtomic._ctor2} (ctor2)
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisAtomic._ctor2)
}

// user defined constructor
void MPQC::GaussianBasisAtomic_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisAtomic._ctor)
  // Insert-Code-Here {MPQC.GaussianBasisAtomic._ctor} (constructor)
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisAtomic._ctor)
}

// user defined destructor
void MPQC::GaussianBasisAtomic_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisAtomic._dtor)

  // JPK: problems here
  //for(int i=0; i<nshell_; ++i)
  //  delete &shell_array_[i];
  //delete shell_array_;

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisAtomic._dtor)
}

// static class initializer
void MPQC::GaussianBasisAtomic_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisAtomic._load)
  // Insert-Code-Here {MPQC.GaussianBasisAtomic._load} (class initialization)
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisAtomic._load)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 * Method:  initialize[]
 */
void
MPQC::GaussianBasisAtomic_impl::initialize_impl (
  /* in */void* scbasis,
  /* in */int32_t atomnum ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisAtomic.initialize)

  atomnum_ = atomnum;
  
  gbs_ptr_ = static_cast<GaussianBasisSet*>(scbasis);
  sc_gbs_.assign_pointer( gbs_ptr_ );  
  if(sc_gbs_.null())
    cout << "Atomic: sc::GaussianBasisSet is null" << endl;

  // create shell array
  nshell_ = sc_gbs_->nshell_on_center(atomnum_);
  shell_array_ = new GaussianBasisShell[nshell_];
  for(int i=0; i<nshell_; ++i) {
    shell_array_[i] = GaussianBasisShell::_create();
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

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisAtomic.initialize)
}

/**
 *  Get the canonical basis set name. 
 * @return canonical name 
 */
::std::string
MPQC::GaussianBasisAtomic_impl::get_name_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisAtomic.get_name)

  return sc_gbs_->name();

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisAtomic.get_name)
}

/**
 *  Get the number of basis functions.
 * @return number of functions 
 */
int32_t
MPQC::GaussianBasisAtomic_impl::get_n_basis_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisAtomic.get_n_basis)

  return sc_gbs_->nbasis_on_center(atomnum_);

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisAtomic.get_n_basis)
}

/**
 *  Get the number of shells.
 * @return number of shells 
 */
int32_t
MPQC::GaussianBasisAtomic_impl::get_n_shell_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisAtomic.get_n_shell)

  return nshell_;

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisAtomic.get_n_shell)
}

/**
 *  Get the max angular momentum for any shell on the atom.
 * @return max angular momentum 
 */
int32_t
MPQC::GaussianBasisAtomic_impl::get_max_angular_momentum_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisAtomic.get_max_angular_momentum)

  return max_am_;

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisAtomic.get_max_angular_momentum)
}

/**
 *  Get the angular type for the atom.
 * @return enum AngularType 
 */
::Chemistry::QC::GaussianBasis::AngularType
MPQC::GaussianBasisAtomic_impl::get_angular_type_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisAtomic.get_angular_type)

  return angular_type_;

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisAtomic.get_angular_type)
}

/**
 *  Get a gaussian shell. 
 * @param shellnum shell number
 * @return Shell 
 */
::Chemistry::QC::GaussianBasis::ShellInterface
MPQC::GaussianBasisAtomic_impl::get_shell_impl (
  /* in */int32_t shellnum ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisAtomic.get_shell)

  return shell_array_[shellnum];

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisAtomic.get_shell)
}

/**
 *  Print the atomic basis data. 
 */
void
MPQC::GaussianBasisAtomic_impl::print_atomic_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisAtomic.print_atomic)

  std::cout << "\n  Atomic basis set:";
  for( int i=0; i<nshell_; ++i ) 
    shell_array_[i].print_shell();

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisAtomic.print_atomic)
}


// DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisAtomic._misc)
// Insert-Code-Here {MPQC.GaussianBasisAtomic._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(MPQC.GaussianBasisAtomic._misc)

