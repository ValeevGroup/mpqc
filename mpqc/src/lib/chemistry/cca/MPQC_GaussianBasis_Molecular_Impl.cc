// 
// File:          MPQC_GaussianBasis_Molecular_Impl.cc
// Symbol:        MPQC.GaussianBasis_Molecular-v0.2
// Symbol Type:   class
// Babel Version: 0.10.2
// Description:   Server-side implementation for MPQC.GaussianBasis_Molecular
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.2
// 
#include "MPQC_GaussianBasis_Molecular_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Molecular._includes)
#include <chemistry/molecule/molecule.h>
// DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Molecular._includes)

// user-defined constructor.
void MPQC::GaussianBasis_Molecular_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Molecular._ctor)
  // add construction details here
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Molecular._ctor)
}

// user-defined destructor.
void MPQC::GaussianBasis_Molecular_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Molecular._dtor)

  // JK: problems here
  //for(int i=0; i<natom_; ++i)
  //  delete &atomic_array_[i];
  //delete atomic_array_;

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Molecular._dtor)
}

// static class initializer.
void MPQC::GaussianBasis_Molecular_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Molecular._load)
  // guaranteed to be called at most once before any other method in this class
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Molecular._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Method:  initialize[]
 */
void
MPQC::GaussianBasis_Molecular_impl::initialize (
  /* in */ void* scbasis,
  /* in */ const ::std::string& label ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Molecular.initialize)
  
  label_ = label;

  gbs_ptr_ = static_cast< GaussianBasisSet* >(scbasis);
  sc_gbs_.assign_pointer( gbs_ptr_ );  
  if(sc_gbs_.null())
    cout << "Molecular: sc::GaussianBasisSet is null" << endl;

  // determine angular type
  int has_pure = 0;
  int has_cartesian = 0;
  for(int i=0; i<sc_gbs_->nshell(); ++i) {
    for(int j=0; j<sc_gbs_->shell(i).ncontraction(); ++j) {
      if( sc_gbs_->shell(i).is_cartesian(j) )
	++has_cartesian;
      if( sc_gbs_->shell(i).is_pure(j) )
	++has_pure;
    }
  }

  if(has_pure && has_cartesian)
    angular_type_ = AngularType_MIXED;
  else if(has_pure)
    angular_type_ = AngularType_SPHERICAL;
  else if(has_cartesian)
    angular_type_ = AngularType_CARTESIAN;

  // create a CCA molecule
  Ref<sc::Molecule> scmol = sc_gbs_->molecule();
  natom_ = scmol->natom();
  molecule_ = Chemistry_Molecule::_create();
  molecule_.initialize(natom_, "bohr");
  for( int i=0; i<natom_; ++i) {
    molecule_.set_atomic_number(i,scmol->Z(i));
    for( int j=0; j<3; ++j) 
      molecule_.set_cart_coor( i, j, scmol->r(i,j) );
  }

  // create array of atomic basis sets
  atomic_array_ = new MPQC::GaussianBasis_Atomic[sc_gbs_->ncenter()];
  for( int i=0; i<sc_gbs_->ncenter(); ++i) {
    atomic_array_[i] = MPQC::GaussianBasis_Atomic::_create();
    atomic_array_[i].initialize(sc_gbs_.pointer(),i);
  }
  
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Molecular.initialize)
}

/**
 * Method:  sc_gbs_pointer[]
 */
void*
MPQC::GaussianBasis_Molecular_impl::sc_gbs_pointer ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Molecular.sc_gbs_pointer)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Molecular.sc_gbs_pointer)
}

/**
 * Get the user specified name.
 * @return name 
 */
::std::string
MPQC::GaussianBasis_Molecular_impl::get_label ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Molecular.get_label)
  return label_;
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Molecular.get_label)
}

/**
 * Get the number of basis functions.
 * @return number of functions 
 */
int64_t
MPQC::GaussianBasis_Molecular_impl::get_n_basis ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Molecular.get_n_basis)
  return sc_gbs_->nbasis();
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Molecular.get_n_basis)
}

/**
 * Get the number of shells.
 * @return number of shells 
 */
int64_t
MPQC::GaussianBasis_Molecular_impl::get_n_shell ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Molecular.get_n_shell)
  return sc_gbs_->nshell();
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Molecular.get_n_shell)
}

/**
 * Get the max angular momentum for any contraction in the 
 * basis set.
 * @return max angular momentum 
 */
int32_t
MPQC::GaussianBasis_Molecular_impl::get_max_angular_momentum ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Molecular.get_max_angular_momentum)
  return sc_gbs_->max_angular_momentum();
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Molecular.get_max_angular_momentum)
}

/**
 * Get the angular type.
 * @return enum AngularType 
 */
::Chemistry::QC::GaussianBasis::AngularType
MPQC::GaussianBasis_Molecular_impl::get_angular_type ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Molecular.get_angular_type)
  return angular_type_;
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Molecular.get_angular_type)
}

/**
 * Get an atomic basis set.
 * @param atomnum atom number 
 * @return Atomic 
 */
::Chemistry::QC::GaussianBasis::Atomic
MPQC::GaussianBasis_Molecular_impl::get_atomic (
  /* in */ int64_t atomnum ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Molecular.get_atomic)
  return atomic_array_[atomnum];
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Molecular.get_atomic)
}

/**
 * Get the molecule.
 * @return Molecule 
 */
::Chemistry::Molecule
MPQC::GaussianBasis_Molecular_impl::get_molecule ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Molecular.get_molecule)
  return molecule_;
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Molecular.get_molecule)
}

/**
 * Print the molecular basis data. 
 */
void
MPQC::GaussianBasis_Molecular_impl::print_molecular ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Molecular.print_molecular)
  std::cout << "\nMolecular Basis Set:";
  for( int i=0; i<natom_; ++i) {
    atomic_array_[i].print_atomic();
  }
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Molecular.print_molecular)
}


// DO-NOT-DELETE splicer.begin(MPQC.GaussianBasis_Molecular._misc)
// Put miscellaneous code here
// DO-NOT-DELETE splicer.end(MPQC.GaussianBasis_Molecular._misc)

