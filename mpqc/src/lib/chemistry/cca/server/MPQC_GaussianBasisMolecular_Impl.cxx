// 
// File:          MPQC_GaussianBasisMolecular_Impl.cxx
// Symbol:        MPQC.GaussianBasisMolecular-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.GaussianBasisMolecular
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 
#include "MPQC_GaussianBasisMolecular_Impl.hxx"

// 
// Includes for all method dependencies.
// 
#ifndef included_Chemistry_MoleculeInterface_hxx
#include "Chemistry_MoleculeInterface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_AngularType_hxx
#include "Chemistry_QC_GaussianBasis_AngularType.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_AtomicInterface_hxx
#include "Chemistry_QC_GaussianBasis_AtomicInterface.hxx"
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
// DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisMolecular._includes)

#include <chemistry/molecule/molecule.h>

using namespace std;
using namespace Chemistry::QC::GaussianBasis;
using namespace Chemistry;
using namespace sc;

// DO-NOT-DELETE splicer.end(MPQC.GaussianBasisMolecular._includes)

// speical constructor, used for data wrapping(required).  Do not put code here unless you really know what you're doing!
MPQC::GaussianBasisMolecular_impl::GaussianBasisMolecular_impl() : StubBase(
  reinterpret_cast< void*>(::MPQC::GaussianBasisMolecular::_wrapObj(
  reinterpret_cast< void*>(this))),false) , _wrapped(true){ 
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisMolecular._ctor2)
  // Insert-Code-Here {MPQC.GaussianBasisMolecular._ctor2} (ctor2)
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisMolecular._ctor2)
}

// user defined constructor
void MPQC::GaussianBasisMolecular_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisMolecular._ctor)
  // Insert-Code-Here {MPQC.GaussianBasisMolecular._ctor} (constructor)
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisMolecular._ctor)
}

// user defined destructor
void MPQC::GaussianBasisMolecular_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisMolecular._dtor)

  // JK: problems here
  //for(int i=0; i<natom_; ++i)
  //  delete &atomic_array_[i];
  //delete atomic_array_;

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisMolecular._dtor)
}

// static class initializer
void MPQC::GaussianBasisMolecular_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisMolecular._load)
  // Insert-Code-Here {MPQC.GaussianBasisMolecular._load} (class initialization)
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisMolecular._load)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 * Method:  initialize[]
 */
void
MPQC::GaussianBasisMolecular_impl::initialize_impl (
  /* in */void* scbasis,
  /* in */const ::std::string& label ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisMolecular.initialize)

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
  molecule_ = ChemistryCXX::Molecule::_create();
  molecule_.initialize(natom_,0,"bohr");
  for( int i=0; i<natom_; ++i) {
    molecule_.set_atomic_number(i,scmol->Z(i));
    for( int j=0; j<3; ++j) 
      molecule_.set_cart_coor( i, j, scmol->r(i,j) );
  }

  // create array of atomic basis sets
  atomic_array_ = new MPQC::GaussianBasisAtomic[sc_gbs_->ncenter()];
  for( int i=0; i<sc_gbs_->ncenter(); ++i) {
    atomic_array_[i] = MPQC::GaussianBasisAtomic::_create();
    atomic_array_[i].initialize(sc_gbs_.pointer(),i);
  }

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisMolecular.initialize)
}

/**
 * Method:  sc_gbs_pointer[]
 */
void*
MPQC::GaussianBasisMolecular_impl::sc_gbs_pointer_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisMolecular.sc_gbs_pointer)
  // Insert-Code-Here {MPQC.GaussianBasisMolecular.sc_gbs_pointer} (sc_gbs_pointer method)
  // 
  // This method has not been implemented
  // 
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "sc_gbs_pointer");
    throw ex;
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisMolecular.sc_gbs_pointer)
}

/**
 *  Get the user specified name.
 * @return name 
 */
::std::string
MPQC::GaussianBasisMolecular_impl::get_label_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisMolecular.get_label)

  return label_;

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisMolecular.get_label)
}

/**
 *  Get the number of basis functions.
 * @return number of functions 
 */
int64_t
MPQC::GaussianBasisMolecular_impl::get_n_basis_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisMolecular.get_n_basis)

  return sc_gbs_->nbasis();

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisMolecular.get_n_basis)
}

/**
 *  Get the number of shells.
 * @return number of shells 
 */
int64_t
MPQC::GaussianBasisMolecular_impl::get_n_shell_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisMolecular.get_n_shell)

  return sc_gbs_->nshell();

  // Insert-Code-Here {MPQC.GaussianBasisMolecular.get_n_shell} (get_n_shell method)
  // 
  // This method has not been implemented
  // 
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "get_n_shell");
    throw ex;
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisMolecular.get_n_shell)
}

/**
 *  Get the max angular momentum for any contraction in the 
 * basis set.
 * @return max angular momentum 
 */
int32_t
MPQC::GaussianBasisMolecular_impl::get_max_angular_momentum_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisMolecular.get_max_angular_momentum)

  return sc_gbs_->max_angular_momentum();

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisMolecular.get_max_angular_momentum)
}

/**
 *  Get the angular type.
 * @return enum AngularType 
 */
::Chemistry::QC::GaussianBasis::AngularType
MPQC::GaussianBasisMolecular_impl::get_angular_type_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisMolecular.get_angular_type)

  return angular_type_;

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisMolecular.get_angular_type)
}

/**
 *  Get an atomic basis set.
 * @param atomnum atom number 
 * @return Atomic 
 */
::Chemistry::QC::GaussianBasis::AtomicInterface
MPQC::GaussianBasisMolecular_impl::get_atomic_impl (
  /* in */int64_t atomnum ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisMolecular.get_atomic)

  return atomic_array_[atomnum];

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisMolecular.get_atomic)
}

/**
 *  Get the molecule.
 * @return Molecule 
 */
::Chemistry::MoleculeInterface
MPQC::GaussianBasisMolecular_impl::get_molecule_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisMolecular.get_molecule)

  return molecule_;

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisMolecular.get_molecule)
}

/**
 *  Print the molecular basis data. 
 */
void
MPQC::GaussianBasisMolecular_impl::print_molecular_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisMolecular.print_molecular)

  std::cout << "\nMolecular Basis Set:";
  for( int i=0; i<natom_; ++i) {
    atomic_array_[i].print_atomic();
  }

  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisMolecular.print_molecular)
}


// DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisMolecular._misc)
// Insert-Code-Here {MPQC.GaussianBasisMolecular._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(MPQC.GaussianBasisMolecular._misc)

