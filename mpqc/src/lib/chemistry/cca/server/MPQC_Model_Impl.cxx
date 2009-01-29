// 
// File:          MPQC_Model_Impl.cxx
// Symbol:        MPQC.Model-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.Model
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 
#include "MPQC_Model_Impl.hxx"

// 
// Includes for all method dependencies.
// 
#ifndef included_Chemistry_MoleculeInterface_hxx
#include "Chemistry_MoleculeInterface.hxx"
#endif
#ifndef included_gov_cca_TypeMap_hxx
#include "gov_cca_TypeMap.hxx"
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
// DO-NOT-DELETE splicer.begin(MPQC.Model._includes)

static int ex_counter = 0;

#include "ChemistryCXX_Molecule.hxx"
#include <iostream>
using namespace std;

sidl::array<double>
vector_to_array(const sc::RefSCVector &v);

sc::RefSCVector
array_to_vector(sidl::array<double>, const sc::RefSCVector &v);

sc::RefSCVector
array_to_vector(sidl::array<double>, const sc::RefSCVector &v, double conv);

sidl::array<double>
matrix_to_array(const sc::RefSymmSCMatrix &v);

sc::RefSymmSCMatrix
array_to_matrix(sidl::array<double>, const sc::RefSymmSCMatrix &v);

// DO-NOT-DELETE splicer.end(MPQC.Model._includes)

// speical constructor, used for data wrapping(required).  Do not put code here unless you really know what you're doing!
MPQC::Model_impl::Model_impl() : StubBase(reinterpret_cast< void*>(
  ::MPQC::Model::_wrapObj(reinterpret_cast< void*>(this))),false) , _wrapped(
  true){ 
  // DO-NOT-DELETE splicer.begin(MPQC.Model._ctor2)
  // Insert-Code-Here {MPQC.Model._ctor2} (ctor2)
  // DO-NOT-DELETE splicer.end(MPQC.Model._ctor2)
}

// user defined constructor
void MPQC::Model_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.Model._ctor)
  // Insert-Code-Here {MPQC.Model._ctor} (constructor)
  // DO-NOT-DELETE splicer.end(MPQC.Model._ctor)
}

// user defined destructor
void MPQC::Model_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.Model._dtor)
  // Insert-Code-Here {MPQC.Model._dtor} (destructor)
  // DO-NOT-DELETE splicer.end(MPQC.Model._dtor)
}

// static class initializer
void MPQC::Model_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.Model._load)
  // Insert-Code-Here {MPQC.Model._load} (class initialization)
  // DO-NOT-DELETE splicer.end(MPQC.Model._load)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 * Method:  initialize_parsedkeyval[]
 */
void
MPQC::Model_impl::initialize_parsedkeyval_impl (
  /* in */const ::std::string& keyword,
  /* in */const ::std::string& input ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.initialize_parsedkeyval)

  sc::Ref<sc::ParsedKeyVal> kv = new sc::ParsedKeyVal();
  //sc::ClassDesc::list_all_classes();
  kv->parse_string(input.c_str());
  sc::Ref<sc::DescribedClass> dc;
  try {
    dc = kv->describedclassvalue(keyword.c_str());
  }
  catch (std::exception &e) {
    e.what();
  }
  wfn_ = dynamic_cast<sc::Wavefunction*>(dc.pointer());

  sc::Molecule* scMol = wfn_->molecule().pointer();
  molecule_ = ChemistryCXX::Molecule::_create();
  molecule_.initialize( scMol->n_non_q_atom(), scMol->n_q_atom(), "bohr");

  // DO-NOT-DELETE splicer.end(MPQC.Model.initialize_parsedkeyval)
}

/**
 * Method:  initialize_parsedkeyval_file[]
 */
void
MPQC::Model_impl::initialize_parsedkeyval_file_impl (
  /* in */const ::std::string& keyword,
  /* in */const ::std::string& filename ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.initialize_parsedkeyval_file)

  std::cout << "reading " << keyword << " from " << filename << std::endl;
  sc::Ref<sc::ParsedKeyVal> kv = new sc::ParsedKeyVal(filename.c_str());
  sc::Ref<sc::DescribedClass> dc = kv->describedclassvalue(keyword.c_str());
  if (dc.null()) {
      std::cout << "WARNING: dc is null" << std::endl;
  }
  wfn_ = dynamic_cast<sc::Wavefunction*>(dc.pointer());
  if (wfn_.null()) {
      std::cout << "WARNING: wfn is null" << std::endl;
  }
  
  // DO-NOT-DELETE splicer.end(MPQC.Model.initialize_parsedkeyval_file)
}

/**
 * Method:  initialize_aggregatekeyval[]
 */
void
MPQC::Model_impl::initialize_aggregatekeyval_impl (
  /* in */const ::std::string& keyword,
  /* in */const ::std::string& input,
  /* in */void* describedclass ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.initialize_aggregatekeyval)

  // this doesn't seem to work, though JPK and CLJ both think it should

  std::cout << "Initializing MPQC model using aggregate keyval\n";
  
  sc::Ref<sc::ParsedKeyVal> kv = new sc::ParsedKeyVal;
  kv->parse_string(input.c_str());
  
  sc::Ref<sc::DescribedClass>* dc_ptr = 
    static_cast<sc::Ref<sc::DescribedClass>*>(describedclass);
  sc::Ref<sc::DescribedClass> dc = *dc_ptr;
  sc::Ref<sc::AssignedKeyVal> akv = new sc::AssignedKeyVal;
  akv->assign("model:integrals",dc);

  sc::Ref<sc::AggregateKeyVal> aggkv = new sc::AggregateKeyVal(kv,akv);

  sc::Ref<sc::DescribedClass> dc2 = aggkv->describedclassvalue(keyword.c_str());
  wfn_ = dynamic_cast<sc::Wavefunction*>(dc2.pointer());

  // DO-NOT-DELETE splicer.end(MPQC.Model.initialize_aggregatekeyval)
}

/**
 * Method:  initialize_pointer[]
 */
void
MPQC::Model_impl::initialize_pointer_impl (
  /* in */void* ptr ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.initialize_pointer)

  wfn_ = reinterpret_cast<sc::Wavefunction*>(ptr);

  // DO-NOT-DELETE splicer.end(MPQC.Model.initialize_pointer)
}

/**
 * Method:  get_energy[]
 */
double
MPQC::Model_impl::get_energy_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.get_energy)

 // if ((++ex_counter)%3 == 0) {
 //     sidl::BaseException e = sidl::BaseException::_create();
 //     e.setNote("Simulated Numerical Error");
 //     throw e;
 //   }
  return wfn_->energy();

  // DO-NOT-DELETE splicer.end(MPQC.Model.get_energy)
}

/**
 *  Sets the accuracy for subsequent energy calculations.
 * @param acc The new accuracy. 
 */
void
MPQC::Model_impl::set_energy_accuracy_impl (
  /* in */double acc ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.set_energy_accuracy)

  wfn_->set_desired_value_accuracy(acc);

  // DO-NOT-DELETE splicer.end(MPQC.Model.set_energy_accuracy)
}

/**
 *  Returns the accuracy to which the energy is already computed.
 * The result is undefined if the energy has not already
 * been computed.
 * @return The energy accuracy. 
 */
double
MPQC::Model_impl::get_energy_accuracy_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.get_energy_accuracy)

  return wfn_->actual_value_accuracy();

  // DO-NOT-DELETE splicer.end(MPQC.Model.get_energy_accuracy)
}

/**
 *  This allows a programmer to request that if any result
 * is computed,
 * then the energy is computed too.  This allows, say, for a request
 * for a gradient to cause the energy to be computed.  This computed
 * energy is cached and returned when the get_energy() member
 * is called.
 * @param doit Whether or not to compute the energy.
 */
void
MPQC::Model_impl::set_do_energy_impl (
  /* in */bool doit ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.set_do_energy)

  wfn_->do_value(1);

  // DO-NOT-DELETE splicer.end(MPQC.Model.set_do_energy)
}

/**
 *  Returns the Cartesian gradient.  
 */
::sidl::array<double>
MPQC::Model_impl::get_gradient_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.get_gradient)

  return vector_to_array(wfn_->gradient());

  // DO-NOT-DELETE splicer.end(MPQC.Model.get_gradient)
}

/**
 *  Sets the accuracy for subsequent gradient calculations
 * @param acc The new accuracy for gradients. 
 */
void
MPQC::Model_impl::set_gradient_accuracy_impl (
  /* in */double acc ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.set_gradient_accuracy)

  wfn_->set_desired_gradient_accuracy(acc);

  // DO-NOT-DELETE splicer.end(MPQC.Model.set_gradient_accuracy)
}

/**
 *  Returns the accuracy to which the gradient is already computed.
 * The result is undefined if the gradient has not already
 * been computed.
 * @return The current gradient accuracy. 
 */
double
MPQC::Model_impl::get_gradient_accuracy_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.get_gradient_accuracy)

  return wfn_->actual_gradient_accuracy();

  // DO-NOT-DELETE splicer.end(MPQC.Model.get_gradient_accuracy)
}

/**
 *  Returns the Cartesian Hessian. @return The Hessian. 
 */
::sidl::array<double>
MPQC::Model_impl::get_hessian_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.get_hessian)

  return matrix_to_array(wfn_->hessian());

  // DO-NOT-DELETE splicer.end(MPQC.Model.get_hessian)
}

/**
 *  Sets the accuracy for subsequent Hessian calculations.
 * @param acc The new accuracy for Hessians. 
 */
void
MPQC::Model_impl::set_hessian_accuracy_impl (
  /* in */double acc ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.set_hessian_accuracy)

  wfn_->set_desired_hessian_accuracy(acc);

  // DO-NOT-DELETE splicer.end(MPQC.Model.set_hessian_accuracy)
}

/**
 *  Returns the accuracy to which the Hessian is already computed.
 * The result is undefined if the Hessian has not already
 * been computed. 
 */
double
MPQC::Model_impl::get_hessian_accuracy_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.get_hessian_accuracy)

  return wfn_->actual_hessian_accuracy();

  // DO-NOT-DELETE splicer.end(MPQC.Model.get_hessian_accuracy)
}

/**
 *  Returns a Cartesian guess Hessian. 
 */
::sidl::array<double>
MPQC::Model_impl::get_guess_hessian_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.get_guess_hessian)

  sc::RefSymmSCMatrix hess(wfn_->dimension(), wfn_->matrixkit());
  wfn_->guess_hessian(hess);
  return matrix_to_array(hess);

  // DO-NOT-DELETE splicer.end(MPQC.Model.get_guess_hessian)
}

/**
 *  Sets the accuracy for subsequent guess Hessian calculations.
 * @param acc The new accuracy for guess Hessians. 
 */
void
MPQC::Model_impl::set_guess_hessian_accuracy_impl (
  /* in */double acc ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.set_guess_hessian_accuracy)
  // Insert-Code-Here {MPQC.Model.set_guess_hessian_accuracy} (set_guess_hessian_accuracy method)
  // 
  // This method has not been implemented
  // 
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "set_guess_hessian_accuracy");
    throw ex;
  // DO-NOT-DELETE splicer.end(MPQC.Model.set_guess_hessian_accuracy)
}

/**
 *  Returns the accuracy to which the guess Hessian is
 * already computed.  The result is undefined if the guess Hessian
 * has not already been computed.
 * @return The guess hessian accuracy.  
 */
double
MPQC::Model_impl::get_guess_hessian_accuracy_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.get_guess_hessian_accuracy)

  return DBL_EPSILON;

  // DO-NOT-DELETE splicer.end(MPQC.Model.get_guess_hessian_accuracy)
}

/**
 *  This should be called when the object is no longer needed.
 * No other members may be called after finalize. 
 */
int32_t
MPQC::Model_impl::finalize_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.finalize)
  // DO-NOT-DELETE splicer.end(MPQC.Model.finalize)
}

/**
 *  Set the molecule. @param molecule The new molecule. 
 */
void
MPQC::Model_impl::set_molecule_impl (
  /* in */::Chemistry::MoleculeInterface molecule ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.set_molecule)
  
  int i,j;
  double conv = molecule.get_units().convert_to("bohr");
  sc::Molecule* scMol = wfn_->molecule().pointer();
  int n_non_q = molecule_.get_n_atom();
  int n_q = molecule_.get_n_pcharge();

  for( i=0; i<n_non_q; ++i) {
     int nqid = scMol->non_q_atom(i);
     for( j=0; j<3; ++j) 
        scMol->r(nqid)[j] = molecule.get_cart_coor(i,j);
  }
  for( i=0; i<n_q; ++i) {
     int qid = scMol->q_atom(i);
     for( j=0; j<3; ++j)
        scMol->r(qid)[j] = molecule.get_pcharge_cart_coor(i,j);
  }
  wfn_->set_x(array_to_vector(molecule.get_coor(), 
              wfn_->matrixkit()->vector(wfn_->dimension()),1.0));

  return;

  // DO-NOT-DELETE splicer.end(MPQC.Model.set_molecule)
}

/**
 *  Returns the molecule.  @return The Molecule object. 
 */
::Chemistry::MoleculeInterface
MPQC::Model_impl::get_molecule_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.get_molecule)

  sc::Molecule* scMol = wfn_->molecule().pointer();
  scMol->print();
  
  int i,j;
  int n_non_q = molecule_.get_n_atom();
  int n_q = molecule_.get_n_pcharge();
  for( i=0; i<n_non_q; ++i ) {
    int nqid = scMol->non_q_atom(i);
    molecule_.set_atomic_number(i,scMol->Z(nqid));
    for( j=0; j<3; ++j )
      molecule_.set_cart_coor( i, j, scMol->r(nqid)[j] );
  }
  for( i=0; i<n_q; ++i ) {
    int qid = scMol->q_atom(i);
    for( j=0; j<3; ++j )
      molecule_.set_pcharge_cart_coor( i, j, scMol->r(qid)[j] );
    molecule_.set_point_charge( i, scMol->charge(qid) );
  } 

  return molecule_;

  // DO-NOT-DELETE splicer.end(MPQC.Model.get_molecule)
}

/**
 *  Sets the initial CQoS metadata typemap.
 * The model may augment this typemap.
 * @param typemap The initial typemap. 
 */
void
MPQC::Model_impl::set_metadata_impl (
  /* in */::gov::cca::TypeMap typemap ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.set_metadata)

  cqos_tm_ = typemap;

  // DO-NOT-DELETE splicer.end(MPQC.Model.set_metadata)
}

/**
 *  Returns CQoS metadata typemap.
 * @return Metadata typemap. 
 */
::gov::cca::TypeMap
MPQC::Model_impl::get_metadata_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Model.get_metadata)

  cqos_tm_.putInt("nelectron",wfn_->nelectron());
  cqos_tm_.putInt("AtomicOrbitals",wfn_->ao_dimension());
  cqos_tm_.putInt("OrthSymmOrbitals",wfn_->oso_dimension());
  return cqos_tm_;
    
  // DO-NOT-DELETE splicer.end(MPQC.Model.get_metadata)
}


// DO-NOT-DELETE splicer.begin(MPQC.Model._misc)

sidl::array<double>
vector_to_array(const sc::RefSCVector &v)
{
  sidl::array<double> a = sidl::array<double>::create1d(v.dim().n());

  for (int i=0,ai=a.lower(0); i<v.dim().n(); i++,ai++)
      a.set(ai, v(i));

  return a;
}

sc::RefSCVector
array_to_vector(sidl::array<double> a, const sc::RefSCVector &v)
{
  for (int i=0,ai=a.lower(0); i<v.dim().n(); i++,ai++) {
      v(i) = a.get(ai);
    }

  return v;
}

sc::RefSCVector
array_to_vector(sidl::array<double> a, const sc::RefSCVector &v, double conv)
{
  for (int i=0,ai=a.lower(0); i<v.dim().n(); i++,ai++) {
      v(i) = a.get(ai) * conv;
    }

  return v;
}

sidl::array<double>
matrix_to_array(const sc::RefSymmSCMatrix &v)
{
  sidl::array<double> a = sidl::array<double>::create2dCol(v.dim().n(), v.dim().n());

  for (int i=0,ai=a.lower(0); i<v.dim().n(); i++,ai++) {
      for (int j=0,aj=a.lower(1); j<v.dim().n(); j++,aj++) {
          a.set(ai, aj, v(i,j));
        }
    }

  return a;
}

sc::RefSymmSCMatrix
array_to_matrix(sidl::array<double> a, const sc::RefSymmSCMatrix &v)
{
  for (int i=0,ai=a.lower(0); i<v.dim().n(); i++,ai++) {
      for (int j=0,aj=a.lower(1); j<v.dim().n(); j++,aj++) {
          v(i,j) = a.get(ai,aj);
        }
    }

  return v;
}

// DO-NOT-DELETE splicer.end(MPQC.Model._misc)

