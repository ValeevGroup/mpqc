// 
// File:          MPQC_Chemistry_QC_Model_Impl.cc
// Symbol:        MPQC.Chemistry_QC_Model-v0.2
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for MPQC.Chemistry_QC_Model
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// 
#include "MPQC_Chemistry_QC_Model_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model._includes)

static int ex_counter = 0;

#include "Chemistry_Chemistry_Molecule.hh"
//#include <util/class/class.h>
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

// DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model._includes)

// user-defined constructor.
void MPQC::Chemistry_QC_Model_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model._ctor)
  // add construction details here
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model._ctor)
}

// user-defined destructor.
void MPQC::Chemistry_QC_Model_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model._dtor)
  // add destruction details here
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model._dtor)
}

// static class initializer.
void MPQC::Chemistry_QC_Model_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model._load)
  // guaranteed to be called at most once before any other method in this class
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Method:  initialize_parsedkeyval[]
 */
void
MPQC::Chemistry_QC_Model_impl::initialize_parsedkeyval (
  /* in */ const ::std::string& keyword,
  /* in */ const ::std::string& input ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model.initialize_parsedkeyval)

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

  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model.initialize_parsedkeyval)
}

/**
 * Method:  initialize_parsedkeyval_file[]
 */
void
MPQC::Chemistry_QC_Model_impl::initialize_parsedkeyval_file (
  /* in */ const ::std::string& keyword,
  /* in */ const ::std::string& filename ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model.initialize_parsedkeyval_file)
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
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model.initialize_parsedkeyval_file)
}

/**
 * Method:  initialize_aggregatekeyval[]
 */
void
MPQC::Chemistry_QC_Model_impl::initialize_aggregatekeyval (
  /* in */ const ::std::string& keyword,
  /* in */ const ::std::string& input,
  /* in */ void* describedclass ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model.initialize_aggregatekeyval)

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
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model.initialize_aggregatekeyval)
}

/**
 * Method:  initialize_pointer[]
 */
void
MPQC::Chemistry_QC_Model_impl::initialize_pointer (
  /* in */ void* ptr ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model.initialize_pointer)
  wfn_ = reinterpret_cast<sc::Wavefunction*>(ptr);
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model.initialize_pointer)
}

/**
 * Set the molecule. @param molecule The new molecule. 
 */
void
MPQC::Chemistry_QC_Model_impl::set_molecule (
  /* in */ ::Chemistry::Molecule molecule ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model.set_molecule)
  molecule_ = molecule;
  double conv = molecule_.get_units().convert_to("bohr");
  wfn_->molecule()->print();
  sc::Molecule* scMol = wfn_->molecule().pointer();
  for( int i=0; i<molecule_.get_n_atom(); ++i)
     for( int j=0; j<3; ++j) 
        scMol->r(i)[j] = molecule_.get_cart_coor(i,j)*conv;
  wfn_->set_x(array_to_vector(molecule_.get_coor(), 
              wfn_->matrixkit()->vector(wfn_->dimension()),conv));
  return;
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model.set_molecule)
}

/**
 * Returns the molecule.  @return The Molecule object. 
 */
::Chemistry::Molecule
MPQC::Chemistry_QC_Model_impl::get_molecule ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model.get_molecule)
  return molecule_;
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model.get_molecule)
}

/**
 * Method:  get_energy[]
 */
double
MPQC::Chemistry_QC_Model_impl::get_energy ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model.get_energy)
 // if ((++ex_counter)%3 == 0) {
 //     sidl::BaseException e = sidl::BaseException::_create();
 //     e.setNote("Simulated Numerical Error");
 //     throw e;
 //   }
  return wfn_->energy();
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model.get_energy)
}

/**
 * Sets the accuracy for subsequent energy calculations.
 * @param acc The new accuracy. 
 */
void
MPQC::Chemistry_QC_Model_impl::set_energy_accuracy (
  /* in */ double acc ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model.set_energy_accuracy)
  wfn_->set_desired_value_accuracy(acc);
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model.set_energy_accuracy)
}

/**
 * Returns the accuracy to which the energy is already computed.
 * The result is undefined if the energy has not already 
 * been computed.
 * @return The energy accuracy. 
 */
double
MPQC::Chemistry_QC_Model_impl::get_energy_accuracy ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model.get_energy_accuracy)
  return wfn_->actual_value_accuracy();
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model.get_energy_accuracy)
}

/**
 * This allows a programmer to request that if any result 
 * is computed,
 * then the energy is computed too.  This allows, say, for a request
 * for a gradient to cause the energy to be computed.  This computed
 * energy is cached and returned when the get_energy() member 
 * is called.
 * @param doit Whether or not to compute the energy.
 */
void
MPQC::Chemistry_QC_Model_impl::set_do_energy (
  /* in */ bool doit ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model.set_do_energy)
  wfn_->do_value(1);
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model.set_do_energy)
}

/**
 * Returns the Cartesian gradient.  
 */
::sidl::array<double>
MPQC::Chemistry_QC_Model_impl::get_gradient ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model.get_gradient)
  return vector_to_array(wfn_->gradient());
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model.get_gradient)
}

/**
 * Sets the accuracy for subsequent gradient calculations
 * @param acc The new accuracy for gradients. 
 */
void
MPQC::Chemistry_QC_Model_impl::set_gradient_accuracy (
  /* in */ double acc ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model.set_gradient_accuracy)
  wfn_->set_desired_gradient_accuracy(acc);
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model.set_gradient_accuracy)
}

/**
 * Returns the accuracy to which the gradient is already computed.
 * The result is undefined if the gradient has not already 
 * been computed.
 * @return The current gradient accuracy. 
 */
double
MPQC::Chemistry_QC_Model_impl::get_gradient_accuracy ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model.get_gradient_accuracy)
  return wfn_->actual_gradient_accuracy();
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model.get_gradient_accuracy)
}

/**
 * Returns the Cartesian Hessian. @return The Hessian. 
 */
::sidl::array<double>
MPQC::Chemistry_QC_Model_impl::get_hessian ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model.get_hessian)
  return matrix_to_array(wfn_->hessian());
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model.get_hessian)
}

/**
 * Sets the accuracy for subsequent Hessian calculations.
 * @param acc The new accuracy for Hessians. 
 */
void
MPQC::Chemistry_QC_Model_impl::set_hessian_accuracy (
  /* in */ double acc ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model.set_hessian_accuracy)
  wfn_->set_desired_hessian_accuracy(acc);
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model.set_hessian_accuracy)
}

/**
 * Returns the accuracy to which the Hessian is already computed.
 * The result is undefined if the Hessian has not already 
 * been computed. 
 */
double
MPQC::Chemistry_QC_Model_impl::get_hessian_accuracy ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model.get_hessian_accuracy)
  return wfn_->actual_hessian_accuracy();
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model.get_hessian_accuracy)
}

/**
 * Returns a Cartesian guess Hessian. 
 */
::sidl::array<double>
MPQC::Chemistry_QC_Model_impl::get_guess_hessian ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model.get_guess_hessian)
  sc::RefSymmSCMatrix hess(wfn_->dimension(), wfn_->matrixkit());
  wfn_->guess_hessian(hess);
  return matrix_to_array(hess);
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model.get_guess_hessian)
}

/**
 * Sets the accuracy for subsequent guess Hessian calculations.
 * @param acc The new accuracy for guess Hessians. 
 */
void
MPQC::Chemistry_QC_Model_impl::set_guess_hessian_accuracy (
  /* in */ double acc ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model.set_guess_hessian_accuracy)
  // noop
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model.set_guess_hessian_accuracy)
}

/**
 * Returns the accuracy to which the guess Hessian is 
 * already computed.  The result is undefined if the guess Hessian 
 * has not already been computed.
 * @return The guess hessian accuracy.  
 */
double
MPQC::Chemistry_QC_Model_impl::get_guess_hessian_accuracy ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model.get_guess_hessian_accuracy)
  return DBL_EPSILON;
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model.get_guess_hessian_accuracy)
}

/**
 * This can be called when this Model object is no longer needed.  
 * No other members may be called after finalize. 
 */
int32_t
MPQC::Chemistry_QC_Model_impl::finalize ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model.finalize)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model.finalize)
}


// DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_Model._misc)

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

// DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_Model._misc)

