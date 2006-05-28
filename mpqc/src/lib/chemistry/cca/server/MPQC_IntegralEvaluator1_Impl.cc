// 
// File:          MPQC_IntegralEvaluator1_Impl.cc
// Symbol:        MPQC.IntegralEvaluator1-v0.2
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for MPQC.IntegralEvaluator1
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// 
#include "MPQC_IntegralEvaluator1_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1._includes)
// Insert-Code-Here {MPQC.IntegralEvaluator1._includes} (additional includes or code)
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1._includes)

// user-defined constructor.
void MPQC::IntegralEvaluator1_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1._ctor)
  // Insert-Code-Here {MPQC.IntegralEvaluator1._ctor} (constructor)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1._ctor)
}

// user-defined destructor.
void MPQC::IntegralEvaluator1_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1._dtor)
  // Insert-Code-Here {MPQC.IntegralEvaluator1._dtor} (destructor)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1._dtor)
}

// static class initializer.
void MPQC::IntegralEvaluator1_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1._load)
  // Insert-Code-Here {MPQC.IntegralEvaluator1._load} (class initialization)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Method:  add_evaluator[]
 */
void
MPQC::IntegralEvaluator1_impl::add_evaluator (
  /* in */ void* eval,
  /* in */ ::Chemistry::QC::GaussianBasis::IntegralDescr desc ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.add_evaluator)

  eval_.add_evaluator(eval,desc);

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.add_evaluator)
}

/**
 * Method:  set_basis[]
 */
void
MPQC::IntegralEvaluator1_impl::set_basis (
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1 ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.set_basis)

  basis_sets_.push_back(bs1);

  eval_.set_basis( basis_sets_ );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.set_basis)
}

/**
 * Method:  set_reorder[]
 */
void
MPQC::IntegralEvaluator1_impl::set_reorder (
  /* in */ int32_t reorder ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.set_reorder)

  eval_.set_reorder( reorder );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.set_reorder)
}

/**
 * Get buffer pointer for given type.
 * @return Buffer pointer. 
 */
void*
MPQC::IntegralEvaluator1_impl::get_buffer (
  /* in */ ::Chemistry::QC::GaussianBasis::IntegralDescr desc ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.get_buffer)
  
  return eval_.get_buffer( desc );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.get_buffer)
}

/**
 * Method:  get_deriv_centers[]
 */
::Chemistry::QC::GaussianBasis::DerivCenters
MPQC::IntegralEvaluator1_impl::get_deriv_centers ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.get_deriv_centers)

  return eval_.get_deriv_centers();

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.get_deriv_centers)
}

/**
 * Method:  get_descriptor[]
 */
::Chemistry::QC::GaussianBasis::CompositeIntegralDescr
MPQC::IntegralEvaluator1_impl::get_descriptor ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.get_descriptor)

  return eval_.get_descriptor();

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.get_descriptor)
}

/**
 * Compute a shell singlet of integrals.
 * @param shellnum1 Gaussian shell number 1.
 * @param deriv_level Derivative level. 
 */
void
MPQC::IntegralEvaluator1_impl::compute (
  /* in */ int64_t shellnum1 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.compute)

  computer_.set_shells( shellnum1 );
  eval_.compute( &computer_ );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.compute)
}

/**
 * Compute a shell singlet of integrals and return as a borrowed
 * sidl array.
 * @param shellnum1 Gaussian shell number 1.
 * @return Borrowed sidl array buffer. 
 */
::sidl::array<double>
MPQC::IntegralEvaluator1_impl::compute_array (
  /* in */ int64_t shellnum1 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.compute_array)

  computer_.set_shells( shellnum1 );
  return eval_.compute_array( &computer_ );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.compute_array)
}


// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1._misc)
// Insert-Code-Here {MPQC.IntegralEvaluator1._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1._misc)

