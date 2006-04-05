// 
// File:          MPQC_IntegralEvaluator2_Impl.cc
// Symbol:        MPQC.IntegralEvaluator2-v0.2
// Symbol Type:   class
// Babel Version: 0.10.2
// Description:   Server-side implementation for MPQC.IntegralEvaluator2
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.2
// 
#include "MPQC_IntegralEvaluator2_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2._includes)
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2._includes)

// user-defined constructor.
void MPQC::IntegralEvaluator2_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2._ctor)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2._ctor)
}

// user-defined destructor.
void MPQC::IntegralEvaluator2_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2._dtor)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2._dtor)
}

// static class initializer.
void MPQC::IntegralEvaluator2_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2._load)
  // guaranteed to be called at most once before any other method in this class
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Method:  add_evaluator[]
 */
void
MPQC::IntegralEvaluator2_impl::add_evaluator (
  /* in */ void* eval,
  /* in */ ::Chemistry::QC::GaussianBasis::IntegralDescr desc ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2.add_evaluator)  

  eval_.add_evaluator(&eval,desc);

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2.add_evaluator)
}

/**
 * Method:  set_basis[]
 */
void
MPQC::IntegralEvaluator2_impl::set_basis (
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2 ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2.set_basis)

  basis_sets_.push_back(bs1);
  basis_sets_.push_back(bs2);

  eval_.set_basis( basis_sets_ );
  
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2.set_basis)
}

/**
 * Method:  set_reorder[]
 */
void
MPQC::IntegralEvaluator2_impl::set_reorder (
  /* in */ int32_t reorder ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2.set_reorder)

  eval_.set_reorder( reorder );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2.set_reorder)
}

/**
 * Get buffer pointer for given type.
 * @return Buffer pointer. 
 */
void*
MPQC::IntegralEvaluator2_impl::get_buffer (
  /* in */ ::Chemistry::QC::GaussianBasis::IntegralDescr desc ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2.get_buffer)

  return eval_.get_buffer( desc );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2.get_buffer)
}

/**
 * Method:  get_deriv_centers[]
 */
::Chemistry::QC::GaussianBasis::DerivCenters
MPQC::IntegralEvaluator2_impl::get_deriv_centers ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2.get_deriv_centers)

  return eval_.get_deriv_centers();

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2.get_deriv_centers)
}

/**
 * Method:  get_descriptor[]
 */
::Chemistry::QC::GaussianBasis::CompositeIntegralDescr
MPQC::IntegralEvaluator2_impl::get_descriptor ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2.get_descriptor)

  return eval_.get_descriptor();

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2.get_descriptor)
}

/**
 * Compute a shell doublet of integrals.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @param deriv_level Derivative level. 
 */
void
MPQC::IntegralEvaluator2_impl::compute (
  /* in */ int64_t shellnum1,
  /* in */ int64_t shellnum2 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2.compute)
  
  computer_.set_shells( shellnum1, shellnum2 );
  eval_.compute( &computer_ );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2.compute)
}

/**
 * Compute a shell doublet of integrals and return as a borrowed
 * sidl array.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @return Borrowed sidl array buffer. 
 */
::sidl::array<double>
MPQC::IntegralEvaluator2_impl::compute_array (
  /* in */ int64_t shellnum1,
  /* in */ int64_t shellnum2 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2.compute_array)
  
  computer_.set_shells( shellnum1, shellnum2 );
  return eval_.compute_array( &computer_ ); 

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2.compute_array)
}


// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2._misc)
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2._misc)

