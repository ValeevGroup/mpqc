// 
// File:          MPQC_IntegralEvaluator3_Impl.cc
// Symbol:        MPQC.IntegralEvaluator3-v0.2
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for MPQC.IntegralEvaluator3
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// xml-url       = /home/jpkenny/src/mpqc-libint2.build-shared/src/lib/chemistry/cca/server/../../../../../lib/cca/repo/MPQC.IntegralEvaluator3-v0.2.xml
// 
#include "MPQC_IntegralEvaluator3_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3._includes)
// Put additional includes or other arbitrary code here...
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3._includes)

// user-defined constructor.
void MPQC::IntegralEvaluator3_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3._ctor)
  // add construction details here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3._ctor)
}

// user-defined destructor.
void MPQC::IntegralEvaluator3_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3._dtor)
  // add destruction details here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3._dtor)
}

// static class initializer.
void MPQC::IntegralEvaluator3_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3._load)
  // guaranteed to be called at most once before any other method in this class
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Method:  add_evaluator[]
 */
void
MPQC::IntegralEvaluator3_impl::add_evaluator (
  /* in */ void* eval,
  /* in */ ::Chemistry::QC::GaussianBasis::IntegralDescr desc ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.add_evaluator)

  eval_.add_evaluator(&eval,desc);

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.add_evaluator)
}

/**
 * Method:  set_basis[]
 */
void
MPQC::IntegralEvaluator3_impl::set_basis (
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs3 ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.set_basis)
 
  basis_sets_.push_back(bs1);
  basis_sets_.push_back(bs2);
  basis_sets_.push_back(bs3);

  eval_.set_basis( basis_sets_ );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.set_basis)
}

/**
 * Method:  set_reorder[]
 */
void
MPQC::IntegralEvaluator3_impl::set_reorder (
  /* in */ int32_t reorder ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.set_reorder)

  eval_.set_reorder( reorder );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.set_reorder)
}

/**
 * Get buffer pointer for given type.
 * @return Buffer pointer. 
 */
void*
MPQC::IntegralEvaluator3_impl::get_buffer (
  /* in */ ::Chemistry::QC::GaussianBasis::IntegralDescr desc ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.get_buffer)

  return eval_.get_buffer( desc );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.get_buffer)
}

/**
 * Method:  get_deriv_centers[]
 */
::Chemistry::QC::GaussianBasis::DerivCenters
MPQC::IntegralEvaluator3_impl::get_deriv_centers ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.get_deriv_centers)

  return eval_.get_deriv_centers();

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.get_deriv_centers)
}

/**
 * Method:  get_descriptor[]
 */
::Chemistry::QC::GaussianBasis::CompositeIntegralDescr
MPQC::IntegralEvaluator3_impl::get_descriptor ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.get_descriptor)

  return eval_.get_descriptor();

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.get_descriptor)
}

/**
 * Compute a shell triplet of integrals.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @param shellnum3 Gaussian shell number 3.
 * @param deriv_level Derivative level. 
 */
void
MPQC::IntegralEvaluator3_impl::compute (
  /* in */ int64_t shellnum1,
  /* in */ int64_t shellnum2,
  /* in */ int64_t shellnum3 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.compute)

  computer_.set_shells( shellnum1, shellnum2, shellnum3 );
  eval_.compute( &computer_ );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.compute)
}

/**
 * Compute a shell triplet of integrals and return as a borrowed
 * sidl array.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @param shellnum3 Gaussian shell number 3.
 * @return Borrowed sidl array buffer. 
 */
::sidl::array<double>
MPQC::IntegralEvaluator3_impl::compute_array (
  /* in */ int64_t shellnum1,
  /* in */ int64_t shellnum2,
  /* in */ int64_t shellnum3 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.compute_array)
  
  computer_.set_shells( shellnum1, shellnum2, shellnum3 );
  eval_.compute( &computer_ );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.compute_array)
}


// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3._misc)
// Put miscellaneous code here
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3._misc)

