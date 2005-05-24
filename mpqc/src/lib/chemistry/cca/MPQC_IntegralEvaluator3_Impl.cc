// 
// File:          MPQC_IntegralEvaluator3_Impl.cc
// Symbol:        MPQC.IntegralEvaluator3-v0.2
// Symbol Type:   class
// Babel Version: 0.10.2
// Description:   Server-side implementation for MPQC.IntegralEvaluator3
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.2
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
 * Method:  set_integral_package[]
 */
void
MPQC::IntegralEvaluator3_impl::set_integral_package (
  /* in */ const ::std::string& label ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.set_integral_package)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.set_integral_package)
}

/**
 * Initialize the evaluator.
 * @param bs1 Molecular basis on center 1.
 * @param bs2 Molecular basis on center 2.
 * @param bs3 Molecular basis on center 3.
 * @param label String specifying integral type.
 * @param max_deriv Max derivative to compute. 
 */
void
MPQC::IntegralEvaluator3_impl::initialize (
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs3,
  /* in */ const ::std::string& label,
  /* in */ int64_t max_deriv ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.initialize)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.initialize)
}

/**
 * Get the buffer pointer
 * @return Buffer pointer 
 */
void*
MPQC::IntegralEvaluator3_impl::get_buffer ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.get_buffer)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.get_buffer)
}

/**
 * Compute a shell triplet of integrals.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @param shellnum3 Gaussian shell number 3.
 * @param deriv_level Derivative level. 
 * @param deriv_ctr Derivative center descriptor. 
 */
void
MPQC::IntegralEvaluator3_impl::compute (
  /* in */ int64_t shellnum1,
  /* in */ int64_t shellnum2,
  /* in */ int64_t shellnum3,
  /* in */ int64_t deriv_level,
  /* in */ ::Chemistry::QC::GaussianBasis::DerivCenters deriv_ctr ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.compute)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.compute)
}

/**
 * Compute a shell triplet of integrals and return as a borrowed
 * sidl array.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @param shellnum3 Gaussian shell number 3.
 * @param deriv_level Derivative level.
 * @param deriv_ctr Derivative center desctiptor.
 * @return Borrowed sidl array buffer. 
 */
::sidl::array<double>
MPQC::IntegralEvaluator3_impl::compute_array (
  /* in */ int64_t shellnum1,
  /* in */ int64_t shellnum2,
  /* in */ int64_t shellnum3,
  /* in */ int64_t deriv_level,
  /* in */ ::Chemistry::QC::GaussianBasis::DerivCenters deriv_ctr ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.compute_array)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.compute_array)
}


// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3._misc)
// Put miscellaneous code here
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3._misc)

