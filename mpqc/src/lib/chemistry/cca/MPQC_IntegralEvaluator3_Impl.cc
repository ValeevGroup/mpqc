// 
// File:          MPQC_IntegralEvaluator3_Impl.cc
// Symbol:        MPQC.IntegralEvaluator3-v0.2
// Symbol Type:   class
// Babel Version: 0.9.8
// Description:   Server-side implementation for MPQC.IntegralEvaluator3
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.9.8
// 
#include "MPQC_IntegralEvaluator3_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3._includes)
// Put additional includes or other arbitrary code here...
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3._includes)

// user defined constructor
void MPQC::IntegralEvaluator3_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3._ctor)
  // add construction details here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3._ctor)
}

// user defined destructor
void MPQC::IntegralEvaluator3_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3._dtor)
  // add destruction details here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3._dtor)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 * Method:  set_integral_package[]
 */
void
MPQC::IntegralEvaluator3_impl::set_integral_package (
  /*in*/ const ::std::string& label ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.set_integral_package)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.set_integral_package)
}

/**
 * Method:  initialize_by_name[]
 */
void
MPQC::IntegralEvaluator3_impl::initialize_by_name (
  /*in*/ ::Chemistry::Molecule molecule,
  /*in*/ const ::std::string& basis_name,
  /*in*/ const ::std::string& evaluator_label,
  /*in*/ int64_t max_deriv ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.initialize_by_name)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.initialize_by_name)
}

/**
 * Method:  initialize[]
 */
void
MPQC::IntegralEvaluator3_impl::initialize (
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs3,
  /*in*/ const ::std::string& evaluator_label,
  /*in*/ int64_t max_deriv ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.initialize)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.initialize)
}

/**
 * Method:  buffer[]
 */
void*
MPQC::IntegralEvaluator3_impl::buffer () 
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.buffer)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.buffer)
}

/**
 * Method:  compute[]
 */
void
MPQC::IntegralEvaluator3_impl::compute (
  /*in*/ int64_t shellnum1,
  /*in*/ int64_t shellnum2,
  /*in*/ int64_t shellnum3,
  /*in*/ int64_t deriv_level ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.compute)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.compute)
}

/**
 * Method:  compute_array[]
 */
::sidl::array<double>
MPQC::IntegralEvaluator3_impl::compute_array (
  /*in*/ int64_t shellnum1,
  /*in*/ int64_t shellnum2,
  /*in*/ int64_t shellnum3,
  /*in*/ int64_t deriv_level ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.compute_array)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.compute_array)
}


// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3._misc)
// Put miscellaneous code here
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3._misc)

