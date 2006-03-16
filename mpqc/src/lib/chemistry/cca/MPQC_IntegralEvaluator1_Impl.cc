// 
// File:          MPQC_IntegralEvaluator1_Impl.cc
// Symbol:        MPQC.IntegralEvaluator1-v0.2
// Symbol Type:   class
// Babel Version: 0.10.2
// Description:   Server-side implementation for MPQC.IntegralEvaluator1
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.2
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
 * Method:  set_integral_package[]
 */
void
MPQC::IntegralEvaluator1_impl::set_integral_package (
  /* in */ ::Chemistry::QC::GaussianBasis::Package type ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.set_integral_package)
  // Insert-Code-Here {MPQC.IntegralEvaluator1.set_integral_package} (set_integral_package method)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.set_integral_package)
}

/**
 * Initialize the evaluator.
 * @param bs1 Molecular basis on center 1.
 * @param type ObIntEvalType specifying eval type.
 * @param max_deriv Max derivative to compute.
 * @param storage Available storage in bytes.
 * @param deriv_ctr Derivative center descriptor. 
 */
void
MPQC::IntegralEvaluator1_impl::obint_initialize (
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::ObIntEvalType type,
  /* in */ int32_t max_deriv,
  /* in */ int64_t storage,
  /* in */ ::Chemistry::QC::GaussianBasis::DerivCenters deriv_ctr ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.obint_initialize)
  // Insert-Code-Here {MPQC.IntegralEvaluator1.obint_initialize} (obint_initialize method)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.obint_initialize)
}

/**
 * Get one body int buffer pointer.
 * @return Buffer pointer. 
 */
void*
MPQC::IntegralEvaluator1_impl::get_obint_buffer ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.get_obint_buffer)
  // Insert-Code-Here {MPQC.IntegralEvaluator1.get_obint_buffer} (get_obint_buffer method)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.get_obint_buffer)
}

/**
 * Compute a shell singlet of integrals.  deriv_atom must
 * be used for nuclear derivatives if the operator contains
 * nuclear coordinates, otherwise, set to -1 and use deriv_ctr.
 * @param shellnum1 Gaussian shell number 1.
 * @param deriv_level Derivative level.
 * @param deriv_atom Atom number for derivative
 * (-1 if using DerivCenter). 
 */
void
MPQC::IntegralEvaluator1_impl::compute (
  /* in */ int64_t shellnum1,
  /* in */ int32_t deriv_level,
  /* in */ int64_t deriv_atom ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.compute)
  // Insert-Code-Here {MPQC.IntegralEvaluator1.compute} (compute method)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.compute)
}

/**
 * Compute a shell singlet of integrals and return as a borrowed
 * sidl array.  deriv_atom must be used for nuclear derivatives if
 * the operator contains nuclear coordinates, otherwise, set to -1
 * and use deriv_ctr.
 * @param shellnum1 Gaussian shell number 1.
 * @param deriv_level Derivative level.
 * @param deriv_atom Atom number for derivative
 * (-1 if using DerivCenter).
 * @return Borrowed sidl array buffer. 
 */
::sidl::array<double>
MPQC::IntegralEvaluator1_impl::compute_array (
  /* in */ int64_t shellnum1,
  /* in */ int32_t deriv_level,
  /* in */ int64_t deriv_atom ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.compute_array)
  // Insert-Code-Here {MPQC.IntegralEvaluator1.compute_array} (compute_array method)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.compute_array)
}


// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1._misc)
// Insert-Code-Here {MPQC.IntegralEvaluator1._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1._misc)

