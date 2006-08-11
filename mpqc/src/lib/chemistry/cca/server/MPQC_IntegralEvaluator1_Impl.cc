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

sc::Ref<sc::GaussianBasisSet>
basis_cca_to_sc( Chemistry::QC::GaussianBasis::Molecular& );

// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1._includes)

// user-defined constructor.
void MPQC::IntegralEvaluator1_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1._ctor)
  reorder_ = false;
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
  buffer_size_.update( desc.get_deriv_lvl(), desc.get_n_segment() );

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

  bs1_ = basis_cca_to_sc( bs1 );

  buffer_size_.init( 1, bs1_, bs2_, bs3_, bs4_ );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.set_basis)
}

/**
 * Method:  init_reorder[]
 */
void
MPQC::IntegralEvaluator1_impl::init_reorder ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.init_reorder)

  reorder_ = true;
  reorder_engine_.init( 4, bs1_, bs2_, bs3_, bs4_ );
  CompositeIntegralDescr desc = eval_.get_descriptor();
  for( int i=0; i < desc.get_n_descr(); ++i) {
    IntegralDescr idesc = desc.get_descr(i);
    reorder_engine_.add_buffer( eval_.get_buffer( idesc ), idesc );
  }

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.init_reorder)
}

/**
 * Method:  set_opaque_deriv_centers[]
 */
void
MPQC::IntegralEvaluator1_impl::set_opaque_deriv_centers (
  /* in */ void* dc ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.set_opaque_deriv_centers)
  // Insert-Code-Here {MPQC.IntegralEvaluator1.set_opaque_deriv_centers} (set_opaque_deriv_centers method)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.set_opaque_deriv_centers)
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
  if( reorder_ )
    reorder_engine_.do_it( shellnum1, -1, -1, -1 ); // -1 ignored

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
  /* in */ const ::std::string& type,
  /* in */ int32_t deriv_lvl,
  /* in */ int64_t shellnum1 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.compute_array)

  computer_.set_shells( shellnum1 );
  sidl::array<double> array = 
    eval_.compute_array( &computer_, type, deriv_lvl, buffer_size_.size() );
  if( reorder_ )
    reorder_engine_.do_it( shellnum1, -1, -1 , -1 );

  return array;

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.compute_array)
}

/**
 * Compute integral bounds.
 * @param shellnum1 Gaussian shell number 1. 
 */
double
MPQC::IntegralEvaluator1_impl::compute_bounds (
  /* in */ int64_t shellnum1 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.compute_bounds)

  sidl::SIDLException ex = sidl::SIDLException::_create();
  try {
    ex.setNote("function not implemented yet");
    ex.add(__FILE__, __LINE__,"");
  }
  catch(...) { }
  throw ex;

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.compute_bounds)
}

/**
 * Compute array of integral bounds.
 * @param shellnum1 Gaussian shell number 1. 
 */
::sidl::array<double>
MPQC::IntegralEvaluator1_impl::compute_bounds_array (
  /* in */ int64_t shellnum1 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.compute_bounds_array)

  sidl::SIDLException ex = sidl::SIDLException::_create();
  try {
    ex.setNote("function not implemented yet");
    ex.add(__FILE__, __LINE__,"");
  }
  catch(...) { }
  throw ex;

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.compute_bounds_array)
}


// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1._misc)
// Insert-Code-Here {MPQC.IntegralEvaluator1._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1._misc)

