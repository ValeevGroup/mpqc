// 
// File:          MPQC_IntegralEvaluator1_Impl.cxx
// Symbol:        MPQC.IntegralEvaluator1-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.IntegralEvaluator1
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 
#include "MPQC_IntegralEvaluator1_Impl.hxx"

// 
// Includes for all method dependencies.
// 
#ifndef included_Chemistry_QC_GaussianBasis_CompositeIntegralDescrInterface_hxx
#include "Chemistry_QC_GaussianBasis_CompositeIntegralDescrInterface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralDescrInterface_hxx
#include "Chemistry_QC_GaussianBasis_IntegralDescrInterface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_MolecularInterface_hxx
#include "Chemistry_QC_GaussianBasis_MolecularInterface.hxx"
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
// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1._includes)

sc::Ref<sc::GaussianBasisSet>
basis_cca_to_sc( Chemistry::QC::GaussianBasis::MolecularInterface& );

// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1._includes)

// speical constructor, used for data wrapping(required).  Do not put code here unless you really know what you're doing!
MPQC::IntegralEvaluator1_impl::IntegralEvaluator1_impl() : StubBase(
  reinterpret_cast< void*>(::MPQC::IntegralEvaluator1::_wrapObj(
  reinterpret_cast< void*>(this))),false) , _wrapped(true){ 
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1._ctor2)
  // Insert-Code-Here {MPQC.IntegralEvaluator1._ctor2} (ctor2)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1._ctor2)
}

// user defined constructor
void MPQC::IntegralEvaluator1_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1._ctor)
  reorder_ = false;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1._ctor)
}

// user defined destructor
void MPQC::IntegralEvaluator1_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1._dtor)
  // Insert-Code-Here {MPQC.IntegralEvaluator1._dtor} (destructor)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1._dtor)
}

// static class initializer
void MPQC::IntegralEvaluator1_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1._load)
  // Insert-Code-Here {MPQC.IntegralEvaluator1._load} (class initialization)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1._load)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 * Method:  add_evaluator[]
 */
void
MPQC::IntegralEvaluator1_impl::add_evaluator_impl (
  /* in */void* eval,
  /* in */::Chemistry::QC::GaussianBasis::IntegralDescrInterface desc ) 
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
MPQC::IntegralEvaluator1_impl::set_basis_impl (
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs1 ) 
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
MPQC::IntegralEvaluator1_impl::init_reorder_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.init_reorder)

  reorder_ = true;
  reorder_engine_.init( 4, bs1_, bs2_, bs3_, bs4_ );
  CompositeIntegralDescrInterface desc = eval_.get_descriptor();
  for( int i=0; i < desc.get_n_descr(); ++i) {
    IntegralDescrInterface idesc = desc.get_descr(i);
    reorder_engine_.add_buffer( eval_.get_buffer( idesc ), idesc );
  }

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.init_reorder)
}

/**
 *  Get buffer pointer for given type.
 * @return Buffer pointer. 
 */
void*
MPQC::IntegralEvaluator1_impl::get_buffer_impl (
  /* in */::Chemistry::QC::GaussianBasis::IntegralDescrInterface desc ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.get_buffer)
  
  return eval_.get_buffer( desc );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.get_buffer)
}

/**
 *  Get sidl array buffer for given type.
 * @param desc Integral descriptor.
 * @return Sidl array. 
 */
::sidl::array<double>
MPQC::IntegralEvaluator1_impl::get_array_impl (
  /* in */::Chemistry::QC::GaussianBasis::IntegralDescrInterface desc ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.get_array)

  return eval_.get_array( desc, buffer_size_.size() );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.get_array)
}

/**
 * Method:  get_descriptor[]
 */
::Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface
MPQC::IntegralEvaluator1_impl::get_descriptor_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.get_descriptor)

  return eval_.get_descriptor();

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.get_descriptor)
}

/**
 *  This should be called when the object is no longer needed.
 * No other members may be called after finalize. 
 */
int32_t
MPQC::IntegralEvaluator1_impl::finalize_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.finalize)
  // Insert-Code-Here {MPQC.IntegralEvaluator1.finalize} (finalize method)
  // 
  // This method has not been implemented
  // 
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "finalize");
    throw ex;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.finalize)
}

/**
 *  Compute a shell singlet of integrals.
 * @param shellnum1 Gaussian shell number 1.
 * @param deriv_level Derivative level. 
 */
void
MPQC::IntegralEvaluator1_impl::compute_impl (
  /* in */int64_t shellnum1 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1.compute)

  computer_.set_shells( shellnum1 );
  eval_.compute( &computer_ );
  if( reorder_ )
    reorder_engine_.do_it( shellnum1, -1, -1, -1 ); // -1 ignored

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1.compute)
}

/**
 *  Compute a shell singlet of integrals and return as a borrowed
 * sidl array.
 * @param shellnum1 Gaussian shell number 1.
 * @return Borrowed sidl array buffer. 
 */
::sidl::array<double>
MPQC::IntegralEvaluator1_impl::compute_array_impl (
  /* in */const ::std::string& type,
  /* in */int32_t deriv_lvl,
  /* in */int64_t shellnum1 ) 
{
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
 *  Compute integral bounds.
 * @param shellnum1 Gaussian shell number 1. 
 */
double
MPQC::IntegralEvaluator1_impl::compute_bounds_impl (
  /* in */int64_t shellnum1 ) 
{
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
 *  Compute array of integral bounds.
 * @param shellnum1 Gaussian shell number 1. 
 */
::sidl::array<double>
MPQC::IntegralEvaluator1_impl::compute_bounds_array_impl (
  /* in */int64_t shellnum1 ) 
{
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

