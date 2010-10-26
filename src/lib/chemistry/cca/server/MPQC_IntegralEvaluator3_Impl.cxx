// 
// File:          MPQC_IntegralEvaluator3_Impl.cxx
// Symbol:        MPQC.IntegralEvaluator3-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.IntegralEvaluator3
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 
#include "MPQC_IntegralEvaluator3_Impl.hxx"

// 
// Includes for all method dependencies.
// 
#ifndef included_Chemistry_QC_GaussianBasis_CompositeDescrInterface_hxx
#include "Chemistry_QC_GaussianBasis_CompositeDescrInterface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_DescrInterface_hxx
#include "Chemistry_QC_GaussianBasis_DescrInterface.hxx"
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
// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3._includes)

#include "basis_cca_to_sc.h"

// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3._includes)

// special constructor, used for data wrapping(required).  Do not put code here unless you really know what you're doing!
MPQC::IntegralEvaluator3_impl::IntegralEvaluator3_impl() : StubBase(
  reinterpret_cast< void*>(::MPQC::IntegralEvaluator3::_wrapObj(
  reinterpret_cast< void*>(this))),false) , _wrapped(true){ 
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3._ctor2)
  // Insert-Code-Here {MPQC.IntegralEvaluator3._ctor2} (ctor2)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3._ctor2)
}

// user defined constructor
void MPQC::IntegralEvaluator3_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3._ctor)
  reorder_ = false;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3._ctor)
}

// user defined destructor
void MPQC::IntegralEvaluator3_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3._dtor)
  // add destruction details here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3._dtor)
}

// static class initializer
void MPQC::IntegralEvaluator3_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3._load)
  // guaranteed to be called at most once before any other method in this class
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3._load)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 * Method:  add_evaluator[]
 */
void
MPQC::IntegralEvaluator3_impl::add_evaluator_impl (
  /* in */void* eval,
  /* in */::Chemistry::QC::GaussianBasis::DescrInterface& desc ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.add_evaluator)

  eval_.add_evaluator(eval,desc);
  buffer_size_.update( desc.get_deriv_lvl(), desc.get_n_segment() );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.add_evaluator)
}

/**
 * Method:  set_basis[]
 */
void
MPQC::IntegralEvaluator3_impl::set_basis_impl (
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface& bs1,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface& bs2,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface& bs3 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.set_basis)

  bs1_ = basis_cca_to_sc( bs1 );
  if( bs2.isSame(bs1) )
    bs2_ = bs1_;
  else
    bs2_ = basis_cca_to_sc( bs2 );
  if( bs3.isSame(bs2) )
    bs3_ = bs2_;
  else
    bs3_ = basis_cca_to_sc( bs3 );

  buffer_size_.init( 3, bs1_, bs2_, bs3_, bs4_ );
 
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.set_basis)
}

/**
 * Method:  init_reorder[]
 */
void
MPQC::IntegralEvaluator3_impl::init_reorder_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.init_reorder)

  reorder_ = true;
  reorder_engine_.init( 3, bs1_, bs2_, bs3_, bs4_ );
  CompositeDescrInterface desc = eval_.get_descriptor();
  for( int i=0; i < desc.get_n_descr(); ++i) {
    DescrInterface idesc = desc.get_descr(i);
    reorder_engine_.add_buffer( eval_.get_buffer( idesc ), idesc );
  }

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.init_reorder)
}

/**
 * Method:  get_descriptor[]
 */
::Chemistry::QC::GaussianBasis::CompositeDescrInterface
MPQC::IntegralEvaluator3_impl::get_descriptor_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.get_descriptor)

  return eval_.get_descriptor();

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.get_descriptor)
}

/**
 *  Get sidl array buffer for given type.
 * @param desc Integral descriptor.
 * @return Sidl array. 
 */
::sidl::array<double>
MPQC::IntegralEvaluator3_impl::get_array_impl (
  /* in */::Chemistry::QC::GaussianBasis::DescrInterface& desc ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.get_array)

  return eval_.get_array( desc, buffer_size_.size() );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.get_array)
}

/**
 *  Returns array of integral bounds.  When multiple integral
 * types are supported within an evaluator, the ordering
 * matches the ordering of descriptors returned by 
 * get_descriptor().
 * @return Integral bounds array. 
 */
::sidl::array<double>
MPQC::IntegralEvaluator3_impl::get_bounds_array_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.get_bounds_array)
  // Insert-Code-Here {MPQC.IntegralEvaluator3.get_bounds_array} (get_bounds_array method)
  // 
  // This method has not been implemented
  // 
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "get_bounds_array");
    throw ex;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.get_bounds_array)
}

/**
 *  This should be called when the object is no longer needed.
 * No other members may be called after finalize. 
 */
int32_t
MPQC::IntegralEvaluator3_impl::finalize_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.finalize)
  // Insert-Code-Here {MPQC.IntegralEvaluator3.finalize} (finalize method)
  // 
  // This method has not been implemented
  // 
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "finalize");
    throw ex;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.finalize)
}

/**
 *  Compute a shell triplet of integrals.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @param shellnum3 Gaussian shell number 3.
 * @param deriv_level Derivative level. 
 */
void
MPQC::IntegralEvaluator3_impl::compute_impl (
  /* in */int64_t shellnum1,
  /* in */int64_t shellnum2,
  /* in */int64_t shellnum3 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.compute)

  computer_.set_shells( shellnum1, shellnum2, shellnum3 );
  eval_.compute( &computer_ );
  if( reorder_ )
    reorder_engine_.do_it( shellnum1, shellnum2, shellnum3, -1 ); // -1 ignored

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.compute)
}

/**
 *  Compute array of integral bounds.  -1 indicates a wild card
 * and the largest possible bound for given non-wild
 * card values is returned.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @param shellnum3 Gaussian shell number 3. 
 */
void
MPQC::IntegralEvaluator3_impl::compute_bounds_impl (
  /* in */int64_t shellnum1,
  /* in */int64_t shellnum2,
  /* in */int64_t shellnum3 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3.compute_bounds)

  sidl::SIDLException ex = sidl::SIDLException::_create();
  try {
    ex.setNote("function not implemented yet");
    ex.add(__FILE__, __LINE__,"");
  }
  catch(...) { }
  throw ex;

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3.compute_bounds)
}


// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3._misc)
// Put miscellaneous code here
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3._misc)

