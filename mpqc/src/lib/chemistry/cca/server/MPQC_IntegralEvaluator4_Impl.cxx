// 
// File:          MPQC_IntegralEvaluator4_Impl.cxx
// Symbol:        MPQC.IntegralEvaluator4-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.IntegralEvaluator4
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 
#include "MPQC_IntegralEvaluator4_Impl.hxx"

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
// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4._includes)

#include <algorithm>

#include <ChemistryIntegralDescrCXX_Eri4IntegralDescr.hxx>
#include <ChemistryIntegralDescrCXX_R12IntegralDescr.hxx>
#include <ChemistryIntegralDescrCXX_R12T1IntegralDescr.hxx>
#include <ChemistryIntegralDescrCXX_R12T2IntegralDescr.hxx>

sc::Ref<sc::GaussianBasisSet>
basis_cca_to_sc( Chemistry::QC::GaussianBasis::MolecularInterface& );

using namespace ChemistryIntegralDescrCXX;

// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4._includes)

// speical constructor, used for data wrapping(required).  Do not put code here unless you really know what you're doing!
MPQC::IntegralEvaluator4_impl::IntegralEvaluator4_impl() : StubBase(
  reinterpret_cast< void*>(::MPQC::IntegralEvaluator4::_wrapObj(
  reinterpret_cast< void*>(this))),false) , _wrapped(true){ 
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4._ctor2)
  // Insert-Code-Here {MPQC.IntegralEvaluator4._ctor2} (ctor2)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4._ctor2)
}

// user defined constructor
void MPQC::IntegralEvaluator4_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4._ctor)

  reorder_ = false;

  IntegralDescrInterface desc = Eri4IntegralDescr::_create();
  descr_to_tbint_type_[desc.get_type()] = sc::TwoBodyInt::eri;
  desc = R12IntegralDescr::_create();
  descr_to_tbint_type_[desc.get_type()] = sc::TwoBodyInt::r12;
  desc = R12T1IntegralDescr::_create();
  descr_to_tbint_type_[desc.get_type()] = sc::TwoBodyInt::r12t1;
  desc = R12T2IntegralDescr::_create();
  descr_to_tbint_type_[desc.get_type()] = sc::TwoBodyInt::r12t2;

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4._ctor)
}

// user defined destructor
void MPQC::IntegralEvaluator4_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4._dtor)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4._dtor)
}

// static class initializer
void MPQC::IntegralEvaluator4_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4._load)
  // guaranteed to be called at most once before any other method in this class
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4._load)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 * Method:  add_evaluator[]
 */
void
MPQC::IntegralEvaluator4_impl::add_evaluator_impl (
  /* in */void* eval,
  /* in */::Chemistry::QC::GaussianBasis::IntegralDescrInterface desc ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.add_evaluator)

  if( desc.get_deriv_lvl() == 0 )
    eval_.add_evaluator(eval,desc);
  else
    deriv_eval_.add_evaluator(eval,desc);

  buffer_size_.update( desc.get_deriv_lvl(), desc.get_n_segment() );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.add_evaluator)
}

/**
 * Method:  add_composite_evaluator[]
 */
void
MPQC::IntegralEvaluator4_impl::add_composite_evaluator_impl (
  /* in */void* eval,
  /* in */::Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface cdesc 
    ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.add_composite_evaluator)
  comp_eval_.add_evaluator(eval,cdesc);

  for( int i=0; i<cdesc.get_n_descr(); ++i ) {
    Chemistry::QC::GaussianBasis::IntegralDescrInterface desc = 
      cdesc.get_descr(i);
    buffer_size_.update( desc.get_deriv_lvl(), desc.get_n_segment() );
    std::pair< std::string, int > p(desc.get_type(),desc.get_deriv_lvl());
    comp_ids_.push_back( p );
  } 
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.add_composite_evaluator)
}

/**
 * Method:  set_basis[]
 */
void
MPQC::IntegralEvaluator4_impl::set_basis_impl (
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs1,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs2,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs3,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs4 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.set_basis)

  bs1_ = basis_cca_to_sc( bs1 );
  if( bs2.isSame(bs1) )
    bs2_ = bs1_;
  else
    bs2_ = basis_cca_to_sc( bs2 );
  if( bs3.isSame(bs2) )
    bs3_ = bs2_;
  else
    bs3_ = basis_cca_to_sc( bs3 );
  if( bs4.isSame(bs3) )
    bs4_ = bs3_;
  else
    bs4_ = basis_cca_to_sc( bs4 );

  buffer_size_.init( 4, bs1_, bs2_, bs3_, bs4_ );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.set_basis)
}

/**
 * Method:  init_reorder[]
 */
void
MPQC::IntegralEvaluator4_impl::init_reorder_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.init_reorder)

  reorder_ = true;
  reorder_engine_.init( 4, bs1_, bs2_, bs3_, bs4_ );

  CompositeIntegralDescrInterface desc = eval_.get_descriptor();
  for( int i=0; i < desc.get_n_descr(); ++i) {
    IntegralDescrInterface idesc = desc.get_descr(i);
    reorder_engine_.add_buffer( eval_.get_buffer( idesc ), idesc );
  }

  CompositeIntegralDescrInterface deriv_desc = deriv_eval_.get_descriptor();
  for( int i=0; i < deriv_desc.get_n_descr(); ++i) {
    IntegralDescrInterface idesc = deriv_desc.get_descr(i);
    reorder_engine_.add_buffer( deriv_eval_.get_buffer( idesc ), idesc );
  }

  CompositeIntegralDescrInterface cdesc = comp_eval_.get_descriptor();
  for( int i=0; i < cdesc.get_n_descr(); ++i ) {
    IntegralDescrInterface idesc = deriv_desc.get_descr(i);
    reorder_engine_.add_buffer( 
      comp_eval_.get_buffer( idesc, descr_to_tbint_type_[idesc.get_type()] ), 
      idesc );
  }

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.init_reorder)
}

/**
 *  Get buffer pointer for given type.
 * @return Buffer pointer. 
 */
void*
MPQC::IntegralEvaluator4_impl::get_buffer_impl (
  /* in */::Chemistry::QC::GaussianBasis::IntegralDescrInterface desc ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.get_buffer)

  double *b;
  if( desc.get_deriv_lvl() == 0 )
    b = eval_.get_buffer( desc );
  else
    b = deriv_eval_.get_buffer( desc );

  if( b == NULL )
    b = comp_eval_.get_buffer( desc, descr_to_tbint_type_[desc.get_type()] );

  return b;

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.get_buffer)
}

/**
 *  Get sidl array buffer for given type.
 * @param desc Integral descriptor.
 * @return Sidl array. 
 */
::sidl::array<double>
MPQC::IntegralEvaluator4_impl::get_array_impl (
  /* in */::Chemistry::QC::GaussianBasis::IntegralDescrInterface desc ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.get_array)

  sidl::array<double> sa;
  double *b;
  if( desc.get_deriv_lvl() == 0 ) {
    b = eval_.get_buffer( desc );
    if( b != NULL )
      return eval_.get_array( desc, buffer_size_.size() );
  }
  else {
    b = deriv_eval_.get_buffer( desc );
    if( b != NULL )
      return deriv_eval_.get_array( desc, buffer_size_.size() );
  }

  if( b == NULL )
    return comp_eval_.get_array( desc, 
                                 descr_to_tbint_type_[desc.get_type()], 
                                 buffer_size_.size() );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.get_array)
}

/**
 * Method:  get_descriptor[]
 */
::Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface
MPQC::IntegralEvaluator4_impl::get_descriptor_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.get_descriptor)

  CompositeIntegralDescrInterface cdesc = 
    ChemistryIntegralDescrCXX::CompositeIntegralDescr::_create();
  CompositeIntegralDescrInterface desc =  eval_.get_descriptor();
  CompositeIntegralDescrInterface deriv_desc = deriv_eval_.get_descriptor();
  CompositeIntegralDescrInterface comp_desc = comp_eval_.get_descriptor();
  for( int i=0; i<desc.get_n_descr(); ++i)
    cdesc.add_descr( desc.get_descr(i) );
  for( int i=0; i<deriv_desc.get_n_descr(); ++i)
    cdesc.add_descr( deriv_desc.get_descr(i) );
  for( int i=0; i<comp_desc.get_n_descr(); ++i)
    cdesc.add_descr( comp_desc.get_descr(i) );

  return cdesc;
  
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.get_descriptor)
}

/**
 *  This should be called when the object is no longer needed.
 * No other members may be called after finalize. 
 */
int32_t
MPQC::IntegralEvaluator4_impl::finalize_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.finalize)
  // Insert-Code-Here {MPQC.IntegralEvaluator4.finalize} (finalize method)
  // 
  // This method has not been implemented
  // 
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "finalize");
    throw ex;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.finalize)
}

/**
 *  Compute a shell quartet of integrals.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @param shellnum3 Gaussian shell number 3.
 * @param shellnum4 Gaussian shell number 4.
 * @param deriv_level Derivative level. 
 */
void
MPQC::IntegralEvaluator4_impl::compute_impl (
  /* in */int64_t shellnum1,
  /* in */int64_t shellnum2,
  /* in */int64_t shellnum3,
  /* in */int64_t shellnum4 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.compute)

  computer_.set_shells( shellnum1, shellnum2, shellnum3, shellnum4 );
  deriv_computer_.set_shells( shellnum1, shellnum2, shellnum3, shellnum4 );
  computer2_.set_shells( shellnum1, shellnum2, shellnum3, shellnum4 );
  eval_.compute( &computer_ );
  deriv_eval_.compute( &deriv_computer_ );
  comp_eval_.compute( &computer2_ );
  if( reorder_ )
    reorder_engine_.do_it( shellnum1, shellnum2, shellnum3, shellnum4 );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.compute)
}

/**
 *  Compute a shell quartet of integrals and return as a borrowed
 * sidl array.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @param shellnum3 Gaussian shell number 3.
 * @param shellnum4 Gaussian shell number 4.
 * @return Borrowed sidl array buffer. 
 */
::sidl::array<double>
MPQC::IntegralEvaluator4_impl::compute_array_impl (
  /* in */const ::std::string& type,
  /* in */int32_t deriv_lvl,
  /* in */int64_t shellnum1,
  /* in */int64_t shellnum2,
  /* in */int64_t shellnum3,
  /* in */int64_t shellnum4 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.compute_array)

  sidl::array<double> array;

  bool is_comp = false;
  std::vector< std::pair< std::string, int > >::iterator 
    comp_id_it = comp_ids_.begin();
  for( ; comp_id_it != comp_ids_.end(); ++comp_id_it ) {
    if( type == (*comp_id_it).first && deriv_lvl == (*comp_id_it).second ){
      is_comp = true;
      computer2_.set_shells( shellnum1, shellnum2, shellnum3, shellnum4 );
      comp_eval_.set_shells( shellnum1, shellnum2, shellnum3, shellnum4 );
      array = comp_eval_.compute_array( &computer2_, type,
                                         deriv_lvl, buffer_size_.size() );
      // currently no comp_eval will need to reorder (cints only!)
      // this fact saves lots of headache
    }
  }

  if( !is_comp ) {
    if( deriv_lvl == 0 ) {
      computer_.set_shells( shellnum1, shellnum2, shellnum3, shellnum4 );
      array = eval_.compute_array( &computer_, type, 
                                   deriv_lvl, buffer_size_.size() );
    }
    else {
      deriv_computer_.set_shells( shellnum1, shellnum2, shellnum3, shellnum4 );
      array = deriv_eval_.compute_array( &deriv_computer_, type,
                                         deriv_lvl, buffer_size_.size() );
    }
    // currently bad for multiple evals because EVERY buffer is reordered
    // every time do_it() is called
    if( reorder_ )
      reorder_engine_.do_it( shellnum1, shellnum2, shellnum3, shellnum4 );
  }

  return array;

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.compute_array)
}

/**
 *  Compute integral bounds.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @param shellnum3 Gaussian shell number 3.
 * @param shellnum4 Gaussian shell number 4. 
 */
double
MPQC::IntegralEvaluator4_impl::compute_bounds_impl (
  /* in */int64_t shellnum1,
  /* in */int64_t shellnum2,
  /* in */int64_t shellnum3,
  /* in */int64_t shellnum4 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.compute_bounds)

  computer_.set_shells( shellnum1, shellnum2, shellnum3, shellnum4 );
  deriv_computer_.set_shells( shellnum1, shellnum2, shellnum3, shellnum4 );
  double bnd =  eval_.compute_bounds( &computer_ );
  double d_bnd = deriv_eval_.compute_bounds( &deriv_computer_ );
  double c_bnd = comp_eval_.compute_bounds( &computer2_ );
  bnd = std::max( bnd, d_bnd );
  return std::max( bnd, c_bnd );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.compute_bounds)
}

/**
 *  Compute array of integral bounds.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @param shellnum3 Gaussian shell number 3.
 * @param shellnum4 Gaussian shell number 4. 
 */
::sidl::array<double>
MPQC::IntegralEvaluator4_impl::compute_bounds_array_impl (
  /* in */int64_t shellnum1,
  /* in */int64_t shellnum2,
  /* in */int64_t shellnum3,
  /* in */int64_t shellnum4 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.compute_bounds_array)

  sidl::SIDLException ex = sidl::SIDLException::_create();
  try {
    ex.setNote("function not implemented yet");
    ex.add(__FILE__, __LINE__,"");
  }
  catch(...) { }
  throw ex;

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.compute_bounds_array)
}


// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4._misc)
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4._misc)

