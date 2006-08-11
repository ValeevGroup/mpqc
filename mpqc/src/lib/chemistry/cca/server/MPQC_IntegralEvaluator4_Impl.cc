// 
// File:          MPQC_IntegralEvaluator4_Impl.cc
// Symbol:        MPQC.IntegralEvaluator4-v0.2
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for MPQC.IntegralEvaluator4
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// 
#include "MPQC_IntegralEvaluator4_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4._includes)

#include <algorithm>

#include <Chemistry_Eri4IntegralDescr.hh>
#include <Chemistry_R12IntegralDescr.hh>
#include <Chemistry_R12T1IntegralDescr.hh>
#include <Chemistry_R12T2IntegralDescr.hh>

sc::Ref<sc::GaussianBasisSet>
basis_cca_to_sc( Chemistry::QC::GaussianBasis::Molecular& );

// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4._includes)

// user-defined constructor.
void MPQC::IntegralEvaluator4_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4._ctor)

  reorder_ = false;
  compute_fast_ = false;

  IntegralDescr desc = Chemistry::Eri4IntegralDescr::_create();
  descr_to_tbint_type_[desc.get_type()] = sc::TwoBodyInt::eri;
  desc = Chemistry::R12IntegralDescr::_create();
  descr_to_tbint_type_[desc.get_type()] = sc::TwoBodyInt::r12;
  desc = Chemistry::R12T1IntegralDescr::_create();
  descr_to_tbint_type_[desc.get_type()] = sc::TwoBodyInt::r12t1;
  desc = Chemistry::R12T2IntegralDescr::_create();
  descr_to_tbint_type_[desc.get_type()] = sc::TwoBodyInt::r12t2;

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4._ctor)
}

// user-defined destructor.
void MPQC::IntegralEvaluator4_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4._dtor)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4._dtor)
}

// static class initializer.
void MPQC::IntegralEvaluator4_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4._load)
  // guaranteed to be called at most once before any other method in this class
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Method:  add_evaluator[]
 */
void
MPQC::IntegralEvaluator4_impl::add_evaluator (
  /* in */ void* eval,
  /* in */ ::Chemistry::QC::GaussianBasis::IntegralDescr desc ) 
throw () 
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
MPQC::IntegralEvaluator4_impl::add_composite_evaluator (
  /* in */ void* eval,
  /* in */ ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr cdesc ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.add_composite_evaluator)
  comp_eval_.add_evaluator(eval,cdesc);

  for( int i=0; i<cdesc.get_n_descr(); ++i ) {
    Chemistry::QC::GaussianBasis::IntegralDescr desc = cdesc.get_descr(i);
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
MPQC::IntegralEvaluator4_impl::set_basis (
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs3,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs4 ) 
throw () 
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
MPQC::IntegralEvaluator4_impl::init_reorder ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.init_reorder)

  reorder_ = true;
  reorder_engine_.init( 4, bs1_, bs2_, bs3_, bs4_ );

  CompositeIntegralDescr desc = eval_.get_descriptor();
  for( int i=0; i < desc.get_n_descr(); ++i) {
    IntegralDescr idesc = desc.get_descr(i);
    reorder_engine_.add_buffer( eval_.get_buffer( idesc ), idesc );
  }

  CompositeIntegralDescr deriv_desc = deriv_eval_.get_descriptor();
  for( int i=0; i < deriv_desc.get_n_descr(); ++i) {
    IntegralDescr idesc = deriv_desc.get_descr(i);
    reorder_engine_.add_buffer( deriv_eval_.get_buffer( idesc ), idesc );
  }

  CompositeIntegralDescr cdesc = comp_eval_.get_descriptor();
  for( int i=0; i < cdesc.get_n_descr(); ++i ) {
    IntegralDescr idesc = deriv_desc.get_descr(i);
    reorder_engine_.add_buffer( 
      comp_eval_.get_buffer( idesc, descr_to_tbint_type_[idesc.get_type()] ), 
      idesc );
  }

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.init_reorder)
}

/**
 * Method:  set_opaque_deriv_centers[]
 */
void
MPQC::IntegralEvaluator4_impl::set_opaque_deriv_centers (
  /* in */ void* dc ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.set_opaque_deriv_centers)

  sc_dc_ = static_cast<sc::DerivCenters*>( dc );
  compute_fast_ = true;

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.set_opaque_deriv_centers)
}

/**
 * Get buffer pointer for given type.
 * @return Buffer pointer. 
 */
void*
MPQC::IntegralEvaluator4_impl::get_buffer (
  /* in */ ::Chemistry::QC::GaussianBasis::IntegralDescr desc ) 
throw ( 
  ::sidl::BaseException
){
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
 * Method:  get_descriptor[]
 */
::Chemistry::QC::GaussianBasis::CompositeIntegralDescr
MPQC::IntegralEvaluator4_impl::get_descriptor ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.get_descriptor)

  CompositeIntegralDescr cdesc = Chemistry::CompositeIntegralDescr::_create();
  CompositeIntegralDescr desc =  eval_.get_descriptor();
  CompositeIntegralDescr deriv_desc = deriv_eval_.get_descriptor();
  CompositeIntegralDescr comp_desc = comp_eval_.get_descriptor();
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
 * Compute a shell quartet of integrals.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @param shellnum3 Gaussian shell number 3.
 * @param shellnum4 Gaussian shell number 4.
 * @param deriv_level Derivative level. 
 */
void
MPQC::IntegralEvaluator4_impl::compute (
  /* in */ int64_t shellnum1,
  /* in */ int64_t shellnum2,
  /* in */ int64_t shellnum3,
  /* in */ int64_t shellnum4 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.compute)

  if( compute_fast_ ) {
    deriv_computer_.set_shells( shellnum1, shellnum2, shellnum3, shellnum4 );
    deriv_eval_.compute_fast( &deriv_computer_, sc_dc_ );
  }
  else {
    computer_.set_shells( shellnum1, shellnum2, shellnum3, shellnum4 );
    deriv_computer_.set_shells( shellnum1, shellnum2, shellnum3, shellnum4 );
    computer2_.set_shells( shellnum1, shellnum2, shellnum3, shellnum4 );
    eval_.compute( &computer_ );
    deriv_eval_.compute( &deriv_computer_ );
    comp_eval_.compute( &computer2_ );
  }
  if( reorder_ )
    reorder_engine_.do_it( shellnum1, shellnum2, shellnum3, shellnum4 );

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.compute)
}

/**
 * Compute a shell quartet of integrals and return as a borrowed
 * sidl array.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @param shellnum3 Gaussian shell number 3.
 * @param shellnum4 Gaussian shell number 4.
 * @return Borrowed sidl array buffer. 
 */
::sidl::array<double>
MPQC::IntegralEvaluator4_impl::compute_array (
  /* in */ const ::std::string& type,
  /* in */ int32_t deriv_lvl,
  /* in */ int64_t shellnum1,
  /* in */ int64_t shellnum2,
  /* in */ int64_t shellnum3,
  /* in */ int64_t shellnum4 ) 
throw ( 
  ::sidl::BaseException
){
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
 * Compute integral bounds.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @param shellnum3 Gaussian shell number 3.
 * @param shellnum4 Gaussian shell number 4. 
 */
double
MPQC::IntegralEvaluator4_impl::compute_bounds (
  /* in */ int64_t shellnum1,
  /* in */ int64_t shellnum2,
  /* in */ int64_t shellnum3,
  /* in */ int64_t shellnum4 ) 
throw ( 
  ::sidl::BaseException
){
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
 * Compute array of integral bounds.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @param shellnum3 Gaussian shell number 3.
 * @param shellnum4 Gaussian shell number 4. 
 */
::sidl::array<double>
MPQC::IntegralEvaluator4_impl::compute_bounds_array (
  /* in */ int64_t shellnum1,
  /* in */ int64_t shellnum2,
  /* in */ int64_t shellnum3,
  /* in */ int64_t shellnum4 ) 
throw ( 
  ::sidl::BaseException
){
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

