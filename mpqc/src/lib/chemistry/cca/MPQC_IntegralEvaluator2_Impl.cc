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
#include <iostream>
#include <sstream>
#include <util/class/scexception.h>
#pragma implementation "ccaiter.h"
#include <ccaiter.h>

using namespace std;
using namespace Chemistry::QC::GaussianBasis;

Ref<GaussianBasisSet> basis_cca_to_sc(Molecular&);

// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2._includes)

// user-defined constructor.
void MPQC::IntegralEvaluator2_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2._ctor)
  deriv_level_ = -1;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2._ctor)
}

// user-defined destructor.
void MPQC::IntegralEvaluator2_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2._dtor)
#ifndef INTV3_ORDER
  if( package_ == "intv3") {
    delete temp_buffer_;
    for( int i=0; i<=maxam_; ++i)
      delete [] reorder_[i];
    delete [] reorder_;
  }
#endif
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
 * Method:  set_integral_package[]
 */
void
MPQC::IntegralEvaluator2_impl::set_integral_package (
  /* in */ const ::std::string& label ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2.set_integral_package)
  package_ = label;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2.set_integral_package)
}

/**
 * Initialize the evaluator.
 * @param bs1 Molecular basis on center 1.
 * @param bs2 Molecular basis on center 2.
 * @param label String specifying integral type.
 * @param max_deriv Max derivative to compute.
 * @param storage Available storage in bytes. 
 */
void
MPQC::IntegralEvaluator2_impl::initialize (
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /* in */ const ::std::string& label,
  /* in */ int64_t max_deriv,
  /* in */ int64_t storage ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2.initialize)

  evaluator_label_ = label;

  bs1_ = basis_cca_to_sc( bs1 );
  if( bs1.isSame(bs2) ) 
    bs2_.assign_pointer( bs1_.pointer() );
  else 
    bs2_ = basis_cca_to_sc( bs2 );
  max_nshell2_ = bs1_->max_ncartesian_in_shell() *
    bs2_->max_ncartesian_in_shell();
  maxam_ = max( bs1_->max_angular_momentum(), bs2_->max_angular_momentum() );
  
  std::string is_deriv("");
  if(max_deriv > 0) is_deriv = " derivative";
  std::cout << "  initializing " << package_ << " " << evaluator_label_
            << is_deriv << " integral evaluator\n";
  if ( package_ == "intv3" ) { 
    integral_ = new IntegralV3( bs1_, bs2_ );
  }
#ifdef HAVE_CINTS
  else if ( package_ == "cints" )
    integral_ = new IntegralCints( bs1_, bs2_ );
#endif
  else {
    throw InputError("bad integral package name",
                     __FILE__,__LINE__);
  }

  integral_->set_storage(storage);
  
  int error = 0;
  if(evaluator_label_ == "overlap") 
    switch( max_deriv ) {
    case 0:
      { eval_ = integral_->overlap(); break; }
    case 1:
      { deriv_eval_ = integral_->overlap_deriv(); break; }
    default:
      ++error;
    }

  else if(evaluator_label_ == "kinetic")
    switch( max_deriv ) {
    case 0:
      { eval_ = integral_->kinetic(); break; }
    case 1:
      { deriv_eval_ = integral_->kinetic_deriv(); break; }
    default:
      ++error;
    }
  
  else if(evaluator_label_ == "potential")
    switch( max_deriv ) {
    case 0:
      { eval_ = integral_->nuclear(); break; }
    case 1:
      { deriv_eval_ = integral_->nuclear_deriv(); break; }
    default:
      ++error;
    }
  
  else if(evaluator_label_ == "1eham")
    switch( max_deriv ) {
    case 0:
      { eval_ = integral_->hcore(); break; }
    case 1:
      { deriv_eval_ = integral_->hcore_deriv(); break; }
    default:
      ++error;
    }

  else 
    throw InputError("unrecognized integral type",
                     __FILE__,__LINE__);
  
  if( error ) {
    throw InputError("derivative level not supported",
                     __FILE__,__LINE__);
  }
  
  if( eval_.nonnull() ) { 
    int_type_ = one_body;
    sc_buffer_ = eval_->buffer();
  }
  else if( deriv_eval_.nonnull() ) { 
    int_type_ = one_body_deriv;
    sc_buffer_ = deriv_eval_->buffer();
  }
  else 
    throw ProgrammingError("bad pointer to sc integal evaluator",
                           __FILE__,__LINE__);
  if( !sc_buffer_ )
    throw ProgrammingError("buffer not assigned",
                           __FILE__,__LINE__);
  // get a non-const pointer we can write to
  buf_ = const_cast<double*>( sc_buffer_ );

  if ( package_ == "intv3" ) {
#ifndef INTV3_ORDER
    initialize_reorder_intv3();
#endif
  }


    // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2.initialize)
}

/**
 * Get the buffer pointer
 * @return Buffer pointer 
 */
void*
MPQC::IntegralEvaluator2_impl::get_buffer ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2.get_buffer)
  return const_cast<double*>( sc_buffer_ );
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2.get_buffer)
}

/**
 * Compute a shell doublet of integrals.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @param deriv_level Derivative level. 
 * @param deriv_atom Atom number for derivative (-1 if using DerivCenter).
 * @param deriv_ctr Derivative center descriptor. 
 */
void
MPQC::IntegralEvaluator2_impl::compute (
  /* in */ int64_t shellnum1,
  /* in */ int64_t shellnum2,
  /* in */ int64_t deriv_level,
  /* in */ int64_t deriv_atom,
  /* in */ ::Chemistry::QC::GaussianBasis::DerivCenters deriv_ctr ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2.compute)

  if( int_type_ == one_body )
    eval_->compute_shell( shellnum1, shellnum2 );
  else if( int_type_ == one_body_deriv ) {

    if( deriv_atom == -1) {
      sc_deriv_centers_.clear();
      deriv_eval_->compute_shell( shellnum1, shellnum2, sc_deriv_centers_ );

      deriv_ctr.clear();
      if(sc_deriv_centers_.has_omitted_center())
        deriv_ctr.add_omitted( sc_deriv_centers_.omitted_center(),
                               sc_deriv_centers_.omitted_atom() );
      for( int i=0; i<sc_deriv_centers_.n() ; ++i)
        deriv_ctr.add_center( sc_deriv_centers_.center(i),
                              sc_deriv_centers_.atom(i) );
    }
    else
      deriv_eval_->compute_shell( shellnum1, shellnum2, deriv_atom );

  }
  else 
    throw ProgrammingError("bad evaluator type",
                           __FILE__,__LINE__);

/*
  sc::GaussianShell* s1 = &( bs1_->shell(shellnum1) );
  sc::GaussianShell* s2 = &( bs2_->shell(shellnum2) );
  int nfunc = s1->nfunction() * s2->nfunction();

  if( int_type_ == one_body_deriv ) {
    std::cerr << "cca buffer for shell doublet(" << evaluator_label_ << "):\n";
    std::cerr << "shellnum1: " << shellnum1 << std::endl;
    int nc1 = s1->ncontraction();
    for (int i=0; i<nc1; ++i)
      std::cerr << "am: " << s1->am(i) << std::endl;
    std::cerr << "shellnum2: " << shellnum2 << std::endl;
    int nc2 = s2->ncontraction();
    for (int i=0; i<nc2; ++i)
      std::cerr << "am: " << s2->am(i) << std::endl;
  
    std::cerr << "dx\n";
    for( int i=0; i<nfunc; ++i)
      std::cerr << sc_buffer_[i] << std::endl;
    std::cerr << "dy\n";
    for( int i=nfunc; i<nfunc*2; ++i)
      std::cerr << sc_buffer_[i] << std::endl;
    std::cerr << "dz\n";
    for( int i=nfunc*2; i<nfunc*3; ++i)
      std::cerr << sc_buffer_[i] << std::endl;
  }
*/

#ifndef INTV3_ORDER
  if( package_ == "intv3") reorder_intv3( shellnum1, shellnum2 );
#endif

/*
  if( int_type_ == one_body_deriv ) {
    std::cerr << "buffer for shell doublet (after reorder): " << shellnum1 << " " << shellnum2 << std::endl;
    std::cerr << "dx\n";
    for( int i=0; i<nfunc; ++i)
      std::cerr << sc_buffer_[i] << std::endl;
    std::cerr << "dy\n";
    for( int i=nfunc; i<nfunc*2; ++i)
      std::cerr << sc_buffer_[i] << std::endl;
    std::cerr << "dz\n";
    for( int i=nfunc*2; i<nfunc*3; ++i)
      std::cerr << sc_buffer_[i] << std::endl;
  }

  sc::GaussianShell &s1 = bs1_->shell(shellnum1);
  sc::GaussianShell &s2 = bs2_->shell(shellnum2);
  int nfunc = s1.nfunction() * s2.nfunction();
  cout << "buffer " << shellnum1 << " " << shellnum2 << endl; 
  for( int i=0; i<nfunc; ++i)
    cout << sc_buffer_[i] << endl;
*/

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2.compute)
}

/**
 * Compute a shell doublet of integrals and return as a borrowed
 * sidl array.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @param deriv_level Derivative level.
 * @param deriv_atom Atom number for derivative (-1 if using DerivCenter).
 * @param deriv_ctr Derivative center descriptor.
 * @return Borrowed sidl array buffer. 
 */
::sidl::array<double>
MPQC::IntegralEvaluator2_impl::compute_array (
  /* in */ int64_t shellnum1,
  /* in */ int64_t shellnum2,
  /* in */ int64_t deriv_level,
  /* in */ int64_t deriv_atom,
  /* in */ ::Chemistry::QC::GaussianBasis::DerivCenters deriv_ctr ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2.compute_array)

  compute( shellnum1, shellnum2, deriv_level, deriv_atom, deriv_ctr );

  // create a proxy SIDL array
  int lower[1] = {0};
  int upper[1]; upper[0] = max_nshell2_-1;
  int stride[1] = {1};
  sidl_buffer_.borrow( const_cast<double*>(sc_buffer_), 1, 
                       lower, upper, stride);
  return sidl_buffer_;

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2.compute_array)
}


// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2._misc)

void
MPQC::IntegralEvaluator2_impl::initialize_reorder_intv3() 
{

  if( int_type_ == one_body )
    temp_buffer_ = new double[max_nshell2_];
  else if( int_type_ == one_body_deriv )
    temp_buffer_ = new double[max_nshell2_*3];

  reorder_ = new int*[maxam_+1];
  reorder_[0] = new int[1];
  reorder_[0][0] = 0;

  for( int i=1; i<=maxam_; ++i) {

    sc::CartesianIter *v3iter = integral_->new_cartesian_iter(i);
    MPQC::CartesianIterCCA iter(i);
    MPQC::CartesianIterCCA *ccaiter = &iter;
    ccaiter->start();
    int ncf = ccaiter->n();
    
    reorder_[i] = new int[ncf];
    v3iter->start();
    for( int j=0; j<ncf; ++j) {
      ccaiter->start();
      for( int k=0; k<ncf; ++k) {
        if( v3iter->a() == ccaiter->a() &&
            v3iter->b() == ccaiter->b() &&
            v3iter->c() == ccaiter->c() ) {
          reorder_[i][j] = k;
          k=ncf; //break k loop
        }
        else ccaiter->next();
      }
      v3iter->next();
    }
  }

}


void
MPQC::IntegralEvaluator2_impl::reorder_intv3(int64_t shellnum1,int64_t shellnum2)
{

  sc::GaussianShell* s1 = &( bs1_->shell(shellnum1) );
  sc::GaussianShell* s2 = &( bs2_->shell(shellnum2) );
  int nc1 = s1->ncontraction();
  int nc2 = s2->ncontraction();

  int reorder_needed=0;
  for (int i=0; i<nc1; ++i) {
    if( s1->am(i) == 1) reorder_needed=1;
    else if( s1->am(i) > 1 && s1->is_cartesian(i) ) reorder_needed=1;
  }
  if (!reorder_needed)
    for (int i=0; i<nc2; ++i) {
      if( s2->am(i) == 1) reorder_needed=1;
      else if( s2->am(i) > 1 && s2->is_cartesian(i) ) reorder_needed=1;
    }
  if( !reorder_needed ) return;

  // copy buffer into temp space
  int nfunc = s1->nfunction() * s2->nfunction();
  if( int_type_ == one_body_deriv ) 
    for( int i=0; i<nfunc*3; ++i)
      temp_buffer_[i] = sc_buffer_[i];
  else
    for( int i=0; i<nfunc; ++i)
      temp_buffer_[i] = sc_buffer_[i];

  // a derivative buffer is composed of 3 "doublets"
  int deriv_offset;
  if( int_type_ == one_body )
    reorder_doublet( s1, s2, nc1, nc2, 0 );
  else if( int_type_ == one_body_deriv )
    for(int i=0; i<3; ++i) {
      deriv_offset = i*nfunc;
      reorder_doublet( s1, s2, nc1, nc2, deriv_offset );
    }


}



void
MPQC::IntegralEvaluator2_impl::reorder_doublet( sc::GaussianShell* s1, sc::GaussianShell* s2,
                                                int nc1, int nc2, int deriv_offset )
{

  int index=deriv_offset, con2_offset=0, local2_offset, 
      c1_base, c2_base;

  int temp;
  for( int c2=0; c2<nc2; ++c2 )
    con2_offset += s2->nfunction(c2);

  int s1_is_cart, s2_is_cart, s1_nfunc, s2_nfunc; 

  for( int c1=0; c1<nc1; ++c1 ) {

    c1_base = index;
    s1_is_cart = s1->is_cartesian(c1);
    s1_nfunc = s1->nfunction(c1);

    for( int fc1=0; fc1<s1_nfunc; ++fc1 ) {

      if( s1_is_cart )
        c2_base = c1_base + reorder_[s1->am(c1)][fc1] * con2_offset;
      else
        c2_base = c1_base + fc1 * con2_offset;

      local2_offset = 0;
      for( int c2=0; c2<nc2; ++c2 ) {
        if( c2>0 ) local2_offset += s2->nfunction(c2-1);
        s2_is_cart = s2->is_cartesian(c2);
        s2_nfunc = s2->nfunction(c2);

        if( s2_is_cart )
          for( int fc2=0; fc2<s2_nfunc; ++fc2 ) {
            buf_[ c2_base + local2_offset + reorder_[s2->am(c2)][fc2] ]
              = temp_buffer_[index];
              ++index;
          }
        else
          for( int fc2=0; fc2<s2_nfunc; ++fc2 ) {
            buf_[ c2_base + local2_offset + fc2 ] = temp_buffer_[index];
            ++index;
          }
      }
    }
  }

}

// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2._misc)

