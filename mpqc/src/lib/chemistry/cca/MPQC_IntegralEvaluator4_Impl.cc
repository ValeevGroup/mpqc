// 
// File:          MPQC_IntegralEvaluator4_Impl.cc
// Symbol:        MPQC.IntegralEvaluator4-v0.2
// Symbol Type:   class
// Babel Version: 0.10.2
// Description:   Server-side implementation for MPQC.IntegralEvaluator4
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.2
// 
#include "MPQC_IntegralEvaluator4_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4._includes)
#include <iostream>
#include <sstream>
#include <util/class/scexception.h>
#include <ccaiter.h>

using namespace std;
using namespace Chemistry::QC::GaussianBasis;

Ref<GaussianBasisSet> basis_cca_to_sc(Molecular&);
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4._includes)

// user-defined constructor.
void MPQC::IntegralEvaluator4_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4._ctor)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4._ctor)
}

// user-defined destructor.
void MPQC::IntegralEvaluator4_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4._dtor)
#ifndef INTV3_ORDER
  if( package_ == "intv3") delete [] temp_buffer_;
#endif
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
 * Method:  set_integral_package[]
 */
void
MPQC::IntegralEvaluator4_impl::set_integral_package (
  /* in */ const ::std::string& label ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.set_integral_package)
  package_ = label;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.set_integral_package)
}

/**
 * Initialize the evaluator.
 * @param bs1 Molecular basis on center 1.
 * @param bs2 Molecular basis on center 2.
 * @param bs3 Molecular basis on center 3.
 * @param bs4 Molecular basis on center 4.
 * @param label String specifying integral type.
 * @param max_deriv Max derivative to compute. 
 */
void
MPQC::IntegralEvaluator4_impl::initialize (
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs3,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs4,
  /* in */ const ::std::string& label,
  /* in */ int64_t max_deriv ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.initialize)
  
  bufn_ = 0;

  evaluator_label_ = label;
  int deriv_level = max_deriv;

  bs1_ = basis_cca_to_sc( bs1 );
  bs2_ = basis_cca_to_sc( bs2 );
  bs3_ = basis_cca_to_sc( bs3 );
  bs4_ = basis_cca_to_sc( bs4 );

  max_nshell4_ = bs1_->max_ncartesian_in_shell();
  max_nshell4_ *= bs2_->max_ncartesian_in_shell();
  max_nshell4_ *= bs3_->max_ncartesian_in_shell();
  max_nshell4_ *= bs4_->max_ncartesian_in_shell();

  std::cout << "  initializing " << package_ << " " << evaluator_label_
            << " integral evaluator\n";
  if ( package_ == "intv3" ) {
    integral_ = new IntegralV3( bs1_ );
#ifdef INTV3_ORDER
    std::cout << "  using intv3 ordering" << std::endl;
#else
    initialize_reorder_intv3();
#endif
  }
#ifdef HAVE_CINTS
  else if ( package_ == "cints" )
    integral_ = new IntegralCints( bs1_ );
#endif
  else
    throw InputError("bad integral package name",
                     __FILE__,__LINE__);

  // a proper solution is required here
  integral_->set_storage(200000000);

  int error = 0;
  if(evaluator_label_ == "eri2")
    switch( deriv_level ) {
    case 0:
      { eval_ = integral_->electron_repulsion(); break; }
    case 1:
      { deriv_eval_ = integral_->electron_repulsion_deriv(); 
        break; 
      }
    default:
      ++error;
    }

  else if(evaluator_label_ == "grt")
    switch( deriv_level ) {
    case 0:
        { eval_ = integral_->grt(); break; }
    default:
      ++error;
    }

  else
    throw InputError("unsupported integral type",
                     __FILE__,__LINE__);

  if( error )
    throw InputError("derivative level not supported",
                     __FILE__,__LINE__);

  if( eval_.nonnull() ) {
    int_type_ = two_body;
    sc_buffer_ = eval_->buffer();
  }
  else if( deriv_eval_.nonnull() ) {
    int_type_ = two_body_deriv;
    sc_buffer_ = deriv_eval_->buffer();
  }
  else
    throw ProgrammingError("bad integral evaluator pointer",
                           __FILE__,__LINE__);
  if( !sc_buffer_ ) 
    throw ProgrammingError("buffer not assigned",
                           __FILE__,__LINE__);

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.initialize)
}

/**
 * Get the buffer pointer.
 * @return Buffer pointer. 
 */
void*
MPQC::IntegralEvaluator4_impl::get_buffer ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.get_buffer)
  return const_cast<double*>( sc_buffer_ );
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.get_buffer)
}

/**
 * Compute a shell quartet of integrals.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @param shellnum3 Gaussian shell number 3.
 * @param shellnum4 Gaussian shell number 4.
 * @param deriv_level Derivative level.
 * @param deriv_ctr Derivative center descriptor. 
 */
void
MPQC::IntegralEvaluator4_impl::compute (
  /* in */ int64_t shellnum1,
  /* in */ int64_t shellnum2,
  /* in */ int64_t shellnum3,
  /* in */ int64_t shellnum4,
  /* in */ int64_t deriv_level,
  /* in */ ::Chemistry::QC::GaussianBasis::DerivCenters deriv_ctr ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.compute)
  
  if( int_type_ == two_body ) {
    eval_->compute_shell( shellnum1, shellnum2,
			  shellnum3, shellnum4);
  }

  else if( int_type_ == two_body_deriv ) {
    sc::DerivCenters dc;

    if(deriv_ctr.has_omitted_center() && deriv_ctr.omitted_center() == 0 )
      dc.add_omitted(0,deriv_ctr.atom(0));
    else
      dc.add_center(0,deriv_ctr.atom(0));

    if(deriv_ctr.has_omitted_center() && deriv_ctr.omitted_center() == 1 )
      dc.add_omitted(1,deriv_ctr.atom(1));
    else
      dc.add_center(1,deriv_ctr.atom(1));

    if(deriv_ctr.has_omitted_center() && deriv_ctr.omitted_center() == 2 )
      dc.add_omitted(2,deriv_ctr.atom(2));
    else
      dc.add_center(2,deriv_ctr.atom(2));

    if(deriv_ctr.has_omitted_center() && deriv_ctr.omitted_center() == 3 )
      dc.add_omitted(3,deriv_ctr.atom(3));
    else
      dc.add_center(3,deriv_ctr.atom(3));

    deriv_eval_->compute_shell( shellnum1, shellnum2, 
                                shellnum3, shellnum4, dc );
  }
  else
    throw ProgrammingError("bad evaluator type",
                           __FILE__,__LINE__);

#ifndef INTV3_ORDER
  if( package_ == "intv3" )
    reorder_intv3( shellnum1, shellnum2, shellnum3, shellnum4 );
#endif

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.compute)
}

/**
 * Compute a shell quartet of integrals and return as a borrowed
 * sidl array.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @param shellnum3 Guassian shell number 3.
 * @param shellnum4 Gaussian shell number 4.
 * @param deriv_level Derivative level.
 * @param deriv_ctr Derivative center descriptor.
 * @return Borrowed sidl array buffer. 
 */
::sidl::array<double>
MPQC::IntegralEvaluator4_impl::compute_array (
  /* in */ int64_t shellnum1,
  /* in */ int64_t shellnum2,
  /* in */ int64_t shellnum3,
  /* in */ int64_t shellnum4,
  /* in */ int64_t deriv_level,
  /* in */ ::Chemistry::QC::GaussianBasis::DerivCenters deriv_ctr ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.compute_array)

  compute( shellnum1, shellnum2, shellnum3, shellnum4, deriv_level, deriv_ctr );

  // this creates a proxy SIDL array
  int lower[1] = {0};
  int upper[1]; upper[0] = max_nshell4_-1;
  int stride[1] = {1};
  sidl_buffer_.borrow( const_cast<double*>(sc_buffer_), 
                       1, lower, upper, stride);

  return sidl_buffer_;

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.compute_array)
}


// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4._misc)

void
MPQC::IntegralEvaluator4_impl::initialize_reorder_intv3()
{
  temp_buffer_ = new double[max_nshell4_];

  int max12 = max( bs1_->max_angular_momentum(), bs2_->max_angular_momentum() );
  int max34 = max( bs3_->max_angular_momentum(), bs4_->max_angular_momentum() );
  int maxam = max( max12, max34 );

  reorder_ = new int*[maxam+1];
  reorder_[0] = new int[1];
  reorder_[0][0] = 0;
  if(maxam==0) return;

  for( int i=1; i<=maxam; ++i) {

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
MPQC::IntegralEvaluator4_impl::reorder_intv3(int64_t shellnum1,
                                             int64_t shellnum2,
                                             int64_t shellnum3,
                                             int64_t shellnum4)
{
  double* buf = const_cast<double*>( sc_buffer_ );

  sc::GaussianShell* s1 = &( bs1_->shell(shellnum1) );
  sc::GaussianShell* s2 = &( bs2_->shell(shellnum2) );
  sc::GaussianShell* s3 = &( bs3_->shell(shellnum3) );
  sc::GaussianShell* s4 = &( bs4_->shell(shellnum4) );
  int nc1 = s1->ncontraction();
  int nc2 = s2->ncontraction();
  int nc3 = s3->ncontraction();
  int nc4 = s4->ncontraction();

  // is this really necessary?
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
  if (!reorder_needed)
    for (int i=0; i<nc3; ++i) {
      if( s3->am(i) == 1) reorder_needed=1;
      else if( s3->am(i) > 1 && s3->is_cartesian(i) ) reorder_needed=1;
    }
  if (!reorder_needed)
    for (int i=0; i<nc4; ++i) {
      if( s4->am(i) == 1) reorder_needed=1;
      else if( s4->am(i) > 1 && s4->is_cartesian(i) ) reorder_needed=1;
    }
  if( !reorder_needed ) return;

  // copy buffer into temp space
  int nfunc = s1->nfunction() * s2->nfunction() * 
                s3->nfunction() * s4->nfunction();
  for( int i=0; i<nfunc; ++i) {
    temp_buffer_[i] = sc_buffer_[i];
  }

  int index=0, con2_offset=0, con3_offset=0, con4_offset=0, con_offset,
      local2_offset, local3_offset, local4_offset, 
      c1_base, c2_base, c3_base, c4_base;

  int temp;
  for( int c4=0; c4<nc4; ++c4 )
    con4_offset += s4->nfunction(c4);

  temp = 0;
  con3_offset = con4_offset;
  for( int c3=0; c3<nc3; ++c3 )
    temp += s3->nfunction(c3);
  con3_offset *= temp;

  temp = 0;
  con2_offset = con3_offset;
  for( int c2=0; c2<nc2; ++c2 )
    temp += s2->nfunction(c2);
  con2_offset *= temp;

  int s1_is_cart, s2_is_cart, s3_is_cart, s4_is_cart,
      s1_nfunc, s2_nfunc, s3_nfunc, s4_nfunc;

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

        for( int fc2=0; fc2<s2_nfunc; ++fc2 ) {

          if( s2_is_cart )
            c3_base = c2_base + (local2_offset + reorder_[s2->am(c2)][fc2]) 
                        * con3_offset;
          else
            c3_base = c2_base + (local2_offset + fc2) * con3_offset;

          local3_offset = 0;
          for( int c3=0; c3<nc3; ++c3 ) {

            if( c3>0 ) local3_offset += s3->nfunction(c3-1);
            s3_is_cart = s3->is_cartesian(c3);
            s3_nfunc = s3->nfunction(c3);

            for( int fc3=0; fc3<s3_nfunc; ++fc3 ) {

              if( s3_is_cart )
                c4_base = c3_base + (local3_offset + reorder_[s3->am(c3)][fc3])
                            * con4_offset;
              else
                c4_base = c3_base + (local3_offset + fc3) * con4_offset;

              local4_offset = 0;
              for( int c4=0; c4<nc4; ++c4 ) {

                if( c4>0 ) local4_offset += s4->nfunction(c4-1);
                s4_is_cart = s4->is_cartesian(c4);
                s4_nfunc = s4->nfunction(c4);
          
                if( s4_is_cart ) 
                  for( int fc4=0; fc4<s4_nfunc; ++fc4 ) {
                    buf[ c4_base + local4_offset + reorder_[s4->am(c4)][fc4] ]
		      = temp_buffer_[index];
                    ++index;
                  }
                else 
                  for( int fc4=0; fc4<s4_nfunc; ++fc4 ) {
                    buf[ c4_base + local4_offset + fc4 ] = temp_buffer_[index];
                    ++index;
                  }
              }
            }
          }
        }
      }
    }
  }

}

void
MPQC::IntegralEvaluator4_impl::reorder_intv3_inline(int64_t shellnum1,
                                                    int64_t shellnum2,
                                                    int64_t shellnum3,
                                                    int64_t shellnum4)
{
  buf_ = const_cast<double*>( sc_buffer_ );

  s1_ = &( bs1_->shell(shellnum1) );
  s2_ = &( bs2_->shell(shellnum2) );
  s3_ = &( bs3_->shell(shellnum3) );
  s4_ = &( bs4_->shell(shellnum4) );
  nc1_ = s1_->ncontraction();
  nc2_ = s2_->ncontraction();
  nc3_ = s3_->ncontraction();
  nc4_ = s4_->ncontraction();

  // is this really necessary?
  int reorder_needed=0;
  for (int i=0; i<nc1_; ++i) {
    if( s1_->am(i) == 1) reorder_needed=1;
    else if( s1_->am(i) > 1 && s1_->is_cartesian(i) ) reorder_needed=1;
  }
  if (!reorder_needed)
    for (int i=0; i<nc2_; ++i) {
      if( s2_->am(i) == 1) reorder_needed=1;
      else if( s2_->am(i) > 1 && s2_->is_cartesian(i) ) reorder_needed=1;
    }
  if (!reorder_needed)
    for (int i=0; i<nc3_; ++i) {
      if( s3_->am(i) == 1) reorder_needed=1;
      else if( s3_->am(i) > 1 && s3_->is_cartesian(i) ) reorder_needed=1;
    }
  if (!reorder_needed)
    for (int i=0; i<nc4_; ++i) {
      if( s4_->am(i) == 1) reorder_needed=1;
      else if( s4_->am(i) > 1 && s4_->is_cartesian(i) ) reorder_needed=1;
    }
  if( !reorder_needed ) return;

  // copy buffer into temp space
  int nfunc = s1_->nfunction() * s2_->nfunction() * 
                s3_->nfunction() * s4_->nfunction();
  for( int i=0; i<nfunc; ++i) {
    temp_buffer_[i] = sc_buffer_[i];
  }

  index_ = con2_offset_ = con3_offset_ = con4_offset_ = 0;

  int temp;
  for( int c4=0; c4<nc4_; ++c4 )
    con4_offset_ += s4_->nfunction(c4);

  temp = 0;
  con3_offset_ = con4_offset_;
  for( int c3=0; c3<nc3_; ++c3 )
    temp += s3_->nfunction(c3);
  con3_offset_ *= temp;

  temp = 0;
  con2_offset_ = con3_offset_;
  for( int c2=0; c2<nc2_; ++c2 )
    temp += s2_->nfunction(c2);
  con2_offset_ *= temp;

  int c3_is_cart, c3_not_cart;

  for( int c1=0; c1<nc1_; ++c1 ) {

    c1_base_ = index_;
    s1_is_cart_ = s1_->is_cartesian(c1);
    s1_nfunc_ = s1_->nfunction(c1);

    for( int fc1=0; fc1<s1_nfunc_; ++fc1 ) {

      if( s1_is_cart_ )
        c2_base_ = c1_base_ + reorder_[s1_->am(c1)][fc1] * con2_offset_;
      else
        c2_base_ = c1_base_ + fc1 * con2_offset_;

      local2_offset_ = 0;
      for( int c2=0; c2<nc2_; ++c2 ) {

        if( c2>0 ) local2_offset_ += s2_->nfunction(c2-1);
        s2_is_cart_ = s2_->is_cartesian(c2);
        s2_nfunc_ = s2_->nfunction(c2);

        for( int fc2=0; fc2<s2_nfunc_; ++fc2 ) {

          if( s2_is_cart_ )
            c3_base_ = c2_base_ + (local2_offset_ + reorder_[s2_->am(c2)][fc2]) 
                        * con3_offset_;
          else
            c3_base_ = c2_base_ + (local2_offset_ + fc2) * con3_offset_;

          local3_offset_ = 0;
          c3_is_cart = s3_->is_cartesian(0);
	  c3_not_cart = !c3_is_cart;
          s3_nfunc_ = s3_->nfunction(0);
          for( int fc3=0; fc3<s3_nfunc_; ++fc3 )
            reorder_c4(0, c3_base_, local3_offset_,
		       c3_is_cart, c3_not_cart, fc3);

          for( int c3=1; c3<nc3_; ++c3 ) {

            local3_offset_ += s3_->nfunction(c3-1);
            c3_is_cart = s3_->is_cartesian(c3);
	    c3_not_cart = !c3_is_cart;
	    s3_nfunc_ = s3_->nfunction(c3);
            for( int fc3=0; fc3<s3_nfunc_; ++fc3 ) {
	      reorder_c4( c3, c3_base_, local3_offset_,
			  c3_is_cart, c3_not_cart, fc3);

            }
          }
        }
      }
    }
  }

}

inline void
MPQC::IntegralEvaluator4_impl::reorder_c4(int c3,
					  int c3_base,
					  int local3_offset,
					  int c3_is_cart,
					  int c3_not_cart,
					  int fc3)
{

  int c4_base, local4_offset, c4_is_cart, c4_not_cart, c4_nfunc;

  c4_base = c3_base 
    + ( c3_is_cart * (local3_offset + reorder_[s3_->am(c3)][fc3]) 
	+ c3_not_cart * (local3_offset + fc3) )
    * con4_offset_;

  local4_offset = 0;
  c4_is_cart = s4_->is_cartesian(0);
  c4_not_cart = !c4_is_cart;
  c4_nfunc = s4_->nfunction(0);
  for( int fc4=0; fc4<c4_nfunc; ++fc4) {
    buf_[ c4_base_ 
	  + c4_is_cart * reorder_[s4_->am(0)][fc4]
	  + c4_not_cart * fc4 ] 
      = temp_buffer_[index_];
    ++index_;
  }

  for( int c4=1; c4<nc4_; ++c4 ) {
    local4_offset += s4_->nfunction(c4-1);
    c4_is_cart = s4_->is_cartesian(c4);
    c4_not_cart = !c4_is_cart;
    c4_nfunc = s4_->nfunction(c4);
    for( int fc4=0; fc4<c4_nfunc; ++fc4 ) {
      buf_[ c4_base + local4_offset 
	    + c4_is_cart * reorder_[s4_->am(c4)][fc4]
            + c4_not_cart * fc4 ] 
	= temp_buffer_[index_];
      ++index_;
    }
  }

}

// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4._misc)

