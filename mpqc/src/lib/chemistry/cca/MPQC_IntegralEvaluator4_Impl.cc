// 
// File:          MPQC_IntegralEvaluator4_Impl.cc
// Symbol:        MPQC.IntegralEvaluator4-v0.2
// Symbol Type:   class
// Babel Version: 0.9.8
// Description:   Server-side implementation for MPQC.IntegralEvaluator4
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.9.8
// 
#include "MPQC_IntegralEvaluator4_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4._includes)
#include <iostream>
#include <sstream>
#include <ccaiter.h>

using namespace std;
using namespace Chemistry::QC::GaussianBasis;

Ref<GaussianBasisSet> basis_cca_to_sc(Molecular&);
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4._includes)

// user defined constructor
void MPQC::IntegralEvaluator4_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4._ctor)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4._ctor)
}

// user defined destructor
void MPQC::IntegralEvaluator4_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4._dtor)
  if( package_ == "intv3") delete [] temp_buffer_;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4._dtor)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 * Method:  set_integral_package[]
 */
void
MPQC::IntegralEvaluator4_impl::set_integral_package (
  /*in*/ const ::std::string& label ) 
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
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs3,
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs4,
  /*in*/ const ::std::string& label,
  /*in*/ int64_t max_deriv ) 
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

  std::cout << "  initializing " << package_ << " " << evaluator_label_
            << " integral evaluator\n";
  if ( package_ == "intv3" ) {
    integral_ = new IntegralV3( bs1_ );
    initialize_reorder_intv3();
  }
#ifdef HAVE_CINTS
  else if ( package_ == "cints" )
    integral_ = new IntegralCints( bs1_ );
#endif
  else {
    std::cout << "\nbad integral package name" << std::endl;
    abort();
  }

  int error = 0;

  integral_->set_storage(200000000);

  if(evaluator_label_ == "eri2")
    switch( deriv_level ) {
    case 0:
      { eval_ = integral_->electron_repulsion(); break; }
    case 1:
      { deriv_eval_ = integral_->electron_repulsion_deriv(); break; }
    case 2:
      { deriv_eval_ = integral_->electron_repulsion_deriv(); break; }
    default:
      ++error;
    }

  if(evaluator_label_ == "grt")
    switch( deriv_level ) {
    case 0:
        { eval_ = integral_->grt(); break; }
    default:
      ++error;
    }

  if( error ) {
    std::cerr << "Error in MPQC::integralEvaluator4:\n"
              << "  integral type is either unrecognized or not supported\n";
    abort();
  }

  if( eval_.nonnull() )
    int_type_ = two_body;
  else if( deriv_eval_.nonnull() )
    int_type_ = two_body_deriv;
  else {
    std::cerr << "Error in MPQC::IntegralEvaluator4:\n"
              << "  bad integral evaluator pointer\n";
    abort();
  }

  sc_buffer_ = eval_->buffer();

  max_nshell4_ = bs1_->max_ncartesian_in_shell();
  max_nshell4_ *= bs2_->max_ncartesian_in_shell();
  max_nshell4_ *= bs3_->max_ncartesian_in_shell();
  max_nshell4_ *= bs4_->max_ncartesian_in_shell();

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
 */
void
MPQC::IntegralEvaluator4_impl::compute (
  /*in*/ int64_t shellnum1,
  /*in*/ int64_t shellnum2,
  /*in*/ int64_t shellnum3,
  /*in*/ int64_t shellnum4,
  /*in*/ int64_t deriv_level ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.compute)

  if( int_type_ == two_body ) {
    //eval_->set_redundant(0);
    eval_->compute_shell( (int) shellnum1, (int) shellnum2,
			  (int) shellnum3, (int) shellnum4);
    if( package_ == "intv3") 
      reorder_intv3( shellnum1, shellnum2, shellnum3, shellnum4 );

    sc::GaussianShell &s1 = bs1_->shell(shellnum1);
    sc::GaussianShell &s2 = bs2_->shell(shellnum2);
    sc::GaussianShell &s3 = bs3_->shell(shellnum3);
    sc::GaussianShell &s4 = bs4_->shell(shellnum4);
    int nfunc = s1.nfunction() * s2.nfunction() * s3.nfunction() * s4.nfunction();

    std::cout << "cca buffer:" << bufn_ << std::endl;
    ++bufn_;
    for( int i=0; i<nfunc; ++i)
      std::cout << sc_buffer_[i] << std::endl;
    std::cout << std::endl;

  }
  else {
    std::cout << "Eval4: int_type is " << int_type_ << std::endl
              << " ... aborting\n";
    abort();
  }

  /* deriv wrt what center? interface needs work
  else if( int_type == two_body_deriv ) {
    deriv_eval_ptr_->compute_shell( shellnum1, shellnum2, ??? );
    sc_buffer = deriv_eval_->buffer();
  }
  */

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
 * @return Borrowed sidl array buffer. 
 */
::sidl::array<double>
MPQC::IntegralEvaluator4_impl::compute_array (
  /*in*/ int64_t shellnum1,
  /*in*/ int64_t shellnum2,
  /*in*/ int64_t shellnum3,
  /*in*/ int64_t shellnum4,
  /*in*/ int64_t deriv_level ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.compute_array)

  compute( shellnum1, shellnum2, shellnum3, shellnum4, deriv_level );

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
    ccaiter->start();
    for( int j=0; j<ncf; ++j) {
      v3iter->start();
      for( int k=0; k<ncf; ++k) {
        if( v3iter->a() == ccaiter->a() &&
            v3iter->b() == ccaiter->b() &&
            v3iter->c() == ccaiter->c() ) {
          reorder_[i][j] = k;
          k=ncf; //break k loop
        }
        else v3iter->next();
      }
      ccaiter->next();
    }

    std::cout << "reorder am=" << i << std::endl;
    for( int j=0; j<ncf; ++j)
      std::cout << reorder_[i][j] << " ";
    std::cout << std::endl;

    delete v3iter;
  }

}

void
MPQC::IntegralEvaluator4_impl::reorder_intv3(int64_t shellnum1,
                                             int64_t shellnum2,
                                             int64_t shellnum3,
                                             int64_t shellnum4)
{

  double *buf = const_cast<double*>( sc_buffer_ );

  sc::GaussianShell &s1 = bs1_->shell(shellnum1);
  sc::GaussianShell &s2 = bs2_->shell(shellnum2);
  sc::GaussianShell &s3 = bs3_->shell(shellnum3);
  sc::GaussianShell &s4 = bs4_->shell(shellnum4);
  int nc1 = s1.ncontraction();
  int nc2 = s2.ncontraction();
  int nc3 = s3.ncontraction();
  int nc4 = s4.ncontraction();

  // copy buffer into temp space
  int nfunc = s1.nfunction() * s2.nfunction() * s3.nfunction() * s4.nfunction();
  for( int i=0; i<nfunc; ++i) {
    temp_buffer_[i] = sc_buffer_[i];
  }

  std::cout << "intv3 buffer:" << bufn_ << std::endl;
  for( int i=0; i<nfunc; ++i) 
     std::cout << temp_buffer_[i] << std::endl;
  std::cout << endl;

  int index=0, con2_offset=0, con3_offset=0, con4_offset=0, con_offset,
      local2_offset, local3_offset, local4_offset, 
      c1_base, c2_base, c3_base, c4_base;

  int temp;
  for( int c4=0; c4<nc4; ++c4 )
    con4_offset += s4.nfunction(c4);

  temp = 0;
  con3_offset = con4_offset;
  for( int c3=0; c3<nc3; ++c3 )
    temp += s3.nfunction(c3);
  con3_offset *= temp;

  temp = 0;
  con2_offset = con3_offset;
  for( int c2=0; c2<nc2; ++c2 )
    temp += s2.nfunction(c2);
  con2_offset *= temp;

  for( int c1=0; c1<nc1; ++c1 ) {
    c1_base = index;
    std::cerr << "c1:" << c1 << " am=" << s1.am(c1) << " cartesian="
              << s1.is_cartesian(c1) << std::endl;

    for( int fc1=0; fc1<s1.nfunction(c1); ++fc1 ) {

      if( s1.is_cartesian(c1) )
        c2_base = c1_base + reorder_[s1.am(c1)][fc1] * con2_offset;
      else
        c2_base = c1_base + fc1 * con2_offset;
      std::cerr << "c2_base:" << c2_base << std::endl;

      for( int c2=0; c2<nc2; ++c2 ) {
        std::cerr << "c2:" << c2 << " am=" << s2.am(c2) << std::endl;

        if( c2==0 ) local2_offset = 0;
        else local2_offset += s2.nfunction(c2-1);
        std::cerr << "local2_offset:" << local2_offset << std::endl;

        for( int fc2=0; fc2<s2.nfunction(c2); ++fc2 ) {

          if( s2.is_cartesian(c2) )
            c3_base = c2_base + (local2_offset + reorder_[s2.am(c2)][fc2]) 
                        * con3_offset;
          else
            c3_base = c2_base + (local2_offset + fc2) * con3_offset;
          std::cerr << "c3_base:" << c3_base << std::endl;

          for( int c3=0; c3<nc3; ++c3 ) {
            std::cerr << "c3:" << c3 << " am=" << s3.am(c3) << std::endl;

            if( c3==0 ) local3_offset = 0;
            else local3_offset += s3.nfunction(c3-1);
            std::cerr << "local3_offset:" << local3_offset << std::endl;

            for( int fc3=0; fc3<s3.nfunction(c3); ++fc3 ) {

              if( s3.is_cartesian(c3) )
                c4_base = c3_base + (local3_offset + reorder_[s3.am(c3)][fc3])
                            * con4_offset;
              else
                c4_base = c3_base + (local3_offset + fc3) * con4_offset;
              std::cerr << "c4_base:" << c4_base << std::endl;

              for( int c4=0; c4<nc4; ++c4 ) {
                std::cerr << "c4:" << c4 << " am=" << s4.am(c4) << std::endl;

                if( c4==0 ) local4_offset = 0;
                else local4_offset += s4.nfunction(c4-1);
                std::cerr << "local4_offset:" << local4_offset << std::endl;
          
                for( int fc4=0; fc4<s4.nfunction(c4); ++fc4 ) {

                  if( s4.is_cartesian(c4) ) {
                    buf[index] =
                      temp_buffer_[ c4_base + local4_offset + 
                                   reorder_[s4.am(c4)][fc4] ];
                    std::cerr << "assigning index:" << index
                              << " fc4:" << fc4 
                              << " reorder:" << reorder_[s4.am(c4)][fc4]
                              << std::endl;
                  }
                  else {
                    buf[index] = temp_buffer_[c4_base + local4_offset + fc4];
                    std::cerr << "assigning index:" << index
                              << " fc4:" << fc4 << std::endl;
                  }
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

// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4._misc)

