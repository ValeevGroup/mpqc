// 
// File:          MPQC_IntegralEvaluator2_Impl.cc
// Symbol:        MPQC.IntegralEvaluator2-v0.2
// Symbol Type:   class
// Babel Version: 0.9.8
// Description:   Server-side implementation for MPQC.IntegralEvaluator2
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.9.8
// 
#include "MPQC_IntegralEvaluator2_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2._includes)
#include <iostream>
#include <sstream>

using namespace std;

Ref<GaussianBasisSet> basis_cca_to_sc(Molecular&);

// same as cints, but cints might not be compiled
class CartesianIterCCA : public CartesianIter {
  int *avec, *bvec, *cvec;
  public:
    CartesianIterCCA(int l) : CartesianIter(l) {}
    void start() {
      bfn_=b_=c_=0;
      a_=l_;
    }
    void next() {
      if (c_ < l_ - a_) {
        b_--;
        c_++;
      }
      else {
        a_--;
        c_ = 0;
        b_ = l_ - a_;
      }
      bfn_++;
    }
    operator int() {
      return (a_ >= 0);
    }
};

// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2._includes)

// user defined constructor
void MPQC::IntegralEvaluator2_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2._ctor)
  deriv_level_ = -1;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2._ctor)
}

// user defined destructor
void MPQC::IntegralEvaluator2_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2._dtor)
  // add destruction details here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2._dtor)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 * Method:  set_integral_package[]
 */
void
MPQC::IntegralEvaluator2_impl::set_integral_package (
  /*in*/ const ::std::string& label ) 
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
 */
void
MPQC::IntegralEvaluator2_impl::initialize (
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /*in*/ const ::std::string& label,
  /*in*/ int64_t max_deriv ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2.initialize)

  evaluator_label_ = label;

  bs1_ = basis_cca_to_sc( bs1 );
  if( bs1.isSame(bs2) ) 
    bs2_.assign_pointer( bs1_.pointer() );
  else 
    bs2_ = basis_cca_to_sc( bs2 );
  
  std::cout << "  initializing " << package_ << " " << evaluator_label_
            << " integral evaluator\n";
  if ( package_ == "intv3" ) { 
    integral_ = new IntegralV3( bs1_, bs2_ );
    initialize_reorder_intv3();
  }
#ifdef HAVE_CINTS
  else if ( package_ == "cints" )
    integral_ = new IntegralCints( bs1_, bs2_ );
#endif
  else {
    std::cout << "\nbad integral package name" << std::endl;
    abort();
  }
  
  max_nshell2_ = bs1_->max_ncartesian_in_shell() * 
    bs2_->max_ncartesian_in_shell();

  int error = 0;
  if(evaluator_label_ == "overlap") 
    switch( max_deriv ) {
    case 0:
      { eval_ = integral_->overlap(); break; }
    case 1:
      { deriv_eval_ = integral_->overlap_deriv(); break; }
    case 2:
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
    case 2:
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
    case 2:
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
    case 2:
      { deriv_eval_ = integral_->hcore_deriv(); break; }
    default:
      ++error;
    }
  
  /*else if(evaluator_label_ == "pointcharge")
    switch( deriv_level ) {
    case 0:
    { eval_ = integral_->point_charge(); break; }
    default:
    ++error;
    }
    
    else if(evaluator_label_ == "dipole")
    switch( deriv_level ) {
    case 0:
    { eval_ = integral_->dipole(); break; }
    default:
    ++error;
    }
    
    else if(evaluator_label_ == "quadrupole")
    switch( deriv_level ) {
    case 0:
    { eval_ = integral_->quadrupole(); break; }
    default:
    ++error;
    }*/
  
  if( error ) {
    std::cerr << "Error in MPQC::IntegralEvaluator2:\n"
	      << "  integral type is either unrecognized or not supported\n";
    abort();
  }
  
  if( eval_.nonnull() ) 
    int_type_ = one_body;
  else if( deriv_eval_.nonnull() ) 
    int_type_ = one_body_deriv;
  else {
    std::cerr << "Error in MPQC::IntegralEvaluator2:\n"
	      << "  bad pointer to sc integal evaluator";
    abort();
  }
  
  sc_buffer_ = eval_->buffer();

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
 */
void
MPQC::IntegralEvaluator2_impl::compute (
  /*in*/ int64_t shellnum1,
  /*in*/ int64_t shellnum2,
  /*in*/ int64_t deriv_level ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2.compute)

  if( int_type_ == one_body ) {
    eval_->compute_shell( (int) shellnum1,  (int) shellnum2 );
    if( package_ == "intv3") reorder_intv3( shellnum1, shellnum2 );
  }
  else 
    abort();

  sc::GaussianShell &s1 = bs1_->shell(shellnum1);
  sc::GaussianShell &s2 = bs2_->shell(shellnum2);
  int nfunc = s1.nfunction() * s2.nfunction();
 
  /* deriv wrt what center?  interface needs work
  else if( int_type_ == one_body_deriv ) {
    deriv_eval_->compute_shell( (int) shellnum1, (int) shellnum2,??? );
    sc_buffer = deriv_eval_->buffer();
  }
  */

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2.compute)
}

/**
 * Compute a shell doublet of integrals and return as a borrowed
 * sidl array.
 * @param shellnum1 Gaussian shell number 1.
 * @param shellnum2 Gaussian shell number 2.
 * @param deriv_level Derivative level.
 * @return Borrowed sidl array buffer. 
 */
::sidl::array<double>
MPQC::IntegralEvaluator2_impl::compute_array (
  /*in*/ int64_t shellnum1,
  /*in*/ int64_t shellnum2,
  /*in*/ int64_t deriv_level ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2.compute_array)

  compute( shellnum1, shellnum2, deriv_level );

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

  int maxam = max( bs1_->max_angular_momentum(), bs2_->max_angular_momentum() );

  reorder_ = new int*[maxam+1];
  reorder_[0] = new int[1];
  reorder_[0][0] = 0;
  if(maxam==0) return;

  for( int i=1; i<=maxam; ++i) {

    sc::CartesianIter *v3iter = integral_->new_cartesian_iter(i);
    CartesianIterCCA iter(i);
    CartesianIterCCA *ccaiter = &iter;
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

    delete v3iter;
  }

}

void
MPQC::IntegralEvaluator2_impl::reorder_intv3(int64_t shellnum1,int64_t shellnum2) 
{

  double *buf = const_cast<double*>( sc_buffer_ );

  sc::GaussianShell &s1 = bs1_->shell(shellnum1);
  sc::GaussianShell &s2 = bs2_->shell(shellnum2);
  int nc1 = s1.ncontraction();
  int nc2 = s2.ncontraction();

  // copy buffer into temp space
  int nfunc = s1.nfunction() * s2.nfunction();
  double *temp_buffer = new double[nfunc];
  for( int i=0; i<nfunc; ++i) { 
    temp_buffer[i] = sc_buffer_[i]; 
  }
   
  int index=0, con_offset=0, local_offset, c1_base, c2_base;

  for( int c2=0; c2<nc2; ++c2 )
    con_offset += s2.nfunction(c2);

  for( int c1=0; c1<nc1; ++c1 ) {
    c1_base = index;

    for( int fc1=0; fc1<s1.nfunction(c1); ++fc1 ) {

      for( int c2=0; c2<nc2; ++c2 ) {

        if( c2==0 ) local_offset = 0;
        else local_offset += s2.nfunction(c2-1);

        if( s1.is_cartesian(c1) ) 
          c2_base = c1_base + reorder_[s1.am(c1)][fc1] * con_offset;
        else
          c2_base = c1_base + fc1 * con_offset;

        for( int fc2=0; fc2<s2.nfunction(c2); ++fc2 ) {

          if( s2.is_cartesian(c2) )
            buf[index] = 
              temp_buffer[c2_base + local_offset + reorder_[s2.am(c2)][fc2]];
          else
            buf[index] = temp_buffer[c2_base + local_offset + fc2];
          ++index;

        }
      }
    }
  }

}

// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2._misc)

