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
  // add destruction details here
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
 * Method:  initialize_opaque[]
 */
void
MPQC::IntegralEvaluator4_impl::initialize_opaque (
  /*in*/ void* integral ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.initialize_opaque)
  //opaque_ = 1;
  //integral_ = static_cast< Ref<Integral> > integral;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.initialize_opaque)
}

/**
 * Method:  initialize_by_name[]
 */
void
MPQC::IntegralEvaluator4_impl::initialize_by_name (
  /*in*/ ::Chemistry::Molecule molecule,
  /*in*/ const ::std::string& basis_name,
  /*in*/ const ::std::string& evaluator_label,
  /*in*/ int64_t max_deriv ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.initialize_by_name)

  std::cout << "Initializing MPQC::IntegralEvaluator4 by name\n"
            << "  evaluator type is " << evaluator_label << std::endl;

  molecule_ = molecule;
  evaluator_label_ = evaluator_label;

  // Create an sc::GaussianBasisSet
  std::ostringstream input;

  double conv = molecule_.get_units().convert_to("bohr");
  input
    << "  molecule<Molecule>: (\n"
    << "    symmetry = auto\n"
    << "    unit = bohr\n"
    << "    {n atoms geometry } = {\n";
  for(int i=0;i<molecule_.get_n_atom();++i) {
    input
      << "\t" << i << "\t" << molecule_.get_atomic_number(i)
      << "\t[  " << molecule_.get_cart_coor(i,0)*conv
      << "  " << molecule_.get_cart_coor(i,1)*conv
      << "  " << molecule_.get_cart_coor(i,2)*conv << "  ]\n";
  }
  input << "    }\n" << "  )\n";

  input << "  basis<GaussianBasisSet>:(\n"
        << "    name = \"" << basis_name << "\"\n"
        << "    molecule = $:molecule\n"
        << "  )\n";

  sc::Ref<sc::ParsedKeyVal> kv = new sc::ParsedKeyVal();
  kv->parse_string(input.str().c_str());
  sc::Ref<sc::DescribedClass> dc = kv->describedclassvalue("basis");
  bs1_ = dynamic_cast< sc::GaussianBasisSet* >(dc.pointer());

  // Initialize the sc::Integral factory
  // not supporting mixed basis sets yet
  std::cout << "  initializing " << package_ << " " << evaluator_label_
            << " integral evaluator by basis name\n";
  if ( package_ == "intv3" )
    integral_ = new IntegralV3( bs1_ );
#ifdef HAVE_CINTS
  else if ( package_ == "cints" )
    integral_ = new IntegralCints( bs1_ );
#endif
  else {
    std::cout << "\nbad integral package name" << std::endl;
    abort();
  }
  //integral_->set_basis(*basis_);
  // need set_storage() too?

  // create a sidl buffer
  int max_nshell = bs1_->max_ncartesian_in_shell();
  max_nshell4_ = max_nshell;
  max_nshell4_ *= max_nshell;
  max_nshell4_ *= max_nshell;
  max_nshell4_ *= max_nshell;
  //sidl_buffer_ = sidl::array<double>::create1d(max_nshell4_);

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.initialize_by_name)
}

/**
 * Method:  initialize[]
 */
void
MPQC::IntegralEvaluator4_impl::initialize (
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs3,
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs4,
  /*in*/ const ::std::string& evaluator_label,
  /*in*/ int64_t max_deriv ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.initialize)

  evaluator_label_ = evaluator_label;
  int deriv_level = max_deriv;

  bs1_ = basis_cca_to_sc( bs1 );
  bs2_ = basis_cca_to_sc( bs2 );
  bs3_ = basis_cca_to_sc( bs3 );
  bs4_ = basis_cca_to_sc( bs4 );

  //if( !opaque_ ) {
  std::cout << "  initializing " << package_ << " " << evaluator_label_
            << " integral evaluator\n";
  if ( package_ == "intv3" ) 
    integral_ = new IntegralV3( bs1_ );
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
  //sidl_buffer_ = sidl::array<double>::create1d(max_nshell4_);

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.initialize)
}

/**
 * Method:  buffer[]
 */
void*
MPQC::IntegralEvaluator4_impl::buffer () 
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator4.buffer)
  return const_cast<double*>( sc_buffer_ );
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.buffer)
}

/**
 * Method:  compute[]
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

/*
  //std::cout << "Evaluator4: computing shell quartet" << std::endl;

  // Do we have the proper evaluator?
  if( (eval_.null() && deriv_eval_.null()) || (deriv_level_ != deriv_level) ) {

    // No, get the evaluator
    eval_.clear();
    deriv_eval_.clear();
    deriv_level_ = deriv_level;

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
  }
*/

  if( int_type_ == two_body ) {
    //eval_->set_redundant(0);
    eval_->compute_shell( (int) shellnum1, (int) shellnum2,
			  (int) shellnum3, (int) shellnum4);

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

  // copy sc_buffer into a sidl array
  // for now, just copy entire buffer
  // worry about efficiency and reordering later

  //for( int i=0; i<max_nshell4_; ++i) 
  //  sidl_buffer_.set(i, sc_buffer[i]);  

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4.compute)
}

/**
 * Method:  compute_array[]
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
// Put miscellaneous code here
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator4._misc)

