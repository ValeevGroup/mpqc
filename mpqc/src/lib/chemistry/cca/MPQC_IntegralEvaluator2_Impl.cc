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
 * Method:  initialize_by_name[]
 */
void
MPQC::IntegralEvaluator2_impl::initialize_by_name (
  /*in*/ ::Chemistry::Molecule molecule,
  /*in*/ const ::std::string& basis_name,
  /*in*/ const ::std::string& evaluator_label,
  /*in*/ int64_t max_deriv ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2.initialize_by_name)

  //cout << "Initializing MPQC::IntegralEvaluator2 by name\n"
  //     << "  evaluator type is " << evaluator_label << endl;

  molecule_ = molecule;
  evaluator_label_ = evaluator_label;

  // Create an sc::GaussianBasisSet
  ostringstream input;

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

  //cout << "Number of Basis Functions: " << basis_->nbasis() << endl;
  //basis_->print_brief();

  // Initialize the sc::Integral factory
  std::cout << "  initializing " << package_ << " " << evaluator_label_
            << " integral evaluator\n";
  if ( package_ == "intv3" )
    integral_ = new IntegralV3( bs1_, bs2_ );
#ifdef HAVE_CINTS
  else if ( package_ == "cints" )
    integral_ = new IntegralCints( bs1_, bs2_ );
#endif
  else {
    std::cout << "\nbad integral package name" << std::endl;
    abort();
  }
  //integral_->set_basis(*basis_);
  // need set_storage() too?

  // create a sidl buffer
  //cout << "Max ncartesian in shell: " 
  //     << basis_->max_ncartesian_in_shell() << endl;
  max_nshell2_ = bs1_->max_ncartesian_in_shell();
  max_nshell2_ *= max_nshell2_;
  //sidl_buffer_ = sidl::array<double>::create1d(max_nshell2_);

  //cout << "MPQC Evaluator2: initialization by name finished\n";

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2.initialize_by_name)
}

/**
 * Method:  initialize[]
 */
void
MPQC::IntegralEvaluator2_impl::initialize (
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /*in*/ const ::std::string& evaluator_label,
  /*in*/ int64_t max_deriv ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2.initialize)

  evaluator_label_ = evaluator_label;

  bs1_ = basis_cca_to_sc( bs1 );
  if( bs1.isSame(bs2) ) 
    bs2_.assign_pointer( bs1_.pointer() );
  else 
    bs2_ = basis_cca_to_sc( bs2 );

  std::cout << "  initializing " << package_ << " " << evaluator_label_
            << " integral evaluator\n";
  if ( package_ == "intv3" ) 
    integral_ = new IntegralV3( bs1_, bs2_ );
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
  //sidl_buffer_ = sidl::array<double>::create1d(max_nshell2_);


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
  std::cout << "\nbuffer pointer transfered" << endl;

    // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2.initialize)
}

/**
 * Method:  buffer[]
 */
void*
MPQC::IntegralEvaluator2_impl::buffer () 
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2.buffer)
  return const_cast<double*>( sc_buffer_ );
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2.buffer)
}

/**
 * Method:  compute[]
 */
void
MPQC::IntegralEvaluator2_impl::compute (
  /*in*/ int64_t shellnum1,
  /*in*/ int64_t shellnum2,
  /*in*/ int64_t deriv_level ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2.compute)

  if( int_type_ == one_body )
    eval_->compute_shell( (int) shellnum1,  (int) shellnum2 );
  else 
    abort();
 
  /* deriv wrt what center?  interface needs work
  else if( int_type_ == one_body_deriv ) {
    deriv_eval_->compute_shell( (int) shellnum1, (int) shellnum2,??? );
    sc_buffer = deriv_eval_->buffer();
  }
  */

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2.compute)
}

/**
 * Method:  compute_array[]
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
  std::cout << "\nborrowing sidl array" << endl;
  int lower[1] = {0};
  int upper[1]; upper[0] = max_nshell2_-1;
  int stride[1] = {1};
  sidl_buffer_.borrow( const_cast<double*>(sc_buffer_), 1, 
                       lower, upper, stride);

  std::cout << "\nreturning borrowed array" << endl;
  return sidl_buffer_;

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2.compute_array)
}


// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2._misc)
// Put miscellaneous code here
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2._misc)

