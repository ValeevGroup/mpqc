// 
// File:          MPQC_ModelFactory_Impl.cxx
// Symbol:        MPQC.ModelFactory-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.ModelFactory
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 
#include "MPQC_ModelFactory_Impl.hxx"

// 
// Includes for all method dependencies.
// 
#ifndef included_Chemistry_MoleculeInterface_hxx
#include "Chemistry_MoleculeInterface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralEvaluatorFactoryInterface_hxx
#include "Chemistry_QC_GaussianBasis_IntegralEvaluatorFactoryInterface.hxx"
#endif
#ifndef included_Chemistry_QC_ModelInterface_hxx
#include "Chemistry_QC_ModelInterface.hxx"
#endif
#ifndef included_gov_cca_CCAException_hxx
#include "gov_cca_CCAException.hxx"
#endif
#ifndef included_gov_cca_Services_hxx
#include "gov_cca_Services.hxx"
#endif
#ifndef included_sidl_BaseInterface_hxx
#include "sidl_BaseInterface.hxx"
#endif
#ifndef included_sidl_ClassInfo_hxx
#include "sidl_ClassInfo.hxx"
#endif
#ifndef included_sidl_RuntimeException_hxx
#include "sidl_RuntimeException.hxx"
#endif
#ifndef included_sidl_NotImplementedException_hxx
#include "sidl_NotImplementedException.hxx"
#endif
// DO-NOT-DELETE splicer.begin(MPQC.ModelFactory._includes)

#include <iostream>
#include <sstream>
#include <iomanip>
#include <MPQC_Model.hxx>
#include <Chemistry_MoleculeFactoryInterface.hxx>
#include <util/misc/ccaenv.h>

using namespace std;
using namespace sc;

// DO-NOT-DELETE splicer.end(MPQC.ModelFactory._includes)

// speical constructor, used for data wrapping(required).  Do not put code here unless you really know what you're doing!
MPQC::ModelFactory_impl::ModelFactory_impl() : StubBase(reinterpret_cast< 
  void*>(::MPQC::ModelFactory::_wrapObj(reinterpret_cast< void*>(this))),false) 
  , _wrapped(true){ 
  // DO-NOT-DELETE splicer.begin(MPQC.ModelFactory._ctor2)
  // Insert-Code-Here {MPQC.ModelFactory._ctor2} (ctor2)
  // DO-NOT-DELETE splicer.end(MPQC.ModelFactory._ctor2)
}

// user defined constructor
void MPQC::ModelFactory_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.ModelFactory._ctor)

  // ccaffeine has main, could use stovepipe to get command-line vars,
  // but for now use environmental variables and fake argc/argv
  int fake_argc=0;
  char** fake_argv=0;

  // always use MPI message group
  //grp_ = new sc::MessageGrp( &fake_argc, &fake_argv);
  //if (grp_.nonnull())
  //  sc::MessageGrp::set_default_messagegrp(grp_);
  grp_ = sc::MessageGrp::get_default_messagegrp();

  // get thread group
  thread_ = sc::ThreadGrp::initial_threadgrp(fake_argc, fake_argv);
  if( thread_.nonnull() )
    sc::ThreadGrp::set_default_threadgrp(thread_);

  // get memory group
  memory_ = sc::MemoryGrp::initial_memorygrp(fake_argc, fake_argv);
  if (memory_.nonnull())
    sc::MemoryGrp::set_default_memorygrp(memory_);

  std::cout << "  Using " << grp_->class_name()
       << " for message passing (number of nodes = " << grp_->n() 
       << ").\n"; 
  if( thread_.nonnull() )
    std::cout << "  Using " << thread_->class_name()
         << " for threading (number of threads = " << thread_->nthread() 
         << ").\n";
  if( memory_.nonnull() )
    std::cout << "  Using " << memory_->class_name()
         << " for distributed shared memory.\n"; 
  if( grp_.nonnull() && thread_.nonnull() )
    std::cout << "  Total number of processors = " 
         << grp_->n() * thread_->nthread() << endl;

  // DO-NOT-DELETE splicer.end(MPQC.ModelFactory._ctor)
}

// user defined destructor
void MPQC::ModelFactory_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.ModelFactory._dtor)

  services_.unregisterUsesPort("ppf");
  services_.unregisterUsesPort("BasisName");
  services_.unregisterUsesPort("TheoryName");
  services_.unregisterUsesPort("MoleculeFile");
  services_.unregisterUsesPort("MoleculeFactoryInterface");

  services_.removeProvidesPort("ModelFactoryInterface");

  // DO-NOT-DELETE splicer.end(MPQC.ModelFactory._dtor)
}

// static class initializer
void MPQC::ModelFactory_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.ModelFactory._load)
  // Insert-Code-Here {MPQC.ModelFactory._load} (class initialization)
  // DO-NOT-DELETE splicer.end(MPQC.ModelFactory._load)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 *  Set a string for package-specific input.
 * @param input A string giving package-specific input.
 */
void
MPQC::ModelFactory_impl::set_input_string_impl (
  /* in */const ::std::string& input ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.ModelFactory.set_input_string)

  input_string_ = input;

  // DO-NOT-DELETE splicer.end(MPQC.ModelFactory.set_input_string)
}

/**
 *  Set package-specific input filename.
 * @param filename Package-specific input filename.
 */
void
MPQC::ModelFactory_impl::set_input_filename_impl (
  /* in */const ::std::string& input ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.ModelFactory.set_input_filename)

  input_filename_ = input;

  // DO-NOT-DELETE splicer.end(MPQC.ModelFactory.set_input_filename)
}

/**
 *  Set the theory name for Model's created with get_model.
 * @param theory A string giving the name of the theory,
 * for example, B3LYP.
 */
void
MPQC::ModelFactory_impl::set_theory_impl (
  /* in */const ::std::string& theory ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.ModelFactory.set_theory)

  theory_ = theory;

  // DO-NOT-DELETE splicer.end(MPQC.ModelFactory.set_theory)
}

/**
 *  Set the basis set name for Model's created with get_model.
 * @param basis The basis set name to use, for example, aug-cc-pVDZ.
 */
void
MPQC::ModelFactory_impl::set_basis_impl (
  /* in */const ::std::string& basis ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.ModelFactory.set_basis)

  basis_ = basis;

  // DO-NOT-DELETE splicer.end(MPQC.ModelFactory.set_basis)
}

/**
 *  Set the Molecule to use for Model's created with get_model.
 * @param molecule An object of type Molecule.
 */
void
MPQC::ModelFactory_impl::set_molecule_impl (
  /* in */::Chemistry::MoleculeInterface molecule ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.ModelFactory.set_molecule)

  molecule_ = molecule;

  // DO-NOT-DELETE splicer.end(MPQC.ModelFactory.set_molecule)
}

/**
 *  Set the object to use to compute integrals for Model's 
 * created with get_model.
 * @param intfact An object of type 
 * GaussianBasis.IntegralEvaluatorFactory.
 */
void
MPQC::ModelFactory_impl::set_integral_factory_impl (
  /* in */::Chemistry::QC::GaussianBasis::IntegralEvaluatorFactoryInterface 
    intfact ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.ModelFactory.set_integral_factory)
  // Insert-Code-Here {MPQC.ModelFactory.set_integral_factory} (set_integral_factory method)
  // 
  // This method has not been implemented
  // 
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "set_integral_factory");
    throw ex;
  // DO-NOT-DELETE splicer.end(MPQC.ModelFactory.set_integral_factory)
}

/**
 *  Returns a newly created Model.  Before get_model can be called, 
 * set_theory, set_basis, and set_molecule must be called.
 * @return The new Model instance.
 */
::Chemistry::QC::ModelInterface
MPQC::ModelFactory_impl::get_model_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.ModelFactory.get_model)

  int i;

  gov::cca::TypeMap tm = pp_.readConfigurationMap();

  /*
   Currently two possibilities for molecule specification:
     1) we are using python GUI, set_molecule() has already been called
        and !molecule evaluates to FALSE
     2) we are using caffeine proper, we execute the following block to get
        a molecule from the molecule factory
   !!MOLECULE NOT ALLOWED IN KEYVAL INPUT FILE!!
  */

  if( !molecule_ ) { 
    molecule_filename_ = 
      std::string( tm.getString("molecule_filename", 
				"failed molecule_filename fetch") );
    molecule_factory_ = 
      sidl::babel_cast<Chemistry::MoleculeFactoryInterface>(
        services_.getPort("MoleculeFactoryInterface") );
    molecule_factory_.set_molecule_filename(molecule_filename_);
    molecule_ = molecule_factory_.get_molecule();
  }

  std::ostringstream input;

  // form molecule section of keyval
  // we do not allow a molecule in keyval input files
  double conv = molecule_.get_units().convert_to("bohr");
  input
    << "  molecule<Molecule>: (" << std::endl
    << "    symmetry = auto" << std::endl
    << "    unit = bohr" << std::endl
    << "    {n atoms geometry } = {" << std::endl;
  for(i=0;i<molecule_.get_n_atom();++i) {
    input << setprecision(16);
    input << "\t" << i << "\t" << molecule_.get_atomic_number(i)
      << "\t[  " << molecule_.get_cart_coor(i,0)*conv
      << "  " << molecule_.get_cart_coor(i,1)*conv
      << "  " << molecule_.get_cart_coor(i,2)*conv << "  ]\n";
  }
  input << "    }\n";
  input << "  )" << std::endl;

  /*
   Currently three possibilities for obtaining model keyval:
     1) keyval string is supplied via set_keyval()
     2) keyval filename is supplied via parameter port
     3) theory and basis are supplied by built-in parameter port 
        and we can construct a simple model keyval input
  */  
  if( input_filename_.size() ==  0 )
    input_filename_ = tm.getString("keyval_filename","");

  // option 1
  if( input_string_.size() > 0 )
    input << input_string_;

  // option 2
  else if( input_filename_.size() > 0 ) {
    std::cout << "opening keyval file: " << input_filename_ << std::endl;
    ifstream infile(input_filename_.c_str());
    if( !infile ) {
      std::cout << "\nerror: could not open keyval file\n";
      abort();
    }
    int i;
    while( (i=infile.get()) && i!=EOF )
      input << char(i);
  }

  // option 3
  else {

    if( theory_.size() == 0 )
      theory_ = std::string( tm.getString("theory",
      					  "failed theory fetch") );
    if( basis_.size() == 0 )
      basis_  = std::string( tm.getString("basis",
					  "failed basis fetch") );

    if (theory_ == "HF") {
      input << "  model<CLHF>:(" << std::endl;
    }
    else if (theory_ == "B3LYP") {
      input << "  model<CLKS>:(" << std::endl;
      input << "    functional<StdDenFunctional>:(name=B3LYP)" << std::endl;
    }
    else {
      std::cout << "bad theory" << std::endl;
      abort();
    }

    input << "    molecule=$:molecule" << std::endl
          << "    basis<GaussianBasisSet>:(" << std::endl
	  << "      name = \"" << basis_ << "\"" << std::endl
	  << "      molecule = $..:molecule" << std::endl
	  << "    )" << std::endl << "  )" << std::endl;
  }    

  /*
    This factory only runs in the absence of mpqc main,
    where the CCAEnv is otherwise initialized.
    IntegralCCA requires CCAEnv, so if cca integrals are desired,
    a handle to the framework (services_) must be used to initialize
    CCAEnv here.  This fact is transparent to IntegralCCA. 
    Only strange side effect is that all code using CCAEnv now
    appears to the framework as MPQC::ModelFactory.
  */
  use_cca_integrals_ = bool( tm.getBool("use_cca_integrals",0) ); 
  if( use_cca_integrals_ )
    if( !CCAEnv::initialized() ) {
      Ref<CCAFramework> fw = new Ext_CCAFramework(services_);
      CCAEnv::init( fw );
    }

  std::cout << "  model input:" << std::endl << input.str() << std::endl;
  
  MPQC::Model model = MPQC::Model::_create();  
  model.initialize_parsedkeyval("model",input.str());
  model.set_molecule(molecule_);
  
  return model;

  // DO-NOT-DELETE splicer.end(MPQC.ModelFactory.get_model)
}

/**
 *  This should be called when the object is no longer needed. 
 * No other members may be called after finalize. 
 */
int32_t
MPQC::ModelFactory_impl::finalize_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.ModelFactory.finalize)

  if(molecule_factory_._not_nil())
    molecule_factory_.finalize();

  services_.releasePort("BasisName");
  services_.releasePort("TheoryName");
  services_.releasePort("MoleculeFile");
  services_.releasePort("MoleculeFactoryInterface");
  services_.releasePort("CONFIG");
  ppf_.removeParameterPort(tm_, services_);
  services_.releasePort("ppf");

  //services_.unregisterUsesPort("IntegralEvaluatorFactory");

  return 0;

  // DO-NOT-DELETE splicer.end(MPQC.ModelFactory.finalize)
}

/**
 *  Starts up a component presence in the calling framework.
 * @param services the component instance's handle on the framework world.
 * Contracts concerning Svc and setServices:
 * 
 * The component interaction with the CCA framework
 * and Ports begins on the call to setServices by the framework.
 * 
 * This function is called exactly once for each instance created
 * by the framework.
 * 
 * The argument Svc will never be nil/null.
 * 
 * Those uses ports which are automatically connected by the framework
 * (so-called service-ports) may be obtained via getPort during
 * setServices.
 */
void
MPQC::ModelFactory_impl::setServices_impl (
  /* in */::gov::cca::Services services ) 
// throws:
//     ::gov::cca::CCAException
//     ::sidl::RuntimeException
{
  // DO-NOT-DELETE splicer.begin(MPQC.ModelFactory.setServices)

  services_ = services;
  if (services_._is_nil()) return;

  try {
      services_.addProvidesPort( *this, "ModelFactoryInterface", 
				"Chemistry.QC.ModelFactoryInterface", 0);
      services_.registerUsesPort("ppf",
				 "gov.cca.ports.ParameterPortFactory", 0);
      services_.registerUsesPort("BasisName", 
				 "Util.StringProvider", 0);
      services_.registerUsesPort("TheoryName", 
				 "Util.StringProvider", 0);
      services_.registerUsesPort("MoleculeFile", 
				 "Util.StringProvider", 0);
      services_.registerUsesPort("MoleculeFactoryInterface", 
				 "Chemistry.MoleculeFactoryInterface", 0);
  //    services_.registerUsesPort("IntegralEvaluatorFactory",
  //		     "Chemistry.QC.GaussianBasis.IntegralEvaluatorFactory",0);
  }
  catch (gov::cca::CCAException e) {
      std::cout << "Error using services: "
                << e.getNote() << std::endl;
  }

  molecule_ = 0;

  // setup parameters
  try {

    tm_ = services_.createTypeMap();
    if(tm_._is_nil()) {
      std::cerr << "TypeMap is nill\n";
      abort();
    }
    ppf_ = sidl::babel_cast<gov::cca::ports::ParameterPortFactory>(
             services_.getPort("ppf") );
    ppf_.initParameterData(tm_, "CONFIG");
    ppf_.setBatchTitle(tm_,"MPQC ModelFactory Options");
    ppf_.setGroupName(tm_,"Job Specification");
    ppf_.addRequestString(tm_, "theory", "Theory name", 
			  "Theory", "HF");
    ppf_.addRequestString(tm_, "basis", "AO basis name", 
			  "Basis", "STO-3G");
    ppf_.addRequestString(tm_, "molecule_filename", 
		          "Full path to molecule file",
			  "Molecule filename", ""); 
    ppf_.addRequestString(tm_, "keyval_filename",  
		          "Full path to keyval input file",
			  "Keyval filename", "");
    ppf_.addRequestBoolean(tm_,"use_cca_integrals",
                           "Use cca integrals factories",
                           "use_cca_integrals",0);

    //ppf_.addRequestString(tm_, "integral_buffer", "Integral buffer approach",
    //			  "Integral buffer", "opaque");

    ppf_.addParameterPort(tm_, services_);
    services_.releasePort("ppf");

    pp_ = sidl::babel_cast<gov::cca::ports::ParameterPort>( 
            services_.getPort("CONFIG") );
    if (pp_._is_nil()) {
      std::cerr << "getport failed\n";
      abort();
    }

  }
  catch(std::exception& e) {
    std::cerr << "Error in parameter port setup: " << e.what() << std::endl;
  }

  // DO-NOT-DELETE splicer.end(MPQC.ModelFactory.setServices)
}


// DO-NOT-DELETE splicer.begin(MPQC.ModelFactory._misc)
// Insert-Code-Here {MPQC.ModelFactory._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(MPQC.ModelFactory._misc)
