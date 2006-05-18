// 
// File:          MPQC_Chemistry_QC_ModelFactory_Impl.cc
// Symbol:        MPQC.Chemistry_QC_ModelFactory-v0.2
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for MPQC.Chemistry_QC_ModelFactory
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// 
#include "MPQC_Chemistry_QC_ModelFactory_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_ModelFactory._includes)

#include <iostream>
#include <sstream>
#include <iomanip>
#include <MPQC_Chemistry_QC_Model.hh>
#include <Chemistry_MoleculeFactory.hh>

using namespace std;
using namespace sc;
// DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_ModelFactory._includes)

// user-defined constructor.
void MPQC::Chemistry_QC_ModelFactory_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_ModelFactory._ctor)

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

  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_ModelFactory._ctor)
}

// user-defined destructor.
void MPQC::Chemistry_QC_ModelFactory_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_ModelFactory._dtor)
  // add destruction details here
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_ModelFactory._dtor)
}

// static class initializer.
void MPQC::Chemistry_QC_ModelFactory_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_ModelFactory._load)
  // guaranteed to be called at most once before any other method in this class
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_ModelFactory._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Set the theory name for Model's created with get_model.
 * @param theory A string giving the name of the theory, 
 * for example, B3LYP.
 */
void
MPQC::Chemistry_QC_ModelFactory_impl::set_theory (
  /* in */ const ::std::string& theory ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_ModelFactory.set_theory)
  theory_ = theory;
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_ModelFactory.set_theory)
}

/**
 * Set the basis set name for Model's created with get_model.
 * @param basis The basis set name to use, for example, aug-cc-pVDZ.
 */
void
MPQC::Chemistry_QC_ModelFactory_impl::set_basis (
  /* in */ const ::std::string& basis ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_ModelFactory.set_basis)
  basis_ = basis;
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_ModelFactory.set_basis)
}

/**
 * Set the Molecule to use for Model's created with get_model.
 * @param molecule An object of type Molecule.
 */
void
MPQC::Chemistry_QC_ModelFactory_impl::set_molecule (
  /* in */ ::Chemistry::Molecule molecule ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_ModelFactory.set_molecule)
  molecule_ = molecule;
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_ModelFactory.set_molecule)
}

/**
 * Set the object to use to compute integrals for Model's 
 * created with get_model.
 * @param intfact An object of type 
 * GaussianBasis.IntegralEvaluatorFactory.
 */
void
MPQC::Chemistry_QC_ModelFactory_impl::set_integral_factory (
  /* in */ ::Chemistry::QC::GaussianBasis::IntegralEvaluatorFactory intfact ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_ModelFactory.set_integral_factory)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_ModelFactory.set_integral_factory)
}

/**
 * Returns a newly created Model.  Before get_model can be called, 
 * set_theory, set_basis, and set_molecule must be called.
 * @return The new Model instance.
 */
::Chemistry::QC::Model
MPQC::Chemistry_QC_ModelFactory_impl::get_model ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_ModelFactory.get_model)
  
  std::cerr << "getting model\n";

  int i;

  gov::cca::TypeMap tm = pp_.readConfigurationMap();

  std::cerr << "got TypeMap\n";

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
    molecule_factory_ = services_.getPort("MoleculeFactory");
    molecule_factory_.set_molecule_filename(molecule_filename_);
    molecule_ = molecule_factory_.get_molecule();
  }
  std::cerr << "got molecule\n";

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
   Currently two possibilities for obtaining model keyval:
     1) theory and basis are supplied by built-in parameter port 
        and we can construct a simple model keyval input
     2) keyval filename is supplied for us to read from
  */  
  std::string keyval_filename = 
    std::string( tm.getString("keyval_filename", 
			      "failed keyval_filename fetch") );
  if( keyval_filename.size() > 0 ) {
    ifstream infile(keyval_filename.c_str());
    if( !infile ) {
      std::cout << "\nerror: could not open keyval file\n";
      abort();
    }
    int i;
    while( (i=infile.get()) && i!=EOF )
      input << char(i);
  }
  else {

    theory_ = std::string( tm.getString("theory",
					"failed theory fetch") );
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

  // currently needed for integrals stuff
  if( basis_.size() == 0 )
    basis_  = std::string( tm.getString("basis",
					"failed basis fetch") );
  
  std::cout << "  model input:" << std::endl << input.str() << std::endl;

  // hook into integrals component (optional)
  try { eval_factory_ = services_.getPort("IntegralEvaluatorFactory"); }
  catch (...) {}
  if( eval_factory_._not_nil() ) {
    bool use_opaque;
    std::string buffer_str = 
      std::string( tm.getString("integral_buffer",
				"failed integral_buffer fetch") );
    if( buffer_str == "opaque") use_opaque=true;
    else if(buffer_str == "array") use_opaque=false;
    else { std::cerr << "\bunrecognized integral buffer option"; abort(); }
    intcca_ = new IntegralCCA(eval_factory_,use_opaque);
    Integral::set_default_integral( Ref<Integral>(intcca_.pointer()) );
  }
  
  MPQC::Chemistry_QC_Model model = MPQC::Chemistry_QC_Model::_create();  
  model.initialize_parsedkeyval("model",input.str());
  model.set_molecule(molecule_);
  
  return model;

  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_ModelFactory.get_model)
}

/**
 * This can be called when this Model object is no longer needed.  
 * No other members may be called after finalize. 
 */
int32_t
MPQC::Chemistry_QC_ModelFactory_impl::finalize ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_ModelFactory.finalize)
  if (molecule_factory_._not_nil())
      services_.releasePort("MoleculeFactory");
  return 0;
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_ModelFactory.finalize)
}

/**
 * Starts up a component presence in the calling framework.
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
MPQC::Chemistry_QC_ModelFactory_impl::setServices (
  /* in */ ::gov::cca::Services services ) 
throw ( 
  ::gov::cca::CCAException
){
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_ModelFactory.setServices)

  services_ = services;
  if (services_._is_nil()) return;

  try {
      services_.addProvidesPort(self, "ModelFactory", 
				"gov.cca.Port", 0);
      services_.registerUsesPort("ppf",
				 "gov.cca.ports.ParameterPortFactory", 0);
      services_.registerUsesPort("BasisName", 
				 "Util.StringProvider", 0);
      services_.registerUsesPort("TheoryName", 
				 "Util.StringProvider", 0);
      services_.registerUsesPort("MoleculeFile", 
				 "Util.StringProvider", 0);
      services_.registerUsesPort("MoleculeFactory", 
				 "Chemistry.MoleculeFactory", 0);
      services_.registerUsesPort("IntegralEvaluatorFactory",
		     "Chemistry.QC.GaussianBasis.IntegralEvaluatorFactory",0);
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
    ppf_ = services_.getPort("ppf");
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
    ppf_.addRequestString(tm_, "integral_buffer", "Integral buffer approach",
			  "Integral buffer", "opaque");

    ppf_.addParameterPort(tm_, services_);
    services_.releasePort("ppf");

    pp_ = services_.getPort("CONFIG");
    if (pp_._is_nil()) {
      std::cerr << "getport failed\n";
      abort();
    }

    std::cerr << "finished parameter port stuff\n";

    /*
    if (services_._not_nil()) {
      gov::cca::TypeMap tm = services_.createTypeMap();
      services_.registerUsesPort("classicParam",
                                 "gov.cca.ParameterPortFactoryService",tm);
      gov::cca::Port p = services_.getPort("classicParam");
      ccaffeine::ports::PortTranslator portX = p;
      if(portX._not_nil()) {
        classic::gov::cca::Port *cp
          =static_cast<classic::gov::cca::Port*>(portX.getClassicPort());
        if(!cp) {
          std::cout << "Couldn't get classic port" << std::endl;
          return;
        }
        ConfigurableParameterFactory *cpf
          = dynamic_cast<ConfigurableParameterFactory *>(cp);
        ConfigurableParameterPort *pp = setup_parameters(cpf);
        classic::gov::cca::Port *clscp
          = dynamic_cast<classic::gov::cca::Port*>(pp);
        if (!clscp) {
          std::cout << "Couldn't cast to classic::gov::cca::Port"
                    << std::endl;
        }
        void *vp = static_cast<void*>(clscp);
        ccaffeine::ports::PortTranslator provideX
          = ccaffeine::ports::PortTranslator::createFromClassic(vp);

        services_.addProvidesPort(provideX,
                                  "configure", "gov.cca.ports.ParameterPort", tm);

        services_.releasePort("classicParam");
        services_.unregisterUsesPort("classicParam");
      }
    }
    */
  }
  catch(std::exception& e) {
    std::cerr << "Error in parameter port setup: " << e.what() << std::endl;
  }
 
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_ModelFactory.setServices)
}


// DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_ModelFactory._misc)

/*
ConfigurableParameterPort *
MPQC::Chemistry_QC_ModelFactory_impl::setup_parameters(ConfigurableParameterFactory *cpf)
{
  ConfigurableParameterPort * pp = cpf->createConfigurableParameterPort();

  pp->setBatchTitle("PortTranslatorStarter Configuration");
  pp->setGroupName("Model Factory Input");

  theory_param_ = new StringParameter("theory", "Theory name",
                                      "theory", "HF");
  basis_param_  = new StringParameter("basis", "AO basis name",
                                      "basis", "STO-3G");
  molecule_filename_param_ =
                  new StringParameter("molecule_filename", 
                                      "Molecule filename",
                                      "molecule_filename", ""); 
  keyval_filename_param_ =
                  new StringParameter("keyval_filename", 
                                      "Keyval input filename",
                                      "keyval_filename", "");
  integral_buffer_param_ = 
                  new StringParameter("integral_buffer",
                                      "Integral buffer method",
                                      "integral_buffer", "opaque");

  pp->addRequest(theory_param_);
  pp->addRequest(basis_param_);
  pp->addRequest(molecule_filename_param_);
  pp->addRequest(keyval_filename_param_);
  pp->addRequest(integral_buffer_param_);

  return pp;
}
*/

// DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_ModelFactory._misc)
