// 
// File:          MPQC_IntegralEvaluatorFactory_Impl.cc
// Symbol:        MPQC.IntegralEvaluatorFactory-v0.2
// Symbol Type:   class
// Babel Version: 0.10.2
// Description:   Server-side implementation for MPQC.IntegralEvaluatorFactory
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.2
// 
#include "MPQC_IntegralEvaluatorFactory_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory._includes)
#include <iostream>
#include <chemistry/qc/intv3/intv3.h>
#ifdef HAVE_CINTS
  #include <chemistry/qc/cints/cints.h>
#endif
#include <MPQC_IntegralEvaluator2.hh>
#include <MPQC_IntegralEvaluator3.hh>
#include <MPQC_IntegralEvaluator4.hh>
#include <Chemistry_QC_GaussianBasis_Package.hh>
using namespace std;
using namespace sc;
using namespace Chemistry::QC;
using namespace Chemistry::QC::GaussianBasis;
sc::Ref<sc::GaussianBasisSet> basis_cca_to_sc(GaussianBasis::Molecular&);
string get_Integral_keyval(string);
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory._includes)

// user-defined constructor.
void MPQC::IntegralEvaluatorFactory_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory._ctor)
  // add construction details here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory._ctor)
}

// user-defined destructor.
void MPQC::IntegralEvaluatorFactory_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory._dtor)
  // add destruction details here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory._dtor)
}

// static class initializer.
void MPQC::IntegralEvaluatorFactory_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory._load)
  // guaranteed to be called at most once before any other method in this class
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Starts up a component presence in the calling framework.
 * @param Svc the component instance's handle on the framework world.
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
MPQC::IntegralEvaluatorFactory_impl::setServices (
  /* in */ ::gov::cca::Services services ) 
throw ( 
  ::gov::cca::CCAException
){
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.setServices)
  
  services_ = services;
  if (services_._is_nil()) return;

  try {
      services_.addProvidesPort(self, "IntegralEvaluatorFactory",
                       "Chemistry.QC.GaussianBasis.IntegralEvaluatorFactory", 0);
  }
  catch (gov::cca::CCAException e) {
    std::cout << "Error using services: "
	      << e.getNote() << std::endl;
  }

  //setup parameters
  try {

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
                                  "configure", "ParameterPort", tm);

        services_.releasePort("classicParam");
        services_.unregisterUsesPort("classicParam");
      }
    }
  }
  catch(std::exception& e) {
    std::cout << "Exception caught: " << e.what() << std::endl;
  }
  
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.setServices)
}

/**
 * Set the ObInt evaluator config
 * @param config ObInt evaluator config 
 */
void
MPQC::IntegralEvaluatorFactory_impl::set_obint_config (
  /* in */ ::Chemistry::QC::GaussianBasis::ObIntEvalConfig config ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.set_obint_config)
  ob_config_ = config;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.set_obint_config)
}

/**
 * Set the TbInt evaluator config
 * @param config TbInt evaluator config 
 */
void
MPQC::IntegralEvaluatorFactory_impl::set_tbint_config (
  /* in */ ::Chemistry::QC::GaussianBasis::TbIntEvalConfig config ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.set_tbint_config)
  tb_config_ = config;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.set_tbint_config)
}

/**
 * Set the molecular basis 
 * @param molbasis The molecular basis 
 */
void
MPQC::IntegralEvaluatorFactory_impl::set_molecular (
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular molbasis ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.set_molecular)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.set_molecular)
}

/**
 * Get the molecular basis
 * @return The molecular basis 
 */
::Chemistry::QC::GaussianBasis::Molecular
MPQC::IntegralEvaluatorFactory_impl::get_molecular ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.get_molecular)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.get_molecular)
}

/**
 * Set the molecule
 * @param The molecule 
 */
void
MPQC::IntegralEvaluatorFactory_impl::set_molecule (
  /* in */ ::Chemistry::Molecule mol ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.set_molecule)
  molecule_ = mol;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.set_molecule)
}

/**
 * Get the molecule
 * @return The molecule 
 */
::Chemistry::Molecule
MPQC::IntegralEvaluatorFactory_impl::get_molecule ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.get_molecule)
  return molecule_;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.get_molecule)
}

/**
 * Set available storage
 * @param storage Available storage in bytes 
 */
void
MPQC::IntegralEvaluatorFactory_impl::set_storage (
  /* in */ int64_t storage ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.set_storage)
  storage_ = storage;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.set_storage)
}

/**
 * Get a 1-body 1-center integral evaluator
 * @param type ObIntEvalType specifying eval type
 * @param max_deriv Maximum derivative that will be computed
 * @param bs1 Molecular basis set on center 1
 * @return 1-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator1
MPQC::IntegralEvaluatorFactory_impl::get_obint_evaluator1 (
  /* in */ ::Chemistry::QC::GaussianBasis::ObIntEvalType type,
  /* in */ int32_t max_deriv,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::DerivCenters deriv_ctr ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.get_obint_evaluator1)
  // Insert-Code-Here {MPQC.IntegralEvaluatorFactory.get_obint_evaluator1} (get_obint_evaluator1 method)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.get_obint_evaluator1)
}

/**
 * Get a 1-body 2-center integral evaluator
 * @param type ObIntEvalType specifying eval type
 * @param max_deriv Maximum derivative that will be computed
 * @param bs1 Molecular basis set on center 1
 * @param bs2 Molecular basis set on center 2
 * @return 2-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator2
MPQC::IntegralEvaluatorFactory_impl::get_obint_evaluator2 (
  /* in */ ::Chemistry::QC::GaussianBasis::ObIntEvalType type,
  /* in */ int32_t max_deriv,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /* in */ ::Chemistry::QC::GaussianBasis::DerivCenters deriv_ctr ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.get_obint_evaluator2)

  // create the evaluator
  MPQC::IntegralEvaluator2 eval = MPQC::IntegralEvaluator2::_create();

  // determine proper integrals package
  bool pkg_set = false;
  for( int i=0; i<ob_config_.get_n_pkg_config(); ++i)
    if( ob_config_.get_pkg_config_type(i) == type ) {
      eval.set_integral_package( ob_config_.get_pkg_config_pkg(i) );
      pkg_set = true;
    }
  if( !pkg_set ) {
    package_ =  package_param_->getValueString();
    if( package_ == "intv3" )
      eval.set_integral_package( Package_INTV3 );
    else if( package_ == "cints" )
      eval.set_integral_package( Package_CINTS );
    else {
      Package pkg = ob_config_.get_default_pkg();
      if( pkg == Package_INTV3  || pkg == Package_CINTS )
        eval.set_integral_package( pkg );
      else {
	std::cout << "Using Intv3 as last resort\n";
        eval.set_integral_package( Package_INTV3 );
      }
    }
  }

  // initialize
  eval.obint_initialize( bs1, bs2, type, 
                         max_deriv, storage_, deriv_ctr );
  return eval;

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.get_obint_evaluator2)
}

/**
 * Get a 2-body 2-center integral evaluator
 * @param type TbIntEvalType specifying eval type
 * @param max_deriv Maximum derivative that will be computed
 * @param bs1 Molecular basis set on center 1
 * @param bs2 Molecular basis set on center 2
 * @return 2-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator2
MPQC::IntegralEvaluatorFactory_impl::get_tbint_evaluator2 (
  /* in */ ::Chemistry::QC::GaussianBasis::TbIntEvalType type,
  /* in */ int32_t max_deriv,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /* in */ ::Chemistry::QC::GaussianBasis::DerivCenters deriv_ctr ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.get_tbint_evaluator2)
  // Insert-Code-Here {MPQC.IntegralEvaluatorFactory.get_tbint_evaluator2} (get_tbint_evaluator2 method)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.get_tbint_evaluator2)
}

/**
 * Get a 2-body 3-center integral evaluator
 * @param type TbIntEvalType specifying eval type
 * @param max_deriv Maximum derivative that will be computed
 * @param bs1 Molecular basis set on center 1
 * @param bs2 Molecular basis set on center 2
 * @param bs3 Molecular basis set on center 3
 * @return 3-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator3
MPQC::IntegralEvaluatorFactory_impl::get_tbint_evaluator3 (
  /* in */ ::Chemistry::QC::GaussianBasis::TbIntEvalType type,
  /* in */ int32_t max_deriv,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs3,
  /* in */ ::Chemistry::QC::GaussianBasis::DerivCenters deriv_ctr ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.get_tbint_evaluator3)

  // create the evaluator
  MPQC::IntegralEvaluator3 eval = MPQC::IntegralEvaluator3::_create();

  // determine proper integrals package
  bool pkg_set = false;
  for( int i=0; i<tb_config_.get_n_pkg_config(); ++i)
    if( tb_config_.get_pkg_config_type(i) == type ) {
      eval.set_integral_package( tb_config_.get_pkg_config_pkg(i) );
      pkg_set = true;
    }
  if( !pkg_set ) {
    package_ =  package_param_->getValueString();
    if( package_ == "intv3" )
      eval.set_integral_package( Package_INTV3 );
    else if( package_ == "cints" )
      eval.set_integral_package( Package_CINTS );
    else {
      Package pkg = tb_config_.get_default_pkg();
      if( pkg == Package_INTV3  || pkg == Package_CINTS )
        eval.set_integral_package( pkg );
      else {
	std::cout << "Using Intv3 as last resort\n";
        eval.set_integral_package( Package_INTV3 );
      }
    }
  }

  // initialize
  eval.tbint_initialize( bs1, bs2, bs3, type, 
                         max_deriv, storage_, deriv_ctr );
  return eval;

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.get_tbint_evaluator3)
}

/**
 * Get a 2-body 4-center integral evaluator
 * @param type TbIntEvalType specifiying eval type
 * @param max_deriv Maximum derivative that will be computed
 * @param bs1 Molecular basis set on center 1
 * @param bs2 Molecular basis set on center 2
 * @param bs3 Molecular basis set on center 3
 * @param bs4 Molecular basis set on center 4
 * @return 4-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator4
MPQC::IntegralEvaluatorFactory_impl::get_tbint_evaluator4 (
  /* in */ ::Chemistry::QC::GaussianBasis::TbIntEvalType type,
  /* in */ int32_t max_deriv,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs3,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs4,
  /* in */ ::Chemistry::QC::GaussianBasis::DerivCenters deriv_ctr ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.get_tbint_evaluator4)

  // create the evaluator
  MPQC::IntegralEvaluator4 eval = MPQC::IntegralEvaluator4::_create();

  // determine proper integrals package
  bool pkg_set = false;
  for( int i=0; i<tb_config_.get_n_pkg_config(); ++i)
    if( tb_config_.get_pkg_config_type(i) == type ) {
      eval.set_integral_package( tb_config_.get_pkg_config_pkg(i) );
      pkg_set = true;
    }
  if( !pkg_set ) {
    package_ =  package_param_->getValueString();
    if( package_ == "intv3" )
      eval.set_integral_package( Package_INTV3 );
    else if( package_ == "cints" )
      eval.set_integral_package( Package_CINTS );
    else {
      Package pkg = tb_config_.get_default_pkg();
      if( pkg == Package_INTV3  || pkg == Package_CINTS )
        eval.set_integral_package( pkg );
      else {
	std::cout << "Using Intv3 as last resort\n";
        eval.set_integral_package( Package_INTV3 );
      }
    }
  }

  // initialize
  eval.tbint_initialize( bs1, bs2, bs3, bs4, type, 
                         max_deriv, storage_, deriv_ctr );
  return eval;

  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.get_tbint_evaluator4)
}


// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory._misc)

ConfigurableParameterPort *
MPQC::IntegralEvaluatorFactory_impl::setup_parameters(
                                              ConfigurableParameterFactory *cpf)
{
  ConfigurableParameterPort * pp = cpf->createConfigurableParameterPort();

  pp->setBatchTitle("PortTranslatorStarter Configuration");
  pp->setGroupName("Model Factory Input");

  package_param_ = new StringParameter("package", "Integral package",
                                      "package", "novalue");
  pp->addRequest(package_param_);

  return pp;
}

// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory._misc)

