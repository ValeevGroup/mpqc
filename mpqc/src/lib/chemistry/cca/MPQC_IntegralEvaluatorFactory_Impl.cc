// 
// File:          MPQC_IntegralEvaluatorFactory_Impl.cc
// Symbol:        MPQC.IntegralEvaluatorFactory-v0.2
// Symbol Type:   class
// Babel Version: 0.9.8
// Description:   Server-side implementation for MPQC.IntegralEvaluatorFactory
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.9.8
// 
#include "MPQC_IntegralEvaluatorFactory_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory._includes)
using namespace std;
string get_Integral_keyval(string);
#include <iostream>
#include "MPQC_IntegralEvaluator2.hh"
#include "MPQC_IntegralEvaluator3.hh"
#include "MPQC_IntegralEvaluator4.hh"
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory._includes)

// user defined constructor
void MPQC::IntegralEvaluatorFactory_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory._ctor)
  // add construction details here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory._ctor)
}

// user defined destructor
void MPQC::IntegralEvaluatorFactory_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory._dtor)
  // add destruction details here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory._dtor)
}

// user defined static methods: (none)

// user defined non-static methods:
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
  /*in*/ ::gov::cca::Services services ) 
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
 * Method:  set_basis_name[]
 */
void
MPQC::IntegralEvaluatorFactory_impl::set_basis_name (
  /*in*/ const ::std::string& basis_name ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.set_basis_name)
  basis_name_ = basis_name;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.set_basis_name)
}

/**
 * Method:  get_basis_name[]
 */
::std::string
MPQC::IntegralEvaluatorFactory_impl::get_basis_name () 
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.get_basis_name)
  return basis_name_;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.get_basis_name)
}

/**
 * Method:  set_molecular[]
 */
void
MPQC::IntegralEvaluatorFactory_impl::set_molecular (
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular molbasis ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.set_molecular)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.set_molecular)
}

/**
 * Method:  get_molecular[]
 */
::Chemistry::QC::GaussianBasis::Molecular
MPQC::IntegralEvaluatorFactory_impl::get_molecular (
  /*in*/ int64_t center ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.get_molecular)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.get_molecular)
}

/**
 * Method:  set_molecule[]
 */
void
MPQC::IntegralEvaluatorFactory_impl::set_molecule (
  /*in*/ ::Chemistry::Molecule mol ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.set_molecule)
  molecule_ = mol;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.set_molecule)
}

/**
 * Method:  get_molecule[]
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
 * Method:  get_integral_evaluator2[]
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator2
MPQC::IntegralEvaluatorFactory_impl::get_integral_evaluator2 (
  /*in*/ const ::std::string& label,
  /*in*/ int64_t max_deriv,
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs2 ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.get_integral_evaluator2)
  MPQC::IntegralEvaluator2 eval = MPQC::IntegralEvaluator2::_create();
  string package( package_param_->getValueString() );
  eval.set_integral_package( package );
  //eval.initialize_by_name( molecule_, basis_name_, label, max_deriv);
  eval.initialize( bs1, bs2, label, max_deriv );
  return eval;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.get_integral_evaluator2)
}

/**
 * Method:  get_integral_evaluator3[]
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator3
MPQC::IntegralEvaluatorFactory_impl::get_integral_evaluator3 (
  /*in*/ const ::std::string& label,
  /*in*/ int64_t max_deriv,
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs3 ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.get_integral_evaluator3)
  MPQC::IntegralEvaluator3 eval = MPQC::IntegralEvaluator3::_create();
  string package = std::string( package_param_->getValueString() );
  eval.set_integral_package( package );
  //eval.initialize_by_name( molecule_, basis_name_, label, max_deriv );
  eval.initialize( bs1, bs2, bs3, label, max_deriv );
  return eval;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.get_integral_evaluator3)
}

/**
 * Method:  get_integral_evaluator4[]
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator4
MPQC::IntegralEvaluatorFactory_impl::get_integral_evaluator4 (
  /*in*/ const ::std::string& label,
  /*in*/ int64_t max_deriv,
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs3,
  /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs4 ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.get_integral_evaluator4)
  MPQC::IntegralEvaluator4 eval = MPQC::IntegralEvaluator4::_create();
  string package = std::string( package_param_->getValueString() );
  eval.set_integral_package( package );
  //eval.initialize_by_name( molecule_, basis_name_, label, max_deriv );
  eval.initialize( bs1, bs2, bs3, bs4, label, max_deriv );
  return eval;
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.get_integral_evaluator4)
}

/**
 * Method:  get_contraction_transform[]
 */
::Chemistry::QC::GaussianBasis::ContractionTransform
MPQC::IntegralEvaluatorFactory_impl::get_contraction_transform () 
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluatorFactory.get_contraction_transform)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory.get_contraction_transform)
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
                                      "package", "intv3");
  pp->addRequest(package_param_);

  return pp;
}

// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluatorFactory._misc)

