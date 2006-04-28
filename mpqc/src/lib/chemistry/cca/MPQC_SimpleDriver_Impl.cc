// 
// File:          MPQC_SimpleDriver_Impl.cc
// Symbol:        MPQC.SimpleDriver-v0.2
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for MPQC.SimpleDriver
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// xml-url       = /home/jpkenny/src/mpqc-libint2.build-shared/src/lib/chemistry/cca/../../../../lib/cca/repo/MPQC.SimpleDriver-v0.2.xml
// 
#include "MPQC_SimpleDriver_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver._includes)
#include <iostream>
#include <iomanip>
#include <fstream>
#include <unistd.h>
#include <cstdio>

using namespace std;
// DO-NOT-DELETE splicer.end(MPQC.SimpleDriver._includes)

// user-defined constructor.
void MPQC::SimpleDriver_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver._ctor)
  // Insert-Code-Here {MPQC.SimpleDriver._ctor} (constructor)
  // DO-NOT-DELETE splicer.end(MPQC.SimpleDriver._ctor)
}

// user-defined destructor.
void MPQC::SimpleDriver_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver._dtor)
  // Insert-Code-Here {MPQC.SimpleDriver._dtor} (destructor)
  // DO-NOT-DELETE splicer.end(MPQC.SimpleDriver._dtor)
}

// static class initializer.
void MPQC::SimpleDriver_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver._load)
  // Insert-Code-Here {MPQC.SimpleDriver._load} (class initialization)
  // DO-NOT-DELETE splicer.end(MPQC.SimpleDriver._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
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
MPQC::SimpleDriver_impl::setServices (
  /* in */ ::gov::cca::Services services ) 
throw ( 
  ::gov::cca::CCAException
){
  // DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver.setServices)

  services_ = services;
  if (services_._is_nil()) return;

  try {
      services_.addProvidesPort(self, "go","gov.cca.ports.GoPort",
                                0);
      services_.registerUsesPort("ModelFactory", "Chemistry.QC.ModelFactory",
                                 0);
      }
  catch (gov::cca::CCAException e) {
      std::cout << "Error using services: "
                << e.getNote() << std::endl;
  }

  // setup parameters
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

  return;

  // DO-NOT-DELETE splicer.end(MPQC.SimpleDriver.setServices)
}

/**
 * Execute some encapsulated functionality on the component. 
 * Return 0 if ok, -1 if internal error but component may be 
 * used further, and -2 if error so severe that component cannot
 * be further used safely.
 */
int32_t
MPQC::SimpleDriver_impl::go ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver.go)

  std::cout << "\nSIMPLE CHEMISTRY COMPONENT DRIVER\n";
  std::cout << "----------------------------------\n";

  std::cout << "SIMPLE DRIVER: getting model factory\n";
  factory_ = services_.getPort("ModelFactory");
  if (factory_._is_nil()) {
      std::cout << "didn't get a model factory\n";
      abort();
    }

  std::cout << "SIMPLE DRIVER: getting model\n";
  model_ = factory_.get_model();
  if (model_._is_nil()) {
      std::cout << "error getting model\n";
      abort();
  }

  std::cout << "SIMPLE DRIVER: getting molecule\n";
  molecule_ = model_.get_molecule();
  if (molecule_._is_nil()) {
      std::cout << "error getting molecule\n";
      abort();
  }

  do_grad_ = grad_param_->value;

  std::cout << "SIMPLE DRIVER: Evaluating energy\n";
  int num_coord = molecule_.get_n_atom() * 3;
  sidl::array<double> coor = molecule_.get_coor();
  double energy = model_.get_energy();

  sidl::array<double> cartg;
  if(do_grad_){
    std::cout << "SIMPLE DRIVER: Evaluating gradient\n";
    cartg = model_.get_gradient();
  }

  // log some results, used for validation
  char *init_dir = getenv("CCACHEM_RESULTS_DIR");
  if (init_dir) {
    chdir(init_dir);
    ofstream outfile( "results.txt" );
    outfile << setprecision(9)
            << setiosflags( ios::showpoint | ios::fixed);
    if(outfile) {
      outfile << "FINAL GEOMETRY:\n";
      for(int i=0; i<num_coord; i+=3) {
        outfile << setw(20) << coor.get(i)
                << setw(20) << coor.get(i+1)
                << setw(20) << coor.get(i+2)
                << endl;
      }
      outfile << endl;

      outfile << "FINAL ENERGY:\n";
      outfile << setw(20) << energy << "\n\n";

      if( do_grad_ ) {
        outfile << "FINAL GRADIENT:\n";
        for(int i=0; i<num_coord; i+=3) {
          outfile << setw(20) << cartg.get(i)
                  << setw(20) << cartg.get(i+1)
                  << setw(20) << cartg.get(i+2)
                  << endl;
        }
        outfile << endl;
      }

    }
    else
      cout << "\n\nCouldn't open results.txt for writing\n\n";
    outfile.close();
  }
    else
      cout << "\n\nCCACHEM_RESULTS_DIR environment variable not found,\n"
           << "  don't know where to log results\n\n";

  if(model_._not_nil())
    model_.finalize();
  if(factory_._not_nil())
    factory_.finalize();
  services_.releasePort("ModelFactory");

  return 0;

  // DO-NOT-DELETE splicer.end(MPQC.SimpleDriver.go)
}


// DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver._misc)

ConfigurableParameterPort *
MPQC::SimpleDriver_impl::setup_parameters(ConfigurableParameterFactory *cpf)
{
  ConfigurableParameterPort * pp = cpf->createConfigurableParameterPort();

  pp->setBatchTitle("PortTranslatorStarter Configuration");
  pp->setGroupName("IntegralTest Input");

  grad_param_  = new BoolParameter("do_gradient",
                                   "Perform gradient evaluation",
                                   "do_gradient",0);
  pp->addRequest(grad_param_);

  return pp;
}

// DO-NOT-DELETE splicer.end(MPQC.SimpleDriver._misc)

