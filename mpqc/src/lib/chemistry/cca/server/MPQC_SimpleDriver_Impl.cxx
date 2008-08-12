// 
// File:          MPQC_SimpleDriver_Impl.cxx
// Symbol:        MPQC.SimpleDriver-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.SimpleDriver
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 
#include "MPQC_SimpleDriver_Impl.hxx"

// 
// Includes for all method dependencies.
// 
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
// DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver._includes)
#include <iostream>
#include <iomanip>
#include <fstream>
#include <unistd.h>
#include <cstdio>

#include "MPQC_CCAException.hxx"

using namespace std;
// DO-NOT-DELETE splicer.end(MPQC.SimpleDriver._includes)

// speical constructor, used for data wrapping(required).  Do not put code here unless you really know what you're doing!
MPQC::SimpleDriver_impl::SimpleDriver_impl() : StubBase(reinterpret_cast< 
  void*>(::MPQC::SimpleDriver::_wrapObj(reinterpret_cast< void*>(this))),false) 
  , _wrapped(true){ 
  // DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver._ctor2)
  // Insert-Code-Here {MPQC.SimpleDriver._ctor2} (ctor2)
  // DO-NOT-DELETE splicer.end(MPQC.SimpleDriver._ctor2)
}

// user defined constructor
void MPQC::SimpleDriver_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver._ctor)
  // Insert-Code-Here {MPQC.SimpleDriver._ctor} (constructor)
  // DO-NOT-DELETE splicer.end(MPQC.SimpleDriver._ctor)
}

// user defined destructor
void MPQC::SimpleDriver_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver._dtor)

  services_.unregisterUsesPort("ModelFactoryInterface");
  services_.unregisterUsesPort("ppf");

  services_.removeProvidesPort("go");

  // DO-NOT-DELETE splicer.end(MPQC.SimpleDriver._dtor)
}

// static class initializer
void MPQC::SimpleDriver_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver._load)
  // Insert-Code-Here {MPQC.SimpleDriver._load} (class initialization)
  // DO-NOT-DELETE splicer.end(MPQC.SimpleDriver._load)
}

// user defined static methods: (none)

// user defined non-static methods:
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
MPQC::SimpleDriver_impl::setServices_impl (
  /* in */::gov::cca::Services services ) 
// throws:
//     ::gov::cca::CCAException
//     ::sidl::RuntimeException
{
  // DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver.setServices)

  services_ = services;
  if (services_._is_nil()) return;

  try {
      services_.addProvidesPort( *this, "go","gov.cca.ports.GoPort",
                                0);
      services_.registerUsesPort("ppf",
		                 "gov.cca.ports.ParameterPortFactory",
				 0);
      services_.registerUsesPort("ModelFactoryInterface", 
                                 "Chemistry.QC.ModelFactoryInterface",
                                 0);
      }
  catch (gov::cca::CCAException e) {
      std::cout << "Error using services: "
                << e.getNote() << std::endl;
  }

  // setup parameters
  try {

    if (services_._not_nil()) {
        tm_ = services_.createTypeMap();
        if(tm_._is_nil()) {
            MPQC::CCAException ex = MPQC::CCAException::_create();
            ex.setNote("tm is nil");
            ex.add(__FILE__,__LINE__,"MPQC::SimpleDriver_impl::setServices");
            throw ex;
          }
        ppf_ = sidl::babel_cast<gov::cca::ports::ParameterPortFactory>(
                 services_.getPort("ppf") );
        if (ppf_._is_nil()) {
            MPQC::CCAException ex = MPQC::CCAException::_create();
            ex.setNote("ppf is nil");
            ex.add(__FILE__,__LINE__,"MPQC::SimpleDriver_impl::setServices");
            throw ex;
          }
        ppf_.initParameterData(tm_,"CONFIG");
        ppf_.setBatchTitle(tm_,"MPQC SimpleDriver Options");
        ppf_.setGroupName(tm_,"Job Options");
        ppf_.addRequestBoolean(tm_,"do_gradient",
                               "Perform gradient evaluation?",
                               "Do gradient?",0);
        ppf_.addParameterPort(tm_, services_);
        services_.releasePort("ppf");

        pp_ = sidl::babel_cast<gov::cca::ports::ParameterPort>(
                services_.getPort("CONFIG") );
        if (pp_._is_nil()) {
            std::cerr << "getport failed\n";
            abort();
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
 *  
 * Execute some encapsulated functionality on the component. 
 * Return 0 if ok, -1 if internal error but component may be 
 * used further, and -2 if error so severe that component cannot
 * be further used safely.
 */
int32_t
MPQC::SimpleDriver_impl::go_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver.go)

  gov::cca::TypeMap tm = pp_.readConfigurationMap();
  do_grad_ = tm.getBool("do_gradient","failed do_gradient fetch");

  std::cout << "\nSIMPLE CHEMISTRY COMPONENT DRIVER\n";
  std::cout << "----------------------------------\n";

  std::cout << "SIMPLE DRIVER: getting model factory\n";
  factory_ = sidl::babel_cast<Chemistry::QC::ModelFactoryInterface>(
               services_.getPort("ModelFactoryInterface") );
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

  std::cout << "SIMPLE DRIVER: Evaluating energy\n";
  int num_coord = ( molecule_.get_n_atom() + molecule_.get_n_pcharge() ) * 3;
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
  services_.releasePort("ModelFactoryInterface");
  services_.releasePort("CONFIG");
  ppf_.removeParameterPort(tm_, services_);
  services_.releasePort("ppf");

  return 0;

  // DO-NOT-DELETE splicer.end(MPQC.SimpleDriver.go)
}


// DO-NOT-DELETE splicer.begin(MPQC.SimpleDriver._misc)

// DO-NOT-DELETE splicer.end(MPQC.SimpleDriver._misc)

