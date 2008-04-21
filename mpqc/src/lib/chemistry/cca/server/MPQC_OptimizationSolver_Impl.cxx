// 
// File:          MPQC_OptimizationSolver_Impl.cxx
// Symbol:        MPQC.OptimizationSolver-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.OptimizationSolver
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 
#include "MPQC_OptimizationSolver_Impl.hxx"

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
// DO-NOT-DELETE splicer.begin(MPQC.OptimizationSolver._includes)
#include <iomanip>
#include <iostream>
#include <sstream>
#include <util/misc/ccaenv.h>
using namespace sc;
// DO-NOT-DELETE splicer.end(MPQC.OptimizationSolver._includes)

// speical constructor, used for data wrapping(required).  Do not put code here unless you really know what you're doing!
MPQC::OptimizationSolver_impl::OptimizationSolver_impl() : StubBase(
  reinterpret_cast< void*>(::MPQC::OptimizationSolver::_wrapObj(
  reinterpret_cast< void*>(this))),false) , _wrapped(true){ 
  // DO-NOT-DELETE splicer.begin(MPQC.OptimizationSolver._ctor2)
  // Insert-Code-Here {MPQC.OptimizationSolver._ctor2} (ctor2)
  // DO-NOT-DELETE splicer.end(MPQC.OptimizationSolver._ctor2)
}

// user defined constructor
void MPQC::OptimizationSolver_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.OptimizationSolver._ctor)

  // DO-NOT-DELETE splicer.end(MPQC.OptimizationSolver._ctor)
}

// user defined destructor
void MPQC::OptimizationSolver_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.OptimizationSolver._dtor)
  // Insert-Code-Here {MPQC.OptimizationSolver._dtor} (destructor)
  // DO-NOT-DELETE splicer.end(MPQC.OptimizationSolver._dtor)
}

// static class initializer
void MPQC::OptimizationSolver_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.OptimizationSolver._load)
  // Insert-Code-Here {MPQC.OptimizationSolver._load} (class initialization)
  // DO-NOT-DELETE splicer.end(MPQC.OptimizationSolver._load)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 * Method:  initialize[]
 */
void
MPQC::OptimizationSolver_impl::initialize_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.OptimizationSolver.initialize)
    
  try {
    opt_model_ = sidl::babel_cast<ChemistryOpt::OptimizationModelInterface>(
                  services_.getPort("OptimizationModelInterface") );
  }
  catch(...) {}
  /* if model_ is assigned, then we need to pass to some sort of sc::Function wrapper */

  /* else we do this */
  if( opt_model_._is_nil() ) {
    coor_model_ = sidl::babel_cast<ChemistryOpt::CoordinateModelInterface>(
                    services_.getPort("CoordinateModelInterface") );
    coor_model_.initialize();
    model_ = coor_model_.get_model();
    molecule_ = model_.get_molecule();
  }

  // form keyval for QNewtonOpt
  std::ostringstream input;
  double conv = molecule_.get_units().convert_to("bohr");
  input
    << "  mole<MolecularEnergyCCA>: ( " << std::endl
    << "    molecule<Molecule>: (" << std::endl
    << "      symmetry = auto" << std::endl
    << "      unit = bohr" << std::endl;
  if( molecule_.get_n_pcharge() ) {
    input << "      charge = [ ";
    for(int i=0; i<molecule_.get_n_pcharge(); ++i )
      input << molecule_.get_point_charge(i) << " ";
    input << "]" << std::endl
          << "      include_q = 0" << std::endl
          << "      include_qq = 0" << std::endl;
  }
  input << "    {atoms geometry } = {" << std::endl;
  input << std::setprecision(16);
  for(int i=0;i<molecule_.get_n_pcharge();++i) {
    input << "\t" << "Q"
      << "\t[  " << molecule_.get_pcharge_cart_coor(i,0)*conv
      << "  " << molecule_.get_pcharge_cart_coor(i,1)*conv
      << "  " << molecule_.get_pcharge_cart_coor(i,2)*conv << "  ]\n";
  }
  for(int i=0;i<molecule_.get_n_atom();++i) {
    input << "\t" << molecule_.get_atomic_number(i)
      << "\t[  " << molecule_.get_cart_coor(i,0)*conv
      << "  " << molecule_.get_cart_coor(i,1)*conv
      << "  " << molecule_.get_cart_coor(i,2)*conv << "  ]\n";
  }
  input << "    }\n";
  input << "    )" << std::endl;
  input << "  )" << std::endl;
  input << "  opt<QNewtonOpt>: (" << std::endl;
  input << "    function = $..:mole" << std::endl;
  input << "    update<BFGSUpdate>: ()" << std::endl;
  input << "    convergence<MolEnergyConvergence>: (" << std::endl;
  input << "      cartesian = yes" << std::endl;
  input << "      energy = $..:..:mole" << std::endl;
  input << "    )" << std::endl;
  input << "  )" << std::endl;

  std::cout << "  opt input:" << std::endl << input.str() << std::endl;

  sc::Ref<sc::ParsedKeyVal> kv = new sc::ParsedKeyVal();
  kv->parse_string(input.str().c_str());
  sc::Ref<sc::DescribedClass> dc;
  try {
    dc = kv->describedclassvalue("opt");
  }
  catch (std::exception &e) {
    e.what();
  }
  opt_ = dynamic_cast<sc::QNewtonOpt*>(dc.pointer());
  
  // DO-NOT-DELETE splicer.end(MPQC.OptimizationSolver.initialize)
}

/**
 * Method:  set_tolerances[]
 */
void
MPQC::OptimizationSolver_impl::set_tolerances_impl (
  /* in */double fatol,
  /* in */double frtol,
  /* in */double catol,
  /* in */double crtol ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.OptimizationSolver.set_tolerances)
  // Insert-Code-Here {MPQC.OptimizationSolver.set_tolerances} (set_tolerances method)
    
    // DO-DELETE-WHEN-IMPLEMENTING exception.begin()
    /*
     * This method has not been implemented
     */
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "set_tolerances");
    throw ex;
    // DO-DELETE-WHEN-IMPLEMENTING exception.end()
    
  // DO-NOT-DELETE splicer.end(MPQC.OptimizationSolver.set_tolerances)
}

/**
 * Method:  set_max_iterations[]
 */
void
MPQC::OptimizationSolver_impl::set_max_iterations_impl (
  /* in */int32_t maxits ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.OptimizationSolver.set_max_iterations)
  // Insert-Code-Here {MPQC.OptimizationSolver.set_max_iterations} (set_max_iterations method)
    
    // DO-DELETE-WHEN-IMPLEMENTING exception.begin()
    /*
     * This method has not been implemented
     */
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "set_max_iterations");
    throw ex;
    // DO-DELETE-WHEN-IMPLEMENTING exception.end()
    
  // DO-NOT-DELETE splicer.end(MPQC.OptimizationSolver.set_max_iterations)
}

/**
 * Method:  solve[]
 */
int32_t
MPQC::OptimizationSolver_impl::solve_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.OptimizationSolver.solve)

  opt_->optimize();
    
  // DO-NOT-DELETE splicer.end(MPQC.OptimizationSolver.solve)
}

/**
 *  Starts up a component presence in the calling framework.
 * @param services the component instance's handle on the framework world.
 * Contracts concerning services and setServices:
 * 
 * The component interaction with the CCA framework
 * and Ports begins on the call to setServices by the framework.
 * 
 * This function is called exactly once for each instance created
 * by the framework.
 * 
 * The argument services will never be nil/null.
 * 
 * Those uses ports which are automatically connected by the framework
 * (so-called service-ports) may be obtained via getPort during
 * setServices.
 */
void
MPQC::OptimizationSolver_impl::setServices_impl (
  /* in */::gov::cca::Services services ) 
// throws:
//     ::gov::cca::CCAException
//     ::sidl::RuntimeException
{
  // DO-NOT-DELETE splicer.begin(MPQC.OptimizationSolver.setServices)
    
  std::cerr << "in solver setServices\n";
  services_ = services;
  if (services_._is_nil()) return;

  try {
      services_.addProvidesPort( *this, "go","gov.cca.ports.GoPort",
                                 0);
      services_.addProvidesPort( *this, "SolverInterface",
        "ChemistryOpt.SolverInterface", 0);
      services_.registerUsesPort("OptimizationModelInterface",
        "ChemistryOpt.OptimizationModel", 0);
      services_.registerUsesPort("CoordinateModelInterface",
        "ChemistryOpt.CoordinateModelInterface", 0);
  }
  catch (gov::cca::CCAException e) {
      std::cout << "Error using services: "
                << e.getNote() << std::endl;
  }

  if( !CCAEnv::initialized() ) {
    Ref<CCAFramework> fw = new ExternalCCAFramework(services_);
    CCAEnv::init( fw );
  }

  std::cerr << "solver setServices complete\n";
    
  // DO-NOT-DELETE splicer.end(MPQC.OptimizationSolver.setServices)
}

/**
 *  
 * Execute some encapsulated functionality on the component. 
 * Return 0 if ok, -1 if internal error but component may be 
 * used further, and -2 if error so severe that component cannot
 * be further used safely.
 */
int32_t
MPQC::OptimizationSolver_impl::go_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.OptimizationSolver.go)

  initialize();
  solve();
  coor_model_.finalize();
    
  // DO-NOT-DELETE splicer.end(MPQC.OptimizationSolver.go)
}


// DO-NOT-DELETE splicer.begin(MPQC.OptimizationSolver._misc)
// Insert-Code-Here {MPQC.OptimizationSolver._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(MPQC.OptimizationSolver._misc)

