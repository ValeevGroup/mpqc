// 
// File:          MPQC_OptimizationSolver_Impl.hxx
// Symbol:        MPQC.OptimizationSolver-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.OptimizationSolver
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 

#ifndef included_MPQC_OptimizationSolver_Impl_hxx
#define included_MPQC_OptimizationSolver_Impl_hxx

#ifndef included_sidl_cxx_hxx
#include "sidl_cxx.hxx"
#endif
#ifndef included_MPQC_OptimizationSolver_IOR_h
#include "MPQC_OptimizationSolver_IOR.h"
#endif
#ifndef included_ChemistryOpt_SolverInterface_hxx
#include "ChemistryOpt_SolverInterface.hxx"
#endif
#ifndef included_MPQC_OptimizationSolver_hxx
#include "MPQC_OptimizationSolver.hxx"
#endif
#ifndef included_gov_cca_CCAException_hxx
#include "gov_cca_CCAException.hxx"
#endif
#ifndef included_gov_cca_Component_hxx
#include "gov_cca_Component.hxx"
#endif
#ifndef included_gov_cca_Port_hxx
#include "gov_cca_Port.hxx"
#endif
#ifndef included_gov_cca_Services_hxx
#include "gov_cca_Services.hxx"
#endif
#ifndef included_gov_cca_ports_GoPort_hxx
#include "gov_cca_ports_GoPort.hxx"
#endif
#ifndef included_sidl_BaseClass_hxx
#include "sidl_BaseClass.hxx"
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


// DO-NOT-DELETE splicer.begin(MPQC.OptimizationSolver._includes)
#include <chemistry/cca/molecule/energy.h>
#include <math/optimize/qnewton.h>
#include <Chemistry_MoleculeInterface.hxx>
#include <Chemistry_QC_ModelInterface.hxx>
#include <ChemistryOpt_CoordinateModelInterface.hxx>
#include <ChemistryOpt_OptimizationModelInterface.hxx>
// DO-NOT-DELETE splicer.end(MPQC.OptimizationSolver._includes)

namespace MPQC { 

  /**
   * Symbol "MPQC.OptimizationSolver" (version 0.2)
   */
  class OptimizationSolver_impl : public virtual ::MPQC::OptimizationSolver 
  // DO-NOT-DELETE splicer.begin(MPQC.OptimizationSolver._inherits)
  // Insert-Code-Here {MPQC.OptimizationSolver._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(MPQC.OptimizationSolver._inherits)
  {

  // All data marked protected will be accessable by 
  // descendant Impl classes
  protected:

    bool _wrapped;

    // DO-NOT-DELETE splicer.begin(MPQC.OptimizationSolver._implementation)
    
    gov::cca::Services services_;
    Chemistry::MoleculeInterface molecule_;
    ChemistryOpt::CoordinateModelInterface coor_model_;
    ChemistryOpt::OptimizationModelInterface opt_model_;
    Chemistry::QC::ModelInterface model_;
    sc::Ref<sc::QNewtonOpt> opt_;

    // DO-NOT-DELETE splicer.end(MPQC.OptimizationSolver._implementation)

  public:
    // default constructor, used for data wrapping(required)
    OptimizationSolver_impl();
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    OptimizationSolver_impl( struct MPQC_OptimizationSolver__object * s ) : 
      StubBase(s,true), _wrapped(false) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~OptimizationSolver_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // true if this object was created by a user newing the impl
    inline bool _isWrapped() {return _wrapped;}

    // static class initializer
    static void _load();

  public:

    /**
     * user defined non-static method.
     */
    void
    initialize_impl() ;
    /**
     * user defined non-static method.
     */
    void
    set_tolerances_impl (
      /* in */double fatol,
      /* in */double frtol,
      /* in */double catol,
      /* in */double crtol
    )
    ;

    /**
     * user defined non-static method.
     */
    void
    set_max_iterations_impl (
      /* in */int32_t maxits
    )
    ;

    /**
     * user defined non-static method.
     */
    int32_t
    solve_impl() ;

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
    setServices_impl (
      /* in */::gov::cca::Services services
    )
    // throws:
    //     ::gov::cca::CCAException
    //     ::sidl::RuntimeException
    ;


    /**
     *  
     * Execute some encapsulated functionality on the component. 
     * Return 0 if ok, -1 if internal error but component may be 
     * used further, and -2 if error so severe that component cannot
     * be further used safely.
     */
    int32_t
    go_impl() ;
  };  // end class OptimizationSolver_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.OptimizationSolver._misc)
// Insert-Code-Here {MPQC.OptimizationSolver._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(MPQC.OptimizationSolver._misc)

#endif
