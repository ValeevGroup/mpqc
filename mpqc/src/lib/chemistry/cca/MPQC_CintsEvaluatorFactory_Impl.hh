// 
// File:          MPQC_CintsEvaluatorFactory_Impl.hh
// Symbol:        MPQC.CintsEvaluatorFactory-v0.2
// Symbol Type:   class
// Babel Version: 0.10.2
// Description:   Server-side implementation for MPQC.CintsEvaluatorFactory
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.2
// 

#ifndef included_MPQC_CintsEvaluatorFactory_Impl_hh
#define included_MPQC_CintsEvaluatorFactory_Impl_hh

#ifndef included_sidl_cxx_hh
#include "sidl_cxx.hh"
#endif
#ifndef included_MPQC_CintsEvaluatorFactory_IOR_h
#include "MPQC_CintsEvaluatorFactory_IOR.h"
#endif
// 
// Includes for all method dependencies.
// 
#ifndef included_Chemistry_QC_GaussianBasis_CompositeIntegralDescr_hh
#include "Chemistry_QC_GaussianBasis_CompositeIntegralDescr.hh"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralEvaluator1_hh
#include "Chemistry_QC_GaussianBasis_IntegralEvaluator1.hh"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralEvaluator2_hh
#include "Chemistry_QC_GaussianBasis_IntegralEvaluator2.hh"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralEvaluator3_hh
#include "Chemistry_QC_GaussianBasis_IntegralEvaluator3.hh"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralEvaluator4_hh
#include "Chemistry_QC_GaussianBasis_IntegralEvaluator4.hh"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_Molecular_hh
#include "Chemistry_QC_GaussianBasis_Molecular.hh"
#endif
#ifndef included_MPQC_CintsEvaluatorFactory_hh
#include "MPQC_CintsEvaluatorFactory.hh"
#endif
#ifndef included_gov_cca_CCAException_hh
#include "gov_cca_CCAException.hh"
#endif
#ifndef included_gov_cca_Services_hh
#include "gov_cca_Services.hh"
#endif
#ifndef included_sidl_BaseException_hh
#include "sidl_BaseException.hh"
#endif
#ifndef included_sidl_BaseInterface_hh
#include "sidl_BaseInterface.hh"
#endif
#ifndef included_sidl_ClassInfo_hh
#include "sidl_ClassInfo.hh"
#endif


// DO-NOT-DELETE splicer.begin(MPQC.CintsEvaluatorFactory._includes)
// Insert-Code-Here {MPQC.CintsEvaluatorFactory._includes} (includes or arbitrary code)
// DO-NOT-DELETE splicer.end(MPQC.CintsEvaluatorFactory._includes)

namespace MPQC { 

  /**
   * Symbol "MPQC.CintsEvaluatorFactory" (version 0.2)
   */
  class CintsEvaluatorFactory_impl
  // DO-NOT-DELETE splicer.begin(MPQC.CintsEvaluatorFactory._inherits)
  // Insert-Code-Here {MPQC.CintsEvaluatorFactory._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(MPQC.CintsEvaluatorFactory._inherits)
  {

  private:
    // Pointer back to IOR.
    // Use this to dispatch back through IOR vtable.
    CintsEvaluatorFactory self;

    // DO-NOT-DELETE splicer.begin(MPQC.CintsEvaluatorFactory._implementation)
    // Insert-Code-Here {MPQC.CintsEvaluatorFactory._implementation} (additional details)
    // DO-NOT-DELETE splicer.end(MPQC.CintsEvaluatorFactory._implementation)

  private:
    // private default constructor (required)
    CintsEvaluatorFactory_impl() 
    {} 

  public:
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    CintsEvaluatorFactory_impl( struct MPQC_CintsEvaluatorFactory__object * s ) 
      : self(s,true) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~CintsEvaluatorFactory_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // static class initializer
    static void _load();

  public:


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
    setServices (
      /* in */ ::gov::cca::Services services
    )
    throw ( 
      ::gov::cca::CCAException
    );

    /**
     * user defined non-static method.
     */
    ::std::string
    get_name() throw ( 
      ::sidl::BaseException
    );
    /**
     * user defined non-static method.
     */
    ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr
    get_descriptor() throw ( 
      ::sidl::BaseException
    );
    /**
     * user defined non-static method.
     */
    ::sidl::array<int32_t>
    get_max_deriv_lvl (
      /* in */ ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr desc
    )
    throw ( 
      ::sidl::BaseException
    );


    /**
     * Set available storage
     * @param storage Available storage in bytes 
     */
    void
    set_storage (
      /* in */ int64_t storage
    )
    throw ( 
      ::sidl::BaseException
    );


    /**
     * Get a 1-center integral evaluator
     * @param desc Integral set descriptor
     * @return 1-center integral evaluator 
     */
    ::Chemistry::QC::GaussianBasis::IntegralEvaluator1
    get_evaluator1 (
      /* in */ ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr desc,
      /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1
    )
    throw ( 
      ::sidl::BaseException
    );


    /**
     * Get a 2-center integral evaluator
     * @param desc Integral set descriptor
     * @return 2-center integral evaluator 
     */
    ::Chemistry::QC::GaussianBasis::IntegralEvaluator2
    get_evaluator2 (
      /* in */ ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr desc,
      /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
      /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2
    )
    throw ( 
      ::sidl::BaseException
    );


    /**
     * Get a 3-center integral evaluator
     * @param desc Integral set descriptor
     * @return 3-center integral evaluator 
     */
    ::Chemistry::QC::GaussianBasis::IntegralEvaluator3
    get_evaluator3 (
      /* in */ ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr desc,
      /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
      /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2,
      /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs3
    )
    throw ( 
      ::sidl::BaseException
    );


    /**
     * Get a 4-center integral evaluator
     * @param desc Integral set descriptor
     * @return 4-center integral evaluator 
     */
    ::Chemistry::QC::GaussianBasis::IntegralEvaluator4
    get_evaluator4 (
      /* in */ ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr desc,
      /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
      /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2,
      /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs3,
      /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs4
    )
    throw ( 
      ::sidl::BaseException
    );

  };  // end class CintsEvaluatorFactory_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.CintsEvaluatorFactory._misc)
// Insert-Code-Here {MPQC.CintsEvaluatorFactory._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(MPQC.CintsEvaluatorFactory._misc)

#endif
