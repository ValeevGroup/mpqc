// 
// File:          MPQC_Libint2EvaluatorFactory_Impl.hxx
// Symbol:        MPQC.Libint2EvaluatorFactory-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.Libint2EvaluatorFactory
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 

#ifndef included_MPQC_Libint2EvaluatorFactory_Impl_hxx
#define included_MPQC_Libint2EvaluatorFactory_Impl_hxx

#ifndef included_sidl_cxx_hxx
#include "sidl_cxx.hxx"
#endif
#ifndef included_MPQC_Libint2EvaluatorFactory_IOR_h
#include "MPQC_Libint2EvaluatorFactory_IOR.h"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_CompositeIntegralDescrInterface_hxx
#include "Chemistry_QC_GaussianBasis_CompositeIntegralDescrInterface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralDescrInterface_hxx
#include "Chemistry_QC_GaussianBasis_IntegralDescrInterface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralEvaluator1Interface_hxx
#include "Chemistry_QC_GaussianBasis_IntegralEvaluator1Interface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralEvaluator2Interface_hxx
#include "Chemistry_QC_GaussianBasis_IntegralEvaluator2Interface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralEvaluator3Interface_hxx
#include "Chemistry_QC_GaussianBasis_IntegralEvaluator3Interface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralEvaluator4Interface_hxx
#include "Chemistry_QC_GaussianBasis_IntegralEvaluator4Interface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralEvaluatorFactoryInterface_hxx
#include "Chemistry_QC_GaussianBasis_IntegralEvaluatorFactoryInterface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_MolecularInterface_hxx
#include "Chemistry_QC_GaussianBasis_MolecularInterface.hxx"
#endif
#ifndef included_MPQC_Libint2EvaluatorFactory_hxx
#include "MPQC_Libint2EvaluatorFactory.hxx"
#endif
#ifndef included_gov_cca_CCAException_hxx
#include "gov_cca_CCAException.hxx"
#endif
#ifndef included_gov_cca_Component_hxx
#include "gov_cca_Component.hxx"
#endif
#ifndef included_gov_cca_Services_hxx
#include "gov_cca_Services.hxx"
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


// DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory._includes)
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/libint2/libint2.h>
// DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory._includes)

namespace MPQC { 

  /**
   * Symbol "MPQC.Libint2EvaluatorFactory" (version 0.2)
   */
  class Libint2EvaluatorFactory_impl : public virtual 
    ::MPQC::Libint2EvaluatorFactory 
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory._inherits)
  // Insert-Code-Here {MPQC.Libint2EvaluatorFactory._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory._inherits)
  {

  // All data marked protected will be accessable by 
  // descendant Impl classes
  protected:

    bool _wrapped;

    // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory._implementation)

    Chemistry::QC::GaussianBasis::MolecularInterface 
      basis1_, basis2_, basis3_, basis4_;
    
    /// registry of known integral types
    Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface 
      cdesc_;
    /// registry of known integral types for which derivatives cannot be computed
    Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface 
      cdesc_no_deriv_;
    
    std::vector< Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface > 
      comps_;
    std::vector< sc::Ref<sc::TwoBodyInt> > grts_;
    gov::cca::Services services_;
    int storage_;
    sc::Ref<sc::IntegralLibint2> integral_;
    sc::Ref<sc::DipoleData> dipole_data_;
    sc::Ref<sc::DipoleData> quad_data_;

    std::vector< std::vector< sc::Ref<sc::OneBodyOneCenterInt> > > 
      obocint_vec_;
    std::vector< std::vector< sc::Ref<sc::OneBodyInt> > > 
      obint_vec_;
    std::vector< std::vector< sc::Ref<sc::OneBodyDerivInt> > > 
      obderivint_vec_;
    std::vector< std::vector< sc::Ref<sc::TwoBodyTwoCenterInt> > > 
      tbtcint_vec_;
    std::vector< std::vector< sc::Ref<sc::TwoBodyThreeCenterInt> > > 
      tb3cint_vec_;
    std::vector< std::vector< sc::Ref<sc::TwoBodyInt> > > 
      tbint_vec_;
    std::vector< std::vector< sc::Ref<sc::TwoBodyDerivInt> > > 
      tbderivint_vec_;
    
    // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory._implementation)

  public:
    // default constructor, used for data wrapping(required)
    Libint2EvaluatorFactory_impl();
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    Libint2EvaluatorFactory_impl( struct MPQC_Libint2EvaluatorFactory__object * 
      s ) : StubBase(s,true), _wrapped(false) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~Libint2EvaluatorFactory_impl() { _dtor(); }

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
    ::std::string
    get_name_impl() ;
    /**
     * user defined non-static method.
     */
    ::Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface
    get_descriptor_impl() ;
    /**
     * user defined non-static method.
     */
    bool
    is_supported_impl (
      /* in */::Chemistry::QC::GaussianBasis::IntegralDescrInterface desc
    )
    ;


    /**
     *  Set available storage
     * @param storage Available storage in bytes 
     */
    void
    set_storage_impl (
      /* in */int64_t storage
    )
    ;


    /**
     *  Get a 1-center integral evaluator
     * @param desc Integral set descriptor
     * @return 1-center integral evaluator 
     */
    ::Chemistry::QC::GaussianBasis::IntegralEvaluator1Interface
    get_evaluator1_impl (
      /* in */::Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface 
        desc,
      /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs1
    )
    ;


    /**
     *  Get a 2-center integral evaluator
     * @param desc Integral set descriptor
     * @return 2-center integral evaluator 
     */
    ::Chemistry::QC::GaussianBasis::IntegralEvaluator2Interface
    get_evaluator2_impl (
      /* in */::Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface 
        desc,
      /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs1,
      /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs2
    )
    ;


    /**
     *  Get a 3-center integral evaluator
     * @param desc Integral set descriptor
     * @return 3-center integral evaluator 
     */
    ::Chemistry::QC::GaussianBasis::IntegralEvaluator3Interface
    get_evaluator3_impl (
      /* in */::Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface 
        desc,
      /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs1,
      /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs2,
      /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs3
    )
    ;


    /**
     *  Get a 4-center integral evaluator
     * @param desc Integral set descriptor
     * @return 4-center integral evaluator 
     */
    ::Chemistry::QC::GaussianBasis::IntegralEvaluator4Interface
    get_evaluator4_impl (
      /* in */::Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface 
        desc,
      /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs1,
      /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs2,
      /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs3,
      /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs4
    )
    ;


    /**
     *  This should be called when the object is no longer needed.
     * No other members may be called after finalize. 
     */
    int32_t
    finalize_impl() ;

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
    setServices_impl (
      /* in */::gov::cca::Services services
    )
    // throws:
    //     ::gov::cca::CCAException
    //     ::sidl::RuntimeException
    ;

  };  // end class Libint2EvaluatorFactory_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory._misc)
// Insert-Code-Here {MPQC.Libint2EvaluatorFactory._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory._misc)

#endif
