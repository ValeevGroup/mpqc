// 
// File:          MPQC_ModelFactory_Impl.hxx
// Symbol:        MPQC.ModelFactory-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.ModelFactory
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 

#ifndef included_MPQC_ModelFactory_Impl_hxx
#define included_MPQC_ModelFactory_Impl_hxx

#ifndef included_sidl_cxx_hxx
#include "sidl_cxx.hxx"
#endif
#ifndef included_MPQC_ModelFactory_IOR_h
#include "MPQC_ModelFactory_IOR.h"
#endif
#ifndef included_Chemistry_MoleculeInterface_hxx
#include "Chemistry_MoleculeInterface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralEvaluatorFactoryInterface_hxx
#include "Chemistry_QC_GaussianBasis_IntegralEvaluatorFactoryInterface.hxx"
#endif
#ifndef included_Chemistry_QC_ModelFactoryInterface_hxx
#include "Chemistry_QC_ModelFactoryInterface.hxx"
#endif
#ifndef included_Chemistry_QC_ModelInterface_hxx
#include "Chemistry_QC_ModelInterface.hxx"
#endif
#ifndef included_MPQC_ModelFactory_hxx
#include "MPQC_ModelFactory.hxx"
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


// DO-NOT-DELETE splicer.begin(MPQC.ModelFactory._includes)

#include <string>
#include "ChemistryCXX_MoleculeFactory.hxx"
#include <util/group/message.h>
#include <util/group/memory.h>
#include <util/group/mstate.h>
#include <util/group/thread.h>
#include <util/group/pregtime.h>
#include <chemistry/cca/int/intcca.h>
#include <chemistry/qc/basis/integral.h>
#include <gov_cca_ports_ParameterPortFactory.hxx>
#include <gov_cca_ports_ParameterPort.hxx>

// DO-NOT-DELETE splicer.end(MPQC.ModelFactory._includes)

namespace MPQC { 

  /**
   * Symbol "MPQC.ModelFactory" (version 0.2)
   */
  class ModelFactory_impl : public virtual ::MPQC::ModelFactory 
  // DO-NOT-DELETE splicer.begin(MPQC.ModelFactory._inherits)
  // Insert-Code-Here {MPQC.ModelFactory._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(MPQC.ModelFactory._inherits)
  {

  // All data marked protected will be accessable by 
  // descendant Impl classes
  protected:

    bool _wrapped;

    // DO-NOT-DELETE splicer.begin(MPQC.ModelFactory._implementation)

      std::string theory_;
      std::string basis_;
      std::string molecule_filename_;
      bool use_cca_integrals_;

      gov::cca::Services services_;
      gov::cca::TypeMap tm_;
      gov::cca::ports::ParameterPortFactory ppf_;
      gov::cca::ports::ParameterPort pp_;

      Chemistry::MoleculeFactoryInterface molecule_factory_;
      Chemistry::MoleculeInterface molecule_;
      Chemistry::QC::GaussianBasis::IntegralEvaluatorFactoryInterface eval_factory_;

      sc::Ref<sc::MessageGrp> grp_;
      sc::Ref<sc::ThreadGrp> thread_;
      sc::Ref<sc::MemoryGrp> memory_;

      sc::Ref<sc::IntegralCCA> intcca_;

    // DO-NOT-DELETE splicer.end(MPQC.ModelFactory._implementation)

  public:
    // default constructor, used for data wrapping(required)
    ModelFactory_impl();
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    ModelFactory_impl( struct MPQC_ModelFactory__object * s ) : StubBase(s,
      true), _wrapped(false) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~ModelFactory_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // true if this object was created by a user newing the impl
    inline bool _isWrapped() {return _wrapped;}

    // static class initializer
    static void _load();

  public:


    /**
     *  Set the theory name for Model's created with get_model.
     * @param theory A string giving the name of the theory, 
     * for example, B3LYP.
     */
    void
    set_theory_impl (
      /* in */const ::std::string& theory
    )
    ;


    /**
     *  Set the basis set name for Model's created with get_model.
     * @param basis The basis set name to use, for example, aug-cc-pVDZ.
     */
    void
    set_basis_impl (
      /* in */const ::std::string& basis
    )
    ;


    /**
     *  Set the Molecule to use for Model's created with get_model.
     * @param molecule An object of type Molecule.
     */
    void
    set_molecule_impl (
      /* in */::Chemistry::MoleculeInterface molecule
    )
    ;


    /**
     *  Set the object to use to compute integrals for Model's 
     * created with get_model.
     * @param intfact An object of type 
     * GaussianBasis.IntegralEvaluatorFactory.
     */
    void
    set_integral_factory_impl (
      /* in */::Chemistry::QC::GaussianBasis::IntegralEvaluatorFactoryInterface 
        intfact
    )
    ;


    /**
     *  Returns a newly created Model.  Before get_model can be called, 
     * set_theory, set_basis, and set_molecule must be called.
     * @return The new Model instance.
     */
    ::Chemistry::QC::ModelInterface
    get_model_impl() ;

    /**
     *  This should be called when the object is no longer needed. 
     * No other members may be called after finalize. 
     */
    int32_t
    finalize_impl() ;

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

  };  // end class ModelFactory_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.ModelFactory._misc)
// Insert-Code-Here {MPQC.ModelFactory._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(MPQC.ModelFactory._misc)

#endif
