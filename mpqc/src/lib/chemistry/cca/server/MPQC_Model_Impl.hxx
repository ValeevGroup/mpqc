// 
// File:          MPQC_Model_Impl.hxx
// Symbol:        MPQC.Model-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.Model
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 

#ifndef included_MPQC_Model_Impl_hxx
#define included_MPQC_Model_Impl_hxx

#ifndef included_sidl_cxx_hxx
#include "sidl_cxx.hxx"
#endif
#ifndef included_MPQC_Model_IOR_h
#include "MPQC_Model_IOR.h"
#endif
#ifndef included_Chemistry_MoleculeInterface_hxx
#include "Chemistry_MoleculeInterface.hxx"
#endif
#ifndef included_Chemistry_QC_ModelInterface_hxx
#include "Chemistry_QC_ModelInterface.hxx"
#endif
#ifndef included_MPQC_Model_hxx
#include "MPQC_Model.hxx"
#endif
#ifndef included_gov_cca_TypeMap_hxx
#include "gov_cca_TypeMap.hxx"
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


// DO-NOT-DELETE splicer.begin(MPQC.Model._includes)

#include <chemistry/qc/wfn/wfn.h>
#include <ChemistryCXX_Molecule.hxx>

// DO-NOT-DELETE splicer.end(MPQC.Model._includes)

namespace MPQC { 

  /**
   * Symbol "MPQC.Model" (version 0.2)
   */
  class Model_impl : public virtual ::MPQC::Model 
  // DO-NOT-DELETE splicer.begin(MPQC.Model._inherits)
  // Insert-Code-Here {MPQC.Model._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(MPQC.Model._inherits)
  {

  // All data marked protected will be accessable by 
  // descendant Impl classes
  protected:

    bool _wrapped;

    // DO-NOT-DELETE splicer.begin(MPQC.Model._implementation)

    sc::Ref<sc::Wavefunction> wfn_;
    Chemistry::MoleculeInterface molecule_;
    gov::cca::TypeMap cqos_tm_;
    double time_get_energy_;

    // DO-NOT-DELETE splicer.end(MPQC.Model._implementation)

  public:
    // default constructor, used for data wrapping(required)
    Model_impl();
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    Model_impl( struct MPQC_Model__object * s ) : StubBase(s,true), _wrapped(
      false) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~Model_impl() { _dtor(); }

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
    initialize_parsedkeyval_impl (
      /* in */const ::std::string& keyword,
      /* in */const ::std::string& input
    )
    ;

    /**
     * user defined non-static method.
     */
    void
    initialize_parsedkeyval_file_impl (
      /* in */const ::std::string& keyword,
      /* in */const ::std::string& filename
    )
    ;

    /**
     * user defined non-static method.
     */
    void
    initialize_aggregatekeyval_impl (
      /* in */const ::std::string& keyword,
      /* in */const ::std::string& input,
      /* in */void* describedclass
    )
    ;

    /**
     * user defined non-static method.
     */
    void
    initialize_pointer_impl (
      /* in */void* ptr
    )
    ;

    /**
     * user defined non-static method.
     */
    double
    get_energy_impl() ;

    /**
     *  Sets the accuracy for subsequent energy calculations.
     * @param acc The new accuracy. 
     */
    void
    set_energy_accuracy_impl (
      /* in */double acc
    )
    ;


    /**
     *  Returns the accuracy to which the energy is already computed.
     * The result is undefined if the energy has not already
     * been computed.
     * @return The energy accuracy. 
     */
    double
    get_energy_accuracy_impl() ;

    /**
     *  This allows a programmer to request that if any result
     * is computed,
     * then the energy is computed too.  This allows, say, for a request
     * for a gradient to cause the energy to be computed.  This computed
     * energy is cached and returned when the get_energy() member
     * is called.
     * @param doit Whether or not to compute the energy.
     */
    void
    set_do_energy_impl (
      /* in */bool doit
    )
    ;


    /**
     *  Returns the Cartesian gradient.  
     */
    ::sidl::array<double>
    get_gradient_impl() ;

    /**
     *  Sets the accuracy for subsequent gradient calculations
     * @param acc The new accuracy for gradients. 
     */
    void
    set_gradient_accuracy_impl (
      /* in */double acc
    )
    ;


    /**
     *  Returns the accuracy to which the gradient is already computed.
     * The result is undefined if the gradient has not already
     * been computed.
     * @return The current gradient accuracy. 
     */
    double
    get_gradient_accuracy_impl() ;

    /**
     *  Returns the Cartesian Hessian. @return The Hessian. 
     */
    ::sidl::array<double>
    get_hessian_impl() ;

    /**
     *  Sets the accuracy for subsequent Hessian calculations.
     * @param acc The new accuracy for Hessians. 
     */
    void
    set_hessian_accuracy_impl (
      /* in */double acc
    )
    ;


    /**
     *  Returns the accuracy to which the Hessian is already computed.
     * The result is undefined if the Hessian has not already
     * been computed. 
     */
    double
    get_hessian_accuracy_impl() ;

    /**
     *  Returns a Cartesian guess Hessian. 
     */
    ::sidl::array<double>
    get_guess_hessian_impl() ;

    /**
     *  Sets the accuracy for subsequent guess Hessian calculations.
     * @param acc The new accuracy for guess Hessians. 
     */
    void
    set_guess_hessian_accuracy_impl (
      /* in */double acc
    )
    ;


    /**
     *  Returns the accuracy to which the guess Hessian is
     * already computed.  The result is undefined if the guess Hessian
     * has not already been computed.
     * @return The guess hessian accuracy.  
     */
    double
    get_guess_hessian_accuracy_impl() ;

    /**
     *  This should be called when the object is no longer needed.
     * No other members may be called after finalize. 
     */
    int32_t
    finalize_impl() ;

    /**
     *  Set the molecule. @param molecule The new molecule. 
     */
    void
    set_molecule_impl (
      /* in */::Chemistry::MoleculeInterface molecule
    )
    ;


    /**
     *  Returns the molecule.  @return The Molecule object. 
     */
    ::Chemistry::MoleculeInterface
    get_molecule_impl() ;

    /**
     *  Sets the initial CQoS metadata typemap.
     * The model may augment this typemap.
     * @param typemap The initial typemap. 
     */
    void
    set_metadata_impl (
      /* in */::gov::cca::TypeMap typemap
    )
    ;


    /**
     *  Returns CQoS metadata typemap.
     * @return Metadata typemap. 
     */
    ::gov::cca::TypeMap
    get_metadata_impl() ;
  };  // end class Model_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.Model._misc)
// Insert-Code-Here {MPQC.Model._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(MPQC.Model._misc)

#endif
