// 
// File:          MPQC_IntegralEvaluator1_Impl.hh
// Symbol:        MPQC.IntegralEvaluator1-v0.2
// Symbol Type:   class
// Babel Version: 0.10.2
// Description:   Server-side implementation for MPQC.IntegralEvaluator1
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.2
// 

#ifndef included_MPQC_IntegralEvaluator1_Impl_hh
#define included_MPQC_IntegralEvaluator1_Impl_hh

#ifndef included_sidl_cxx_hh
#include "sidl_cxx.hh"
#endif
#ifndef included_MPQC_IntegralEvaluator1_IOR_h
#include "MPQC_IntegralEvaluator1_IOR.h"
#endif
// 
// Includes for all method dependencies.
// 
#ifndef included_Chemistry_QC_GaussianBasis_DerivCenters_hh
#include "Chemistry_QC_GaussianBasis_DerivCenters.hh"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_Molecular_hh
#include "Chemistry_QC_GaussianBasis_Molecular.hh"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_ObIntEvalType_hh
#include "Chemistry_QC_GaussianBasis_ObIntEvalType.hh"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_Package_hh
#include "Chemistry_QC_GaussianBasis_Package.hh"
#endif
#ifndef included_MPQC_IntegralEvaluator1_hh
#include "MPQC_IntegralEvaluator1.hh"
#endif
#ifndef included_sidl_BaseInterface_hh
#include "sidl_BaseInterface.hh"
#endif
#ifndef included_sidl_ClassInfo_hh
#include "sidl_ClassInfo.hh"
#endif


// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1._includes)
// Insert-Code-Here {MPQC.IntegralEvaluator1._includes} (includes or arbitrary code)
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1._includes)

namespace MPQC { 

  /**
   * Symbol "MPQC.IntegralEvaluator1" (version 0.2)
   */
  class IntegralEvaluator1_impl
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1._inherits)
  // Insert-Code-Here {MPQC.IntegralEvaluator1._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1._inherits)
  {

  private:
    // Pointer back to IOR.
    // Use this to dispatch back through IOR vtable.
    IntegralEvaluator1 self;

    // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1._implementation)
    // Insert-Code-Here {MPQC.IntegralEvaluator1._implementation} (additional details)
    // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1._implementation)

  private:
    // private default constructor (required)
    IntegralEvaluator1_impl() 
    {} 

  public:
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    IntegralEvaluator1_impl( struct MPQC_IntegralEvaluator1__object * s ) : 
      self(s,true) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~IntegralEvaluator1_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // static class initializer
    static void _load();

  public:

    /**
     * user defined non-static method.
     */
    void
    set_integral_package (
      /* in */ ::Chemistry::QC::GaussianBasis::Package type
    )
    throw () 
    ;


    /**
     * Initialize the evaluator.
     * @param bs1 Molecular basis on center 1.
     * @param type ObIntEvalType specifying eval type.
     * @param max_deriv Max derivative to compute.
     * @param storage Available storage in bytes.
     * @param deriv_ctr Derivative center descriptor. 
     */
    void
    obint_initialize (
      /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
      /* in */ ::Chemistry::QC::GaussianBasis::ObIntEvalType type,
      /* in */ int32_t max_deriv,
      /* in */ int64_t storage,
      /* in */ ::Chemistry::QC::GaussianBasis::DerivCenters deriv_ctr
    )
    throw () 
    ;


    /**
     * Get one body int buffer pointer.
     * @return Buffer pointer. 
     */
    void*
    get_obint_buffer() throw () 
    ;

    /**
     * Compute a shell singlet of integrals.  deriv_atom must
     * be used for nuclear derivatives if the operator contains
     * nuclear coordinates, otherwise, set to -1 and use deriv_ctr.
     * @param shellnum1 Gaussian shell number 1.
     * @param deriv_level Derivative level.
     * @param deriv_atom Atom number for derivative
     * (-1 if using DerivCenter). 
     */
    void
    compute (
      /* in */ int64_t shellnum1,
      /* in */ int32_t deriv_level,
      /* in */ int64_t deriv_atom
    )
    throw () 
    ;


    /**
     * Compute a shell singlet of integrals and return as a borrowed
     * sidl array.  deriv_atom must be used for nuclear derivatives if
     * the operator contains nuclear coordinates, otherwise, set to -1
     * and use deriv_ctr.
     * @param shellnum1 Gaussian shell number 1.
     * @param deriv_level Derivative level.
     * @param deriv_atom Atom number for derivative
     * (-1 if using DerivCenter).
     * @return Borrowed sidl array buffer. 
     */
    ::sidl::array<double>
    compute_array (
      /* in */ int64_t shellnum1,
      /* in */ int32_t deriv_level,
      /* in */ int64_t deriv_atom
    )
    throw () 
    ;

  };  // end class IntegralEvaluator1_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator1._misc)
// Insert-Code-Here {MPQC.IntegralEvaluator1._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator1._misc)

#endif
