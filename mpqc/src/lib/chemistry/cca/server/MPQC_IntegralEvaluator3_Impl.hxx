// 
// File:          MPQC_IntegralEvaluator3_Impl.hxx
// Symbol:        MPQC.IntegralEvaluator3-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.IntegralEvaluator3
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 

#ifndef included_MPQC_IntegralEvaluator3_Impl_hxx
#define included_MPQC_IntegralEvaluator3_Impl_hxx

#ifndef included_sidl_cxx_hxx
#include "sidl_cxx.hxx"
#endif
#ifndef included_MPQC_IntegralEvaluator3_IOR_h
#include "MPQC_IntegralEvaluator3_IOR.h"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_CompositeIntegralDescrInterface_hxx
#include "Chemistry_QC_GaussianBasis_CompositeIntegralDescrInterface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralDescrInterface_hxx
#include "Chemistry_QC_GaussianBasis_IntegralDescrInterface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralEvaluator3Interface_hxx
#include "Chemistry_QC_GaussianBasis_IntegralEvaluator3Interface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_MolecularInterface_hxx
#include "Chemistry_QC_GaussianBasis_MolecularInterface.hxx"
#endif
#ifndef included_MPQC_IntegralEvaluator3_hxx
#include "MPQC_IntegralEvaluator3.hxx"
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


// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3._includes)
#include "integral_evaluator.h"
#include "reorder_engine.h"
using namespace sc;
using namespace Chemistry::QC::GaussianBasis;
using namespace MpqcCca;
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3._includes)

namespace MPQC { 

  /**
   * Symbol "MPQC.IntegralEvaluator3" (version 0.2)
   */
  class IntegralEvaluator3_impl : public virtual ::MPQC::IntegralEvaluator3 
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3._inherits)

  /** IntegralEvaluator3_impl implements a class interface for
      supplying 3-center molecular integrals.

      This is an implementation of a SIDL interface.
      The stub code is generated by the Babel tool.  Do not make
      modifications outside of splicer blocks, as these will be lost.
      This is a server implementation for a Babel class, the Babel
      client code is provided by the cca-chem-generic package.
   */

  // Put additional inheritance here...
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3._inherits)
  {

  // All data marked protected will be accessable by 
  // descendant Impl classes
  protected:

    bool _wrapped;

    // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3._implementation)
 
    IntegralEvaluator< TwoBodyThreeCenterInt, 
		       twobody_threecenter_computer > eval_;
    Ref<GaussianBasisSet> bs1_, bs2_, bs3_, bs4_;
    bool reorder_;
    ReorderEngine reorder_engine_;
    twobody_threecenter_computer computer_;
    BufferSize buffer_size_;

    // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3._implementation)

  public:
    // default constructor, used for data wrapping(required)
    IntegralEvaluator3_impl();
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    IntegralEvaluator3_impl( struct MPQC_IntegralEvaluator3__object * s ) : 
      StubBase(s,true), _wrapped(false) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~IntegralEvaluator3_impl() { _dtor(); }

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
    add_evaluator_impl (
      /* in */void* eval,
      /* in */::Chemistry::QC::GaussianBasis::IntegralDescrInterface desc
    )
    ;

    /**
     * user defined non-static method.
     */
    void
    set_basis_impl (
      /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs1,
      /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs2,
      /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs3
    )
    ;

    /**
     * user defined non-static method.
     */
    void
    init_reorder_impl() ;
    /**
     * user defined non-static method.
     */
    ::Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface
    get_descriptor_impl() ;

    /**
     *  Get sidl array buffer for given type.
     * @param desc Integral descriptor.
     * @return Sidl array. 
     */
    ::sidl::array<double>
    get_array_impl (
      /* in */::Chemistry::QC::GaussianBasis::IntegralDescrInterface desc
    )
    ;


    /**
     *  Returns array of integral bounds.  When multiple integral
     * types are supported within an evaluator, the ordering
     * matches the ordering of descriptors returned by 
     * get_descriptor().
     * @return Integral bounds array. 
     */
    ::sidl::array<double>
    get_bounds_array_impl() ;

    /**
     *  This should be called when the object is no longer needed.
     * No other members may be called after finalize. 
     */
    int32_t
    finalize_impl() ;

    /**
     *  Compute a shell triplet of integrals.
     * @param shellnum1 Gaussian shell number 1.
     * @param shellnum2 Gaussian shell number 2.
     * @param shellnum3 Gaussian shell number 3.
     * @param deriv_level Derivative level. 
     */
    void
    compute_impl (
      /* in */int64_t shellnum1,
      /* in */int64_t shellnum2,
      /* in */int64_t shellnum3
    )
    ;


    /**
     *  Compute array of integral bounds.  -1 indicates a wild card
     * and the largest possible bound for given non-wild
     * card values is returned.
     * @param shellnum1 Gaussian shell number 1.
     * @param shellnum2 Gaussian shell number 2.
     * @param shellnum3 Gaussian shell number 3. 
     */
    void
    compute_bounds_impl (
      /* in */int64_t shellnum1,
      /* in */int64_t shellnum2,
      /* in */int64_t shellnum3
    )
    ;

  };  // end class IntegralEvaluator3_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator3._misc)
// Put miscellaneous things here...
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator3._misc)

#endif
