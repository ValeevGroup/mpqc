// 
// File:          MPQC_IntegralEvaluator2_Impl.hh
// Symbol:        MPQC.IntegralEvaluator2-v0.2
// Symbol Type:   class
// Babel Version: 0.9.8
// Description:   Server-side implementation for MPQC.IntegralEvaluator2
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.9.8
// 

#ifndef included_MPQC_IntegralEvaluator2_Impl_hh
#define included_MPQC_IntegralEvaluator2_Impl_hh

#ifndef included_sidl_cxx_hh
#include "sidl_cxx.hh"
#endif
#ifndef included_MPQC_IntegralEvaluator2_IOR_h
#include "MPQC_IntegralEvaluator2_IOR.h"
#endif
// 
// Includes for all method dependencies.
// 
#ifndef included_Chemistry_Molecule_hh
#include "Chemistry_Molecule.hh"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_Molecular_hh
#include "Chemistry_QC_GaussianBasis_Molecular.hh"
#endif
#ifndef included_MPQC_IntegralEvaluator2_hh
#include "MPQC_IntegralEvaluator2.hh"
#endif
#ifndef included_sidl_BaseInterface_hh
#include "sidl_BaseInterface.hh"
#endif
#ifndef included_sidl_ClassInfo_hh
#include "sidl_ClassInfo.hh"
#endif


// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2._includes)
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/intv3/intv3.h>
#ifdef HAVE_CINTS
  #include <chemistry/qc/cints/cints.h>
#endif
//#include <util/ref/ref.h>
using namespace sc;
using namespace Chemistry::QC::GaussianBasis;
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2._includes)

namespace MPQC { 

  /**
   * Symbol "MPQC.IntegralEvaluator2" (version 0.2)
   */
  class IntegralEvaluator2_impl
  // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2._inherits)
  // Put additional inheritance here...
  // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2._inherits)
  {

  private:
    // Pointer back to IOR.
    // Use this to dispatch back through IOR vtable.
    IntegralEvaluator2 self;

    // DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2._implementation)
    Chemistry::Molecule molecule_;
    std::string evaluator_label_;
    Ref<GaussianBasisSet> bs1_, bs2_;
    Ref<Integral> integral_;
    Ref<OneBodyInt> eval_;
    Ref<OneBodyDerivInt> deriv_eval_;
    int max_nshell2_;
    sidl::array<double> sidl_buffer_;
    const double *sc_buffer_;
    enum { one_body, one_body_deriv};
    int int_type_;
    int deriv_level_;
    std::string package_;
    // DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2._implementation)

  private:
    // private default constructor (required)
    IntegralEvaluator2_impl() {} 

  public:
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    IntegralEvaluator2_impl( struct MPQC_IntegralEvaluator2__object * s ) : 
      self(s,true) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~IntegralEvaluator2_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

  public:

    /**
     * user defined non-static method.
     */
    void
    set_integral_package (
      /*in*/ const ::std::string& label
    )
    throw () 
    ;

    /**
     * user defined non-static method.
     */
    void
    initialize_by_name (
      /*in*/ ::Chemistry::Molecule molecule,
      /*in*/ const ::std::string& basis_name,
      /*in*/ const ::std::string& evaluator_label,
      /*in*/ int64_t max_deriv
    )
    throw () 
    ;

    /**
     * user defined non-static method.
     */
    void
    initialize (
      /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs1,
      /*in*/ ::Chemistry::QC::GaussianBasis::Molecular bs2,
      /*in*/ const ::std::string& evaluator_label,
      /*in*/ int64_t max_deriv
    )
    throw () 
    ;

    /**
     * user defined non-static method.
     */
    void*
    buffer() throw () 
    ;
    /**
     * user defined non-static method.
     */
    void
    compute (
      /*in*/ int64_t shellnum1,
      /*in*/ int64_t shellnum2,
      /*in*/ int64_t deriv_level
    )
    throw () 
    ;

    /**
     * user defined non-static method.
     */
    ::sidl::array<double>
    compute_array (
      /*in*/ int64_t shellnum1,
      /*in*/ int64_t shellnum2,
      /*in*/ int64_t deriv_level
    )
    throw () 
    ;

  };  // end class IntegralEvaluator2_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.IntegralEvaluator2._misc)
// Put miscellaneous things here...
// DO-NOT-DELETE splicer.end(MPQC.IntegralEvaluator2._misc)

#endif
