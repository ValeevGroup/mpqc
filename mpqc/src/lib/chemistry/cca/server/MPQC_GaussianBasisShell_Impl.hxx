// 
// File:          MPQC_GaussianBasisShell_Impl.hxx
// Symbol:        MPQC.GaussianBasisShell-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.GaussianBasisShell
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 

#ifndef included_MPQC_GaussianBasisShell_Impl_hxx
#define included_MPQC_GaussianBasisShell_Impl_hxx

#ifndef included_sidl_cxx_hxx
#include "sidl_cxx.hxx"
#endif
#ifndef included_MPQC_GaussianBasisShell_IOR_h
#include "MPQC_GaussianBasisShell_IOR.h"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_AngularType_hxx
#include "Chemistry_QC_GaussianBasis_AngularType.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_ShellInterface_hxx
#include "Chemistry_QC_GaussianBasis_ShellInterface.hxx"
#endif
#ifndef included_MPQC_GaussianBasisShell_hxx
#include "MPQC_GaussianBasisShell.hxx"
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


// DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisShell._includes)

#include <chemistry/qc/basis/basis.h>
using namespace std;
using namespace Chemistry::QC::GaussianBasis;
using namespace sc;

// DO-NOT-DELETE splicer.end(MPQC.GaussianBasisShell._includes)

namespace MPQC { 

  /**
   * Symbol "MPQC.GaussianBasisShell" (version 0.2)
   */
  class GaussianBasisShell_impl : public virtual ::MPQC::GaussianBasisShell 
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisShell._inherits)
  // Insert-Code-Here {MPQC.GaussianBasisShell._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisShell._inherits)
  {

  // All data marked protected will be accessable by 
  // descendant Impl classes
  protected:

    bool _wrapped;

    // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisShell._implementation)

    GaussianShell *shell_ptr_;
    Ref<GaussianShell> sc_shell_;
    AngularType angular_type_;
    int max_am_;

    // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisShell._implementation)

  public:
    // default constructor, used for data wrapping(required)
    GaussianBasisShell_impl();
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    GaussianBasisShell_impl( struct MPQC_GaussianBasisShell__object * s ) : 
      StubBase(s,true), _wrapped(false) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~GaussianBasisShell_impl() { _dtor(); }

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
    initialize_impl (
      /* in */void* scshell
    )
    ;


    /**
     *  Get the number of contractions in the shell. 
     * @return number of contractions 
     */
    int32_t
    get_n_contraction_impl() ;

    /**
     *  Get the number of primitives in the shell.
     * @return number of primitives 
     */
    int32_t
    get_n_primitive_impl() ;

    /**
     *  Get the coefficient for an unnormalized primitive 
     * in a contraction.
     * @param connum contraction number
     * @param expnum primitive number
     * @return contraction coefficient 
     */
    double
    get_contraction_coef_impl (
      /* in */int32_t connum,
      /* in */int32_t expnum
    )
    ;


    /**
     *  Get the exponent for a primitive.
     * @param expnum primitive id number
     * @return exponent 
     */
    double
    get_exponent_impl (
      /* in */int32_t expnum
    )
    ;


    /**
     *  Get the angular momentum for a single contraction.
     * @param connum contraction id number
     * @return angular momentum value 
     */
    int32_t
    get_angular_momentum_impl (
      /* in */int32_t connum
    )
    ;


    /**
     *  Get the max angular momentum, considering all contractions 
     * in the shell.
     * @return maximum angular momentum value 
     */
    int32_t
    get_max_angular_momentum_impl() ;

    /**
     *  Get the angular type for a single contraction.
     * @param connum contraction number
     * @return enum AngularType 
     */
    ::Chemistry::QC::GaussianBasis::AngularType
    get_contraction_angular_type_impl (
      /* in */int32_t connum
    )
    ;


    /**
     *  Get the angular type.
     * @return enum AngularType 
     */
    ::Chemistry::QC::GaussianBasis::AngularType
    get_angular_type_impl() ;

    /**
     *  Print the shell data. 
     */
    void
    print_shell_impl() ;
  };  // end class GaussianBasisShell_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisShell._misc)
// Insert-Code-Here {MPQC.GaussianBasisShell._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(MPQC.GaussianBasisShell._misc)

#endif
