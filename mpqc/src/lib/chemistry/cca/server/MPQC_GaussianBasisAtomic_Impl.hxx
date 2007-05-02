// 
// File:          MPQC_GaussianBasisAtomic_Impl.hxx
// Symbol:        MPQC.GaussianBasisAtomic-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.GaussianBasisAtomic
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 

#ifndef included_MPQC_GaussianBasisAtomic_Impl_hxx
#define included_MPQC_GaussianBasisAtomic_Impl_hxx

#ifndef included_sidl_cxx_hxx
#include "sidl_cxx.hxx"
#endif
#ifndef included_MPQC_GaussianBasisAtomic_IOR_h
#include "MPQC_GaussianBasisAtomic_IOR.h"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_AngularType_hxx
#include "Chemistry_QC_GaussianBasis_AngularType.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_AtomicInterface_hxx
#include "Chemistry_QC_GaussianBasis_AtomicInterface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_ShellInterface_hxx
#include "Chemistry_QC_GaussianBasis_ShellInterface.hxx"
#endif
#ifndef included_MPQC_GaussianBasisAtomic_hxx
#include "MPQC_GaussianBasisAtomic.hxx"
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


// DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisAtomic._includes)

#include <chemistry/qc/basis/basis.h>
#include <MPQC_GaussianBasisShell.hxx>
using namespace std;
using namespace Chemistry::QC::GaussianBasis;
using namespace sc;
using namespace MPQC;

// DO-NOT-DELETE splicer.end(MPQC.GaussianBasisAtomic._includes)

namespace MPQC { 

  /**
   * Symbol "MPQC.GaussianBasisAtomic" (version 0.2)
   */
  class GaussianBasisAtomic_impl : public virtual ::MPQC::GaussianBasisAtomic 
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisAtomic._inherits)
  // Insert-Code-Here {MPQC.GaussianBasisAtomic._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisAtomic._inherits)
  {

  // All data marked protected will be accessable by 
  // descendant Impl classes
  protected:

    bool _wrapped;

    // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisAtomic._implementation)

    GaussianBasisSet *gbs_ptr_;
    Ref<GaussianBasisSet> sc_gbs_;
    int atomnum_;
    int nshell_;
    int max_am_;
    GaussianBasisShell *shell_array_;
    AngularType angular_type_;

    // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisAtomic._implementation)

  public:
    // default constructor, used for data wrapping(required)
    GaussianBasisAtomic_impl();
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    GaussianBasisAtomic_impl( struct MPQC_GaussianBasisAtomic__object * s ) : 
      StubBase(s,true), _wrapped(false) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~GaussianBasisAtomic_impl() { _dtor(); }

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
      /* in */void* scbasis,
      /* in */int32_t atomnum
    )
    ;


    /**
     *  Get the canonical basis set name. 
     * @return canonical name 
     */
    ::std::string
    get_name_impl() ;

    /**
     *  Get the number of basis functions.
     * @return number of functions 
     */
    int32_t
    get_n_basis_impl() ;

    /**
     *  Get the number of shells.
     * @return number of shells 
     */
    int32_t
    get_n_shell_impl() ;

    /**
     *  Get the max angular momentum for any shell on the atom.
     * @return max angular momentum 
     */
    int32_t
    get_max_angular_momentum_impl() ;

    /**
     *  Get the angular type for the atom.
     * @return enum AngularType 
     */
    ::Chemistry::QC::GaussianBasis::AngularType
    get_angular_type_impl() ;

    /**
     *  Get a gaussian shell. 
     * @param shellnum shell number
     * @return Shell 
     */
    ::Chemistry::QC::GaussianBasis::ShellInterface
    get_shell_impl (
      /* in */int32_t shellnum
    )
    ;


    /**
     *  Print the atomic basis data. 
     */
    void
    print_atomic_impl() ;
  };  // end class GaussianBasisAtomic_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisAtomic._misc)
// Insert-Code-Here {MPQC.GaussianBasisAtomic._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(MPQC.GaussianBasisAtomic._misc)

#endif
