
/* mpqcstate.h -- definition of mpqc state class
 *
 *      THIS SOFTWARE FITS THE DESCRIPTION IN THE U.S. COPYRIGHT ACT OF A
 *      "UNITED STATES GOVERNMENT WORK".  IT WAS WRITTEN AS A PART OF THE
 *      AUTHOR'S OFFICIAL DUTIES AS A GOVERNMENT EMPLOYEE.  THIS MEANS IT
 *      CANNOT BE COPYRIGHTED.  THIS SOFTWARE IS FREELY AVAILABLE TO THE
 *      PUBLIC FOR USE WITHOUT A COPYRIGHT NOTICE, AND THERE ARE NO
 *      RESTRICTIONS ON ITS USE, NOW OR SUBSEQUENTLY.
 *
 *  Author:
 *      E. T. Seidl
 *      Bldg. 12A, Rm. 2033
 *      Computer Systems Laboratory
 *      Division of Computer Research and Technology
 *      National Institutes of Health
 *      Bethesda, Maryland 20892
 *      Internet: seidl@alw.nih.gov
 *      November, 1993
 */

#ifndef _mpqcstate_h
#define _mpqcstate_h

#include <math.h>

#include <util/class/class.h>
#include <util/state/state.h>
#include <math/nihmatrix/nihmatrix.h>
#include <chemistry/molecule/molecule.h>

#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>

#include "input.h"

class MPQC_state : virtual public SavableState {
#   define CLASSNAME MPQC_state
#   define HAVE_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>

  public:
    enum mstate { UNK, NEED_WFN, NEED_GRADIENT, NEED_HESSIAN, NEED_NEXT_GEOM,
                  NEED_PROPER, NEED_SCF_WFN, DONE, ABORT };

  private:
    enum mstate state_;
    int wfn_converged_;
    int geom_converged_;

    DMatrix hessian_;
    DVector coord_;
    DVector deltax_;
    DVector force_;
    Molecule mol_;

    int nsave;
    DVector *coords_;
    DVector *forces_;
    DVector *error_;

    static centers_t *centers;

    void read_input(KeyVal&,MPQC_input&,MPQC_geom_input&,centers_t&);
    int get_symmco(KeyVal&,MPQC_input&,MPQC_geom_input&);
  public:
    MPQC_state();
    MPQC_state(StateIn&);
    MPQC_state(KeyVal&,MPQC_input&,MPQC_geom_input&,centers_t&);
    MPQC_state(const char*,MPQC_input&,MPQC_geom_input&,centers_t&);

    void save_data_state(StateOut&);

    inline enum mstate state() const { return state_; }
    inline int wfn_converged() const { return wfn_converged_; }
    inline int geom_converged() const { return geom_converged_; }

    inline DMatrix& hessian() { return hessian_; }
    inline DVector& coord() { return coord_; }
    inline DVector& deltax() { return deltax_; }
    inline DVector& force() { return force_; }
    inline Molecule& molecule() { return mol_; }

    inline DVector *coords() { return coords_; }
    inline DVector *forces() { return forces_; }
    inline DVector *error() { return error_; }

    inline void converged_wfn() { wfn_converged_=1; }
    inline void converged_geom() { geom_converged_=1; }
    inline void reset() { wfn_converged_=geom_converged_=0; }

    void did_scf(MPQC_input&);
    void did_wfn(MPQC_input&);
    void did_gradient(MPQC_input&);
    void did_hessian(MPQC_input&);
    void did_next_geom(MPQC_input&);
    void did_proper(MPQC_input&);

    void broadcast();
  };

#endif
