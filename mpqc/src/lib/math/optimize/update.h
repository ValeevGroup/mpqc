
#ifndef _math_optimize_update_h
#define _math_optimize_update_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>
#include <util/class/class.h>
#include <math/scmat/matrix.h>
#include <math/optimize/nlp.h>

class RefNonlinearTransform;

////////////////////////////////////////////////////////////////////////
//  hessian update classes.  based on the value of inverse_hessian_
//  x and g may be reversed (see Schlegel, ab initio Methods in Quantum
//  Chemistry I, 1987, p 10

class HessianUpdate: virtual_base public SavableState {
#   define CLASSNAME HessianUpdate
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    int inverse_hessian_;
  public:
    HessianUpdate();
    HessianUpdate(StateIn&);
    HessianUpdate(const RefKeyVal&);
    void save_data_state(StateOut&);
    virtual ~HessianUpdate();
    virtual void update(RefSymmSCMatrix&hessian,RefNLP2&nlp,
                        RefSCVector&xnew,RefSCVector&gnew) = 0;
    void set_inverse();
    virtual void apply_transform(const RefNonlinearTransform&);
};
SavableState_REF_dec(HessianUpdate);

class DFPUpdate: public HessianUpdate {
#   define CLASSNAME DFPUpdate
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    RefSCVector xprev;
    RefSCVector gprev;
  public:
    DFPUpdate();
    DFPUpdate(StateIn&);
    DFPUpdate(const RefKeyVal&);
    void save_data_state(StateOut&);
    ~DFPUpdate();
    void update(RefSymmSCMatrix&ihessian,RefNLP2&nlp,
                RefSCVector&xnew,RefSCVector&gnew);
    void apply_transform(const RefNonlinearTransform&);
};

class BFGSUpdate: public DFPUpdate {
#   define CLASSNAME BFGSUpdate
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    BFGSUpdate();
    BFGSUpdate(StateIn&);
    BFGSUpdate(const RefKeyVal&);
    void save_data_state(StateOut&);
    ~BFGSUpdate();
    void update(RefSymmSCMatrix&ihessian,RefNLP2&nlp,
                RefSCVector&xnew,RefSCVector&gnew);
};

class PowellUpdate: public HessianUpdate {
#   define CLASSNAME PowellUpdate
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    RefSCVector xprev;
    RefSCVector gprev;
  public:
    PowellUpdate();
    PowellUpdate(StateIn&);
    PowellUpdate(const RefKeyVal&);
    void save_data_state(StateOut&);
    ~PowellUpdate();
    void update(RefSymmSCMatrix&ihessian,RefNLP2&nlp,
                RefSCVector&xnew,RefSCVector&gnew);
    void apply_transform(const RefNonlinearTransform&);
};

#endif
