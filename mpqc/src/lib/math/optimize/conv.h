
#ifndef _math_optimize_conv_h
#define _math_optimize_conv_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>
#include <util/class/class.h>
#include <math/scmat/matrix.h>
#include <math/optimize/function.h>

////////////////////////////////////////////////////////////////////////

//. The \clsnm{Convergence} class is used to check for the
//convergence of an optimization.
class Convergence: virtual_base public SavableState {
#   define CLASSNAME Convergence
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    RefSCVector grad_;
    RefSCVector x_;
    RefSCVector nextx_;
    int use_max_disp_;
    double max_disp_;
    int use_max_grad_;
    double max_grad_;
    int use_rms_disp_;
    double rms_disp_;
    int use_rms_grad_;
    double rms_grad_;
    int use_graddisp_;
    double graddisp_;

    void check_conv(const char *heading,
                    double val, double bound,
                    int &pass, int &fail);

    void set_defaults();
  public:
    //. Standard constructors and destructor.
    Convergence();
    Convergence(StateIn&);
    Convergence(const RefKeyVal&);
    virtual ~Convergence();

    void save_data_state(StateOut&);

    //. Set the current gradient and displacement.
    virtual void get_grad(const RefFunction &);
    virtual void get_x(const RefFunction &);
    virtual void get_nextx(const RefFunction &);

    //. Set the current gradient and displacement to null.
    virtual void reset();

    //. Return nonzero if the optimization has converged.
    virtual int converged();
};
SavableState_REF_dec(Convergence);

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
