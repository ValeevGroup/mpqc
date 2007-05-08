//
// mcsearch.h
//
// Based on line search routines found in lbfgs.f on the WWW.
//

#ifndef _math_optimize_mcsearch_h
#define _math_optimize_mcsearch_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/optimize/opt.h>
#include <util/misc/autovec.h>

namespace sc {

/** This performs line searches with cubic steps.  It is based on the
    Fortran MCSRCH and MCSTEP routines produced by: Argonne National
    Laboratory. MINPACK Project. June 1983 Jorge J. More', David
    J. Thuente.
*/
class MCSearch: public LineOpt {
  protected:

    // These are originally from the lb3 common block.
    double gtol_, stpmin_, stpmax_;

    // Local variables in mcsrch
    double dg, fm, fx, fy, dgm, dgx, dgy, fxm, fym, stx, sty, dgxm,
	   dgym;
    int infoc;
    double finit, width, stmin, stmax;
    bool stage1;
    double width1, ftest1;
    bool brackt;
    double dginit, dgtest;

    // Local variables in mcstep
    double p, q, r__, s, sgnd, stpc, stpf, stpq, gamma, theta;
    bool bound;

    // these are saved from call to call
    int info_;
    auto_vec<double> wa_;

    void
    mcstep(double *stx, double *fx, double *dx, 
           double *sty, double *fy, double *dy, double *stp, 
           double *fp, double *dp, bool *brackt, double *stpmin, 
           double *stpmax, int *info);

    void
    mcsrch(int *n, double *x, double *f, 
           double *g, double *s, double *stp, double *ftol, 
           double *xtol, int *maxfev, int *info, int *nfev, 
           double *wa);
    

    void mcinit();
  public:

    /** The MCSearch KeyVal CTOR does not read any input.  See
        the LineOpt KeyVal CTOR for parameters that it takes.
    */
    MCSearch(const Ref<KeyVal>&);
    MCSearch();
    ~MCSearch();
    int update();

    void init(RefSCVector& direction);
    void init(RefSCVector& direction, Ref<Function> function);

    void print(std::ostream& = ExEnv::out0()) const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
