
#include <iostream.h>
#include <math/optimize/nlp.h>
#include <math/optimize/opt.h>
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>
#include <math/scmat/matrix.h>

class Quadratic: public NLP2
{
#   define CLASSNAME Quadratic
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>    
  private:
    RefSCVector x0;
    RefSCVector g0;
    RefSymmSCMatrix h0;
    RefSymmSCMatrix hguess;
  public:
    Quadratic(StateIn&);
    Quadratic(KeyVal&);
    void save_data_state(StateOut&);
    void compute();
    void guess_hessian(RefSymmSCMatrix&);
};
#define CLASSNAME Quadratic
#define PARENTS public NLP2
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
Quadratic::Quadratic(StateIn&s):
  SavableState(s,class_desc_),
  NLP2(s)
{
  x0.restore_state(s);
  g0.restore_state(s);
  h0.restore_state(s);
}
void
Quadratic::save_data_state(StateOut&s)
{
  NLP2::save_data_state(s);
  x0.save_state(s);
  g0.save_state(s);
  h0.save_state(s);
}
void *
Quadratic::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { NLP2::_castdown(cd) };
  return do_castdowns(casts,cd);
}
Quadratic::Quadratic(KeyVal&keyval):
  NLP2(new LocalSCDimension(keyval.count("x0")))
{
  x0 = dimension()->create_vector();
  g0 = dimension()->create_vector();
  h0 = dimension()->create_symmmatrix();
  hguess = dimension()->create_symmmatrix();
  hguess.assign(0.0);
  RefSCElementOp op(new SCElementShiftDiagonal(1.0));
  hguess.element_op(op);
  
  int dim = dimension()->n();
  for (int i=0; i<dim; i++) {
      x0(i) = keyval.doublevalue("x0",i);
      g0(i) = keyval.doublevalue("g0",i);
      for (int j=0; j<=i; j++) {
          h0(i,j) = keyval.doublevalue("h0",i,j);
          hguess(i,j) = keyval.doublevalue("hguess",i,j);
        }
    }
}
// this computes everything, whether or not it was requested
void
Quadratic::compute()
{
  cout << "Quadratic::compute(): entered\n";
  
  // compute the displacement from x0
  RefSCVector d = _x - x0;

  // compute h * d
  RefSCVector h0d = h0 * d;
//   RefSCVector h0d(h0.dim());
//   int n=h0.dim().n();
//   for (int i=0; i<n; i++) {
//       double tmp = 0;
//       for (int j=0; j<n; j++) {
//           tmp += h0(i,j) * d(j);
//         }
//       h0d(i) = tmp;
//     }

  // compute the value
  _value.result_noupdate() =   d.scalar_product(g0)
                             + 0.5 * d.scalar_product(h0d);
  _value.computed() = 1;

  // compute the gradient
  _gradient.result_noupdate() = g0 + h0d;
  _gradient.computed() = 1;

  // compute the hessian
  _hessian.result_noupdate() = h0;
  _hessian.computed() = 1;
}
void
Quadratic::guess_hessian(RefSymmSCMatrix&gh)
{
  gh.assign(hguess);
}

main()
{
//   RefKeyVal kv = new ParsedKeyVal( SRCDIR "/opttest.in");
//   RefKeyVal pkv = new PrefixKeyVal("opt",*kv);
  ParsedKeyVal kv( SRCDIR "/opttest.in");
  kv.unmanage();
  PrefixKeyVal pkv("opt",kv);
  pkv.unmanage();

  for (int i=0; i<pkv.count(); i++) {
      RefOptimize opt(pkv.describedclassvalue(i));
      if (opt.nonnull()) {
          RefSCVector oldx = opt->nlp()->get_x().copy();
          opt->optimize();
          // restore the orginal x, in case the function is used again
          opt->nlp()->set_x(oldx);
        }
    }
}
