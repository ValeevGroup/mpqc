
#include "nlp.h"
#include "opt.h"

#include <util/keyval/keyval.h>
#include <math/nihmatrix/nihmatrix.h>
#include <math/nihmatrix/newmat.h>

// force NoUpdate to be linked in
ClassDesc* NoUpdatelinkage = &NoUpdate::class_desc_;

class Quadratic: public NLP2
{
#   define CLASSNAME Quadratic
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>    
  private:
    int dim;
    DVector x0;
    DVector g0;
    DMatrix h0;
  public:
    Quadratic(StateIn&);
    Quadratic(KeyVal&);
    void save_data_state(StateOut&);
    void compute();
};
#   define CLASSNAME Quadratic
#   define PARENTS public NLP2
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/statei.h>
#   include <util/class/classi.h>
Quadratic::Quadratic(StateIn&s):
  SavableState(s,class_desc_),
  NLP2(s),
  x0(s),
  g0(s),
  h0(s)
{
}
void
Quadratic::save_data_state(StateOut&s)
{
  NLP2::save_data_state(s);
  x0.save_object_state(s);
  g0.save_object_state(s);
  h0.save_object_state(s);
}
void *
Quadratic::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { NLP2::_castdown(cd) };
  return do_castdowns(casts,cd);
}
Quadratic::Quadratic(KeyVal&kv):
  NLP2(kv),
  x0(PrefixKeyVal("x0",kv)),
  g0(PrefixKeyVal("g0",kv))
{
  dim = x0.dim();
  
  // symmetrize h0
  DMatrix h0old(PrefixKeyVal("h0",kv));
  h0.resize(dim,dim);
  for (int i=0; i<dim; i++) {
      for (int j=0; j<dim; j++) {
          h0(i,j) = 0.5*(h0old(i,j) + h0old(j,i));
        }
    }

  grad.result_noupdate().ReDimension(dim);
  Hessian.result_noupdate().ReDimension(dim);
}
// this computes everything, whether or not it was requested
void
Quadratic::compute()
{
  int i;
  
  // compute the displacement from x0
  DVector d(dim);
  for (i=0; i<dim; i++) {
      d[i] = xc.element(i) - x0[i];
    }

  // compute h * d
  DVector h0d(dim);
  h0d.zero();
  for (i=0; i<dim; i++) {
      for (int j=0; j<dim; j++) {
          h0d[i] += h0(i,j)*d[j];
        }
    }

//   printf("h0:\n");
//   h0.print();
// 
//   printf("h0d:\n");
//   h0d.print();

  // compute the value
  fvalue.result_noupdate() = d.dot(g0) + 0.5 * d.dot(h0d);
  fvalue.computed() = 1;

  // compute the gradient
  DVector g = g0 + h0d;
  Convert(g,grad.result_noupdate());
  grad.computed() = 1;

  // compute the hessian
  Convert(h0,Hessian.result_noupdate());
  Hessian.computed() = 1;

//   printf("displacement:\n");
//   d.print();
//   printf("grad:\n");
//   DVector tmp;
//   Convert(grad.result(),tmp);
//   tmp.print();
//   printf("fvalue = %15.10f\n",fvalue.result());
}

main()
{
  ParsedKeyVal kv("opttest.in");
  PrefixKeyVal pkv("opt",kv);

  for (int i=0; i<pkv.count(); i++) {
      RefOptimize opt(pkv.describedclassvalue(i));
      if (opt.nonnull()) opt->optimize();
    }

}
