
#ifdef __GNUC__
#pragma implementation
#endif

#include <stdlib.h>
#include <math.h>

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>

#include <chemistry/qc/basis/gaussshell.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/cartiter.h>

const char* GaussianShell::amtypes = "spdfghijkl";
const char* GaussianShell::AMTYPES = "SPDFGHIJKL";

SavableState_REF_def(GaussianShell);

#define CLASSNAME GaussianShell
#define PARENTS public SavableState
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
GaussianShell::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

int
GaussianShell::max_am() const
{
  int r = 0;
  for (int i=0; i<ncon; i++) if (r<l[i]) r = l[i];
  return r;
}

// this GaussianShell ctor allocates and computes normalization constants
// and computes nfunc
GaussianShell::GaussianShell(
  int ncn,int nprm,double*e,int*am,int*pure,double**c,PrimitiveType pt
  ):
nprim(nprm),
ncon(ncn),
l(am),
puream(pure),
exp(e),
coef(c)
{
  // Compute the number of basis functions in this shell
  compute_nfunc();

  // Convert the coefficients to coefficients for unnormalized primitives,
  // if needed.
  if (pt == Normalized) convert_coef();

  // Compute the normalization constants
  normalize_shell();
}

// this GaussianShell ctor is much like the above except the puream
// array is generated according to the value of pure
GaussianShell::GaussianShell(
  int ncn,int nprm,double*e,int*am,GaussianType pure,double**c,PrimitiveType pt
  ):
  nprim(nprm),
  ncon(ncn),
  l(am),
  exp(e),
  coef(c)
{
  puream = new int [ncontraction()];
  for (int i=0; i<ncontraction(); i++) puream[i] = (pt == Pure);

  // Compute the number of basis functions in this shell
  compute_nfunc();

  // Convert the coefficients to coefficients for unnormalized primitives,
  // if needed.
  if (pt == Normalized) convert_coef();

  // Compute the normalization constants
  normalize_shell();
}

GaussianShell::GaussianShell(const RefKeyVal&keyval)
{
  // read in the shell
  PrimitiveType pt = keyval_init(keyval,0,0);

  // Compute the number of basis functions in this shell
  compute_nfunc();

  // Convert the coefficients to coefficients for unnormalized primitives,
  // if needed.
  if (pt == Normalized) convert_coef();

  // Compute the normalization constants
  normalize_shell();
}

GaussianShell::GaussianShell(StateIn&s):
  SavableState(s)
{
  s.get(nprim);
  s.get(ncon);
  s.get(nfunc);
  s.get(l);
  s.get(puream);
  s.get(exp);
  coef = new double*[ncon];
  for (int i=0; i<ncon; i++) {
      s.get(coef[i]);
    }
}

void
GaussianShell::save_data_state(StateOut&s)
{
  s.put(nprim);
  s.put(ncon);
  s.put(nfunc);
  s.put(l,ncon);
  s.put(puream,ncon);
  s.put(exp,nprim);
  for (int i=0; i<ncon; i++) {
      s.put(coef[i],nprim);
    }
}

GaussianShell::GaussianShell(const RefKeyVal&keyval,int pure)
{
  // read in the shell
  PrimitiveType pt = keyval_init(keyval,1,pure);

  // Compute the number of basis functions in this shell
  compute_nfunc();

  // Convert the coefficients to coefficients for unnormalized primitives,
  // if needed.
  if (pt == Normalized) convert_coef();

  // Compute the normalization constants
  normalize_shell();
}

GaussianShell::PrimitiveType
GaussianShell::keyval_init(const RefKeyVal& keyval,int havepure,int pure)
{
  ncon = keyval->count("type");
  if (keyval->error() != KeyVal::OK) {
      cerr << node0 << indent
           << "GaussianShell couldn't find the \"type\" array:\n";
      keyval->dump(cerr);
      abort();
    }
  nprim = keyval->count("exp");
  if (keyval->error() != KeyVal::OK) {
      cerr << node0 << indent
           << "GaussianShell couldn't find the \"exp\" array:\n";
      keyval->dump(cerr);
      abort();
    }
  int normalized = keyval->booleanvalue("normalized");
  if (keyval->error() != KeyVal::OK) normalized = 1;
  
  l = new int[ncon];
  puream = new int[ncon];
  exp = new double[nprim];
  coef = new double*[ncon];

  int i,j;
  for (i=0; i<nprim; i++) {
      exp[i] = keyval->doublevalue("exp",i);
      if (keyval->error() != KeyVal::OK) {
          cerr << node0 << indent
               << scprintf("GaussianShell: error reading exp:%d: %s\n",
                           i,keyval->errormsg());
          keyval->errortrace(cerr);
          exit(1);
        }
    }
  for (i=0; i<ncon; i++) {
      RefKeyVal prefixkeyval = new PrefixKeyVal("type",keyval,i);
      coef[i] = new double[nprim];
      char* am = prefixkeyval->pcharvalue("am");
      if (prefixkeyval->error() != KeyVal::OK) {
          cerr << node0 << indent
               << scprintf("GaussianShell: error reading am: \"%s\"\n",
                           prefixkeyval->errormsg());
          prefixkeyval->errortrace(cerr);
          exit(1);
        }
      l[i] = -1;
      for (int li=0; amtypes[li] != '\0'; li++) {
	  if (amtypes[li] == am[0] || AMTYPES[li] == am[0]) { l[i] = li; break; }
	}
      if (l[i] == -1 || strlen(am) != 1) {
          cerr << node0 << indent
               << scprintf("GaussianShell: bad angular momentum: \"%s\"\n",
                           am);
          prefixkeyval->errortrace(cerr);
          exit(1);
	}
      if (l[i] <= 1) puream[i] = 0;
      else if (havepure) {
          puream[i] = pure;
        }
      else {
          puream[i] = prefixkeyval->booleanvalue("puream");
          if (prefixkeyval->error() != KeyVal::OK) {
              puream[i] = 0;
              //cerr << node0 << indent
              //     << scprintf("GaussianShell: error reading puream: \"%s\"\n",
              //                 prefixkeyval->errormsg());
              //exit(1);
            }
        }
      for (j=0; j<nprim; j++) {
        coef[i][j] = keyval->doublevalue("coef",i,j);
        if (keyval->error() != KeyVal::OK) {
            cerr << node0 << indent
                 << scprintf("GaussianShell: error reading coef:%d:%d: %s\n",
                             i,j,keyval->errormsg());
            keyval->errortrace(cerr);
            exit(1);
            }
        }
      delete[] am;
    }

  if (normalized) return Normalized;
  else return Unnormalized;
}

int GaussianShell::max_angular_momentum() const
{
  int max = 0;
  for (int i=0; i<ncontraction(); i++) {
      int maxi = l[i];
      if (max < maxi) max = maxi;
    }
  return max;
}

int GaussianShell::nfunction(int con) const
{
  return puream[con]?
           ((l[con]<<1)+1):
           (((l[con]+2)*(l[con]+1))>>1);
}

int GaussianShell::ncartesian(int con) const
{
  return ((l[con]+2)*(l[con]+1))>>1;
}

void GaussianShell::compute_nfunc()
{
  nfunc = 0;
  for (int i=0; i<ncontraction(); i++) nfunc += nfunction(i);
}

/* Compute the norm for ((x^x1)||(x^x2)).  This is slower than need be. */
static double
norm(int x1,int x2,double c,double ss)
{
  if (x1 < x2) return norm(x2,x1,c,ss);
  if (x1 == 1) {
    if (x2 == 1) return c * ss;
    else return 0.0;
    }
  if (x1 == 0) return ss;
  return c * ( (x1-1) * norm(x1-2,x2,c,ss) + (x2 * norm(x1-1,x2-1,c,ss)));
}

void GaussianShell::convert_coef()
{
  int i,gc;
  double c,ss;

  // Convert the contraction coefficients from coefficients over
  // normalized primitives to coefficients over unnormalized primitives
  for (gc=0; gc<ncon; gc++) {
      for (i=0; i<nprim; i++) {
	  c = 0.25/exp[i];
	  ss = pow(3.141592653589793/(exp[i]+exp[i]),1.5);
	  coef[gc][i]
	    *= 1.0/sqrt(::norm(l[gc],l[gc],c,ss));
	}
    }
}

double GaussianShell::coefficient_norm(int con,int prim) const
{
  double c = 0.25/exp[prim];
  double ss = pow(3.141592653589793/(exp[prim]+exp[prim]),1.5);
  return coef[con][prim] * sqrt(::norm(l[con],l[con],c,ss));
}

// Compute the normalization constant for a shell.
// returns 1/sqrt(<(x^l 0 0|(x^l 0 0)>).
// The formula is from Obara and Saika (for the basis functions within
// the shell that have powers of x only (a and b refer to the power
// of x):
// (a||b) = 1/(4 alpha) * ( a (a-1||b) + b (a||b-1) )
double
GaussianShell::shell_normalization(int gc)
{
  int i,j;
  double result,c,ss;

  result = 0.0;
  for (i=0; i<nprim; i++) {
    for (j=0; j<nprim; j++) {
      c = 0.50/(exp[i] + exp[j]);
      ss = pow(3.141592653589793/(exp[i]+exp[j]),1.5);
      result += coef[gc][i] * coef[gc][j] *
               ::norm(l[gc],l[gc],c,ss);
      }
    }

  return 1.0/sqrt(result);
}
 
void GaussianShell::normalize_shell()
{
  int i,gc;

  for (gc=0; gc<ncon; gc++) {
      // Normalize the contraction coefficients
      double normalization = shell_normalization(gc);
      for (i=0; i<nprim; i++) {
	  coef[gc][i] *= normalization;
	}
    }

}

static int
comp_relative_overlap(int i1, int j1, int k1, int i2, int j2, int k2)
{
  int result = 0;

  if (i1) {
      if (i1>1) result += (i1-1)*comp_relative_overlap(i1-2,j1,k1,i2,j2,k2);
      if (i2>0) result += i2*comp_relative_overlap(i1-1,j1,k1,i2-1,j2,k2);
      return result;      
    }                   
  if (j1) {             
      if (j1>1) result += (j1-1)*comp_relative_overlap(i1,j1-2,k1,i2,j2,k2);
      if (j2>0) result += j2*comp_relative_overlap(i1,j1-1,k1,i2,j2-1,k2);
      return result;
    }
  if (k1) {
      if (k1>1) result += (k1-1)*comp_relative_overlap(i1,j1,k1-2,i2,j2,k2);
      if (k2>0) result += k2*comp_relative_overlap(i1,j1,k1-1,i2,j2,k2-1);
      return result;
    }
  
  if (i2) {
      if (i2>1) result += (i2-1)*comp_relative_overlap(i1,j1,k1,i2-2,j2,k2);
      if (i1>0) result += i1*comp_relative_overlap(i1-1,j1,k1,i2-1,j2,k2);
      return result;
    }
  if (j2) {
      if (j2>1) result += (j2-1)*comp_relative_overlap(i1,j1,k1,i2,j2-2,k2);
      if (j1>0) result += j1*comp_relative_overlap(i1,j1-1,k1,i2,j2-1,k2);
      return result;
    }
  if (k2) {
      if (k2>1) result += (k2-1)*comp_relative_overlap(i1,j1,k1,i2,j2,k2-2);
      if (k1>0) result += k1*comp_relative_overlap(i1,j1,k1-1,i2,j2,k2-1);
      return result;
    }

  return 1;
}

double
GaussianShell::relative_overlap(int con,
                                int a1, int b1, int c1,
                                int a2, int b2, int c2) const
{
  int result = comp_relative_overlap(a1,b1,c1,a2,b2,c2);
  return (double) result;
}

double
GaussianShell::relative_overlap(const RefIntegral& ints,
                                int con, int func1, int func2) const
{
  if (puream[con]) {
      // depends on how intv2 currently normalizes things
      cerr << node0 << indent
           << "GaussianShell::relative_overlap "
           << "only implemented for Cartesians\n";
      abort();
    }

  CartesianIter *i1p = ints->new_cartesian_iter(l[con]);
  CartesianIter *i2p = ints->new_cartesian_iter(l[con]);

  CartesianIter& i1 = *i1p;
  CartesianIter& i2 = *i2p;

  int i;
  for (i1.start(), i=0; i<func1; i1.next(), i++);
  for (i2.start(), i=0; i<func2; i2.next(), i++);

  double ret = relative_overlap(con, i1.a(), i1.b(), i1.c(),
                                i2.a(), i2.b(), i2.c());

  delete i1p;
  delete i2p;

  return ret;
}

void
GaussianShell::print(ostream& os) const
{
  int i,j;

  os << node0 << indent << "GaussianShell:\n" << incindent
     << indent << "ncontraction = " << ncon << endl
     << indent << "nprimitive = " << nprim << endl << indent
     << "exponents:";

  for (i=0; i<nprim; i++)
      os << node0 << scprintf(" %f",exp[i]);

  os << node0 << endl << indent << "l:";
  for (i=0; i<ncon; i++)
      os << node0 << scprintf(" %d", l[i]);

  os << node0 << endl << indent << "type:";
  for (i=0; i<ncon; i++)
      os << node0 << scprintf(" %s", puream[i]?"pure":"cart");
  os << node0 << endl;

  for (i=0; i<ncon; i++) {
      os << node0 << indent << scprintf("coef[%d]:",i);
      for (j=0; j<nprim; j++)
          os << node0 << scprintf(" %f",coef[i][j]);
      os << node0 << endl;
    }

  os << node0 << decindent;
}

GaussianShell::~GaussianShell()
{
  delete[] l;
  delete[] puream;
  delete[] exp;

  for (int i=0; i<ncon; i++) {
      delete[] coef[i];
    }

  delete[] coef;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
