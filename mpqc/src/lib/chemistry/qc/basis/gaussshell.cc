//
// gaussshell.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdlib.h>
#include <math.h>

#include <util/misc/formio.h>
#include <util/misc/math.h>
#include <util/keyval/keyval.h>
#include <util/state/stateio.h>

#include <chemistry/qc/basis/gaussshell.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/cartiter.h>

using namespace std;
using namespace sc;

const char* GaussianShell::amtypes = "spdfghiklmnoqrtuvwxyz";
const char* GaussianShell::AMTYPES = "SPDFGHIKLMNOQRTUVWXYZ";

static ClassDesc GaussianShell_cd(
  typeid(GaussianShell),"GaussianShell",2,"public SavableState",
  0, create<GaussianShell>, create<GaussianShell>);

// this GaussianShell ctor allocates and computes normalization constants
// and computes nfunc
GaussianShell::GaussianShell(
  int ncn,int nprm,double*e,int*am,int*pure,double**c,PrimitiveType pt,
  bool do_normalize_shell
  ):
nprim(nprm),
ncon(ncn),
l(am),
puream(pure),
exp(e),
coef(c)
{
  // Compute the number of basis functions in this shell
  init_computed_data();

  // Convert the coefficients to coefficients for unnormalized primitives,
  // if needed.
  if (pt == Normalized) convert_coef();

  // Compute the normalization constants
  if (do_normalize_shell) normalize_shell();
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
  for (int i=0; i<ncontraction(); i++) puream[i] = (pure == Pure);

  // Compute the number of basis functions in this shell
  init_computed_data();

  // Convert the coefficients to coefficients for unnormalized primitives,
  // if needed.
  if (pt == Normalized) convert_coef();

  // Compute the normalization constants
  normalize_shell();
}

GaussianShell::GaussianShell(const Ref<KeyVal>&keyval)
{
  // read in the shell
  PrimitiveType pt = keyval_init(keyval,0,0);

  // Compute the number of basis functions in this shell
  init_computed_data();

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
  if (s.version(::class_desc<GaussianShell>()) < 2) s.get(nfunc);
  s.get(l);
  s.get(puream);
  s.get(exp);
  coef = new double*[ncon];
  for (int i=0; i<ncon; i++) {
      s.get(coef[i]);
    }
  init_computed_data();
}

void
GaussianShell::save_data_state(StateOut&s)
{
  s.put(nprim);
  s.put(ncon);
  s.put(l,ncon);
  s.put(puream,ncon);
  s.put(exp,nprim);
  for (int i=0; i<ncon; i++) {
      s.put(coef[i],nprim);
    }
}

GaussianShell::GaussianShell(const Ref<KeyVal>&keyval,int pure)
{
  // read in the shell
  PrimitiveType pt = keyval_init(keyval,1,pure);

  // Compute the number of basis functions in this shell
  init_computed_data();

  // Convert the coefficients to coefficients for unnormalized primitives,
  // if needed.
  if (pt == Normalized) convert_coef();

  // Compute the normalization constants
  normalize_shell();
}

GaussianShell::PrimitiveType
GaussianShell::keyval_init(const Ref<KeyVal>& keyval,int havepure,int pure)
{
  ncon = keyval->count("type");
  if (keyval->error() != KeyVal::OK) {
      ExEnv::err0() << indent
           << "GaussianShell couldn't find the \"type\" array:\n";
      keyval->dump(ExEnv::err0());
      abort();
    }
  nprim = keyval->count("exp");
  if (keyval->error() != KeyVal::OK) {
      ExEnv::err0() << indent
           << "GaussianShell couldn't find the \"exp\" array:\n";
      keyval->dump(ExEnv::err0());
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
          ExEnv::err0() << indent
               << scprintf("GaussianShell: error reading exp:%d: %s\n",
                           i,keyval->errormsg());
          keyval->errortrace(ExEnv::err0());
          exit(1);
        }
    }
  for (i=0; i<ncon; i++) {
      Ref<KeyVal> prefixkeyval = new PrefixKeyVal(keyval,"type",i);
      coef[i] = new double[nprim];
      char* am = prefixkeyval->pcharvalue("am");
      if (prefixkeyval->error() != KeyVal::OK) {
          ExEnv::err0() << indent
               << scprintf("GaussianShell: error reading am: \"%s\"\n",
                           prefixkeyval->errormsg());
          prefixkeyval->errortrace(ExEnv::err0());
          exit(1);
        }
      l[i] = -1;
      for (int li=0; amtypes[li] != '\0'; li++) {
	  if (amtypes[li] == am[0] || AMTYPES[li] == am[0]) { l[i] = li; break; }
	}
      if (l[i] == -1 || strlen(am) != 1) {
          ExEnv::err0() << indent
               << scprintf("GaussianShell: bad angular momentum: \"%s\"\n",
                           am);
          prefixkeyval->errortrace(ExEnv::err0());
          exit(1);
	}
      if (l[i] < 1) puream[i] = 0;
      else if (havepure) {
          puream[i] = pure;
        }
      else {
          puream[i] = prefixkeyval->booleanvalue("puream");
          if (prefixkeyval->error() != KeyVal::OK) {
              puream[i] = 0;
              //ExEnv::err0() << indent
              //     << scprintf("GaussianShell: error reading puream: \"%s\"\n",
              //                 prefixkeyval->errormsg());
              //exit(1);
            }
        }
      for (j=0; j<nprim; j++) {
        coef[i][j] = keyval->doublevalue("coef",i,j);
        if (keyval->error() != KeyVal::OK) {
            ExEnv::err0() << indent
                 << scprintf("GaussianShell: error reading coef:%d:%d: %s\n",
                             i,j,keyval->errormsg());
            keyval->errortrace(ExEnv::err0());
            exit(1);
            }
        }
      delete[] am;
    }

  if (normalized) return Normalized;
  else return Unnormalized;
}

void
GaussianShell::init_computed_data()
{
  int max = 0;
  int min = 0;
  int nc = 0;
  int nf = 0;
  has_pure_ = 0;
  contr_to_func_ = new int[ncontraction()];
  for (int i=0; i<ncontraction(); i++) {
      int maxi = l[i];
      if (max < maxi) max = maxi;

      int mini = l[i];
      if (min > mini || i == 0) min = mini;

      nc += ncartesian(i);

      contr_to_func_[i] = nf;
      nf += nfunction(i);

      if (is_pure(i)) has_pure_ = 1;
    }
  max_am_ = max;
  min_am_ = min;
  ncart_ = nc;
  nfunc = nf;
  
  // map functions back to contractions
  func_to_contr_ = new int[nfunc];
  int f_offset=0;
  for (int i=0; i<ncontraction(); i++) {
    const int nf = nfunction(i);
    for(int f=0; f<nf; ++f)
      func_to_contr_[f+f_offset] = i;
    f_offset += nf;
  }
}

int GaussianShell::max_cartesian() const
{
  int max = 0;
  for (int i=0; i<ncontraction(); i++) {
      int maxi = ncartesian(i);
      if (max < maxi) max = maxi;
    }
  return max;
}

int GaussianShell::ncartesian_with_aminc(int aminc) const
{
  int ret = 0;
  for (int i=0; i<ncontraction(); i++) {
      ret += (((l[i]+2+aminc)*(l[i]+1+aminc))>>1);
    }
  return ret;
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
	  ss = pow(M_PI/(exp[i]+exp[i]),1.5);
	  coef[gc][i]
	    *= 1.0/sqrt(::norm(l[gc],l[gc],c,ss));
	}
    }
}

double GaussianShell::coefficient_norm(int con,int prim) const
{
  double c = 0.25/exp[prim];
  double ss = pow(M_PI/(exp[prim]+exp[prim]),1.5);
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
      ss = pow(M_PI/(exp[i]+exp[j]),1.5);
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
GaussianShell::relative_overlap(const Ref<Integral>& ints,
                                int con, int func1, int func2) const
{
  if (puream[con]) {
      // depends on how intv2 currently normalizes things
      ExEnv::err0() << indent
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

  os << indent << "GaussianShell:\n" << incindent
     << indent << "ncontraction = " << ncon << endl
     << indent << "nprimitive = " << nprim << endl << indent
     << "exponents:";

  for (i=0; i<nprim; i++)
      os << scprintf(" %f",exp[i]);

  os << endl << indent << "l:";
  for (i=0; i<ncon; i++)
      os << scprintf(" %d", l[i]);

  os << endl << indent << "type:";
  for (i=0; i<ncon; i++)
      os << scprintf(" %s", puream[i]?"pure":"cart");
  os << endl;

  for (i=0; i<ncon; i++) {
      os << indent << scprintf("coef[%d]:",i);
      for (j=0; j<nprim; j++)
          os << scprintf(" %f",coef[i][j]);
      os << endl;
    }

  os << decindent;
}

GaussianShell::~GaussianShell()
{
  delete[] l;
  delete[] puream;
  delete[] exp;
  delete[] contr_to_func_;
  delete[] func_to_contr_;

  for (int i=0; i<ncon; i++) {
      delete[] coef[i];
    }

  delete[] coef;
}

int
GaussianShell::nfunction(int con) const
{
  return puream[con]?
           ((l[con]<<1)+1):
           (((l[con]+2)*(l[con]+1))>>1);
}

int
GaussianShell::equiv(const GaussianShell *s)
{
  if (nprim != s->nprim) return 0;
  if (ncon != s->ncon) return 0;
  for (int i=0; i<ncon; i++) {
      if (l[i] != s->l[i]) return 0;
      if (puream[i] != s->puream[i]) return 0;
      if (fabs((exp[i] - s->exp[i])/exp[i]) > 1.0e-13) return 0;
      for (int j=0; j<nprim; j++) {
        if (coef[i][j] != 0.0) { 
          if (fabs((coef[i][j] - s->coef[i][j])/coef[i][j]) > 1.0e-13) return 0;
        }
        else {
          if (fabs((coef[i][j] - s->coef[i][j])) > 1.0e-13) return 0;
        }
      }
    }
  return 1;
}

double
GaussianShell::extent(double threshold) const
{
  double tol = 0.1;
  double r0 = tol;
  double r1 = 3.0*r0;
  double b0 = monobound(r0);
  double b1 = monobound(r1);
  //ExEnv::outn() << "r0 = " << r0 << " b0 = " << b0 << endl;
  //ExEnv::outn() << "r1 = " << r0 << " b1 = " << b1 << endl;
  if (b0 <= threshold) {
      return r0;
    }
  // step out until r0 and r1 bracket the return value
  while (b1 > threshold) {
      r0 = r1;
      r1 = 3.0*r0;
      b0 = b1;
      b1 = monobound(r1);
      //ExEnv::outn() << "r0 = " << r0 << " b0 = " << b0 << endl;
      //ExEnv::outn() << "r1 = " << r0 << " b1 = " << b1 << endl;
    }
  while (r1 - r0 > 0.1) {
      double rtest = 0.5*(r0+r1);
      double btest = monobound(rtest);
      if (btest <= threshold) {
          b1 = btest;
          r1 = rtest;
          //ExEnv::outn() << "r1 = " << r0 << " b1 = " << b0 << endl;
        }
      else {
          b0 = btest;
          r0 = rtest;
          //ExEnv::outn() << "r0 = " << r0 << " b0 = " << b0 << endl;
        }
    }
  return r1;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
