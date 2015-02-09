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

#include <cstdlib>
#include <cmath>
#include <vector>
#include <cassert>

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
  typeid(GaussianShell),"GaussianShell",3,"public DescribedClass",
  0, create<GaussianShell>, 0);

sc::GaussianShell::GaussianShell(const std::vector<unsigned int>& am,
                                 const std::vector<bool>& pur,
                                 const std::vector<double>& exps,
                                 const std::vector<double>& contr_coefs,
                                 PrimitiveType pt,
                                 bool norm_shell) :
                                     l(am), puream(pur),
                                     exp(exps),
                                     coef_blk(contr_coefs)
{
  // Compute the number of basis functions in this shell
    init_computed_data();

    // Convert the coefficients to coefficients for unnormalized primitives,
    // if needed.
    if (pt == Normalized) convert_coef();

    // Compute the normalization constants
    if (norm_shell) normalize_shell();
}

GaussianShell
GaussianShell::unit() {
  return GaussianShell(std::vector<unsigned int>(1,0u),
                       std::vector<bool>(1,false),
                       std::vector<double>(1,0.0),
                       std::vector<double>(1,1.0),
                       Unnormalized, false);
}

// this GaussianShell constructortor allocates and computes normalization constants
// and computes nfunc
GaussianShell::GaussianShell(
  int ncn,int nprm,double*e,int*am,int*pure,double**c,PrimitiveType pt,
  bool do_normalize_shell
  ):
l(am, am+ncn),
puream(pure, pure+ncn),
exp(e, e+nprm),
coef_blk(nprm * ncn, 0.0) // Not initalizing with pointers.
{
  for(int con=0, cp=0; con<ncn; ++con, cp+=nprm) {
    std::copy(c[con], c[con] + nprm, coef_blk.begin() + cp);
  }

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
    l(am, am+ncn),
    puream(pure == Pure, ncn),
    exp(e, e+nprm),
    coef_blk(nprm*ncn, 0.0) // Not initalizing with pointers.
{
  // Compute the number of basis functions in this shell
  init_computed_data();

  // Convert the coefficients to coefficients for unnormalized primitives,
  // if needed.
  if (pt == Normalized) convert_coef();

  // Compute the normalization constants
  normalize_shell();
}

GaussianShell::GaussianShell(const GaussianShell & other) :
    l(other.l), puream(other.puream), exp(other.exp), coef_blk(other.coef_blk)
{
  init_computed_data();
}

GaussianShell::GaussianShell(const Ref<KeyVal>&keyval, int pure)
{
  PrimitiveType pt = (pure == 0 || pure == 1)
      ? keyval_init(keyval,1,pure)
      : keyval_init(keyval,0,0);

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
  const int ncon = keyval->count("type");
  if (keyval->error() != KeyVal::OK) {
    ExEnv::err0() << indent
        << "GaussianShell couldn't find the \"type\" array:\n";
    keyval->dump(ExEnv::err0());
    throw InputError("GaussianShell KeyVal ctor failed, couldn't find \"type\" array", __FILE__, __LINE__,
                     "type");
  }
  const int nprim = keyval->count("exp");
  if (keyval->error() != KeyVal::OK) {
    ExEnv::err0() << indent
        << "GaussianShell couldn't find the \"exp\" array:\n";
    keyval->dump(ExEnv::err0());
    throw InputError("GaussianShell KeyVal ctor failed, couldn't find \"exp\" array", __FILE__, __LINE__,
                     "exp");
  }
  bool normalized = keyval->booleanvalue("normalized");
  if (keyval->error() != KeyVal::OK) normalized = 1;
  
  l.resize(ncon);
  puream.resize(ncon);
  exp.resize(nprim);
  coef_blk.resize(ncon * nprim);

  for (int i=0; i<nprim; i++) {
      exp[i] = keyval->doublevalue("exp",i);
      if (keyval->error() != KeyVal::OK) {
          std::ostringstream oss;
          oss << scprintf("GaussianShell: error reading exp:%d: %s\n",
                           i,keyval->errormsg());
          ExEnv::err0() << indent << oss.str();
          keyval->errortrace(ExEnv::err0());
          throw InputError(oss.str().c_str(), __FILE__, __LINE__);
        }
    }
  for (int i=0; i<ncon; i++) {
      Ref<KeyVal> prefixkeyval = new PrefixKeyVal(keyval,"type",i);
      std::string am = prefixkeyval->stringvalue("am");
      if (prefixkeyval->error() != KeyVal::OK) {
        std::ostringstream oss;
        oss << scprintf("GaussianShell: error reading am: \"%s\"\n",
                        prefixkeyval->errormsg());
        ExEnv::err0() << indent << oss.str();
        prefixkeyval->errortrace(ExEnv::err0());
        throw InputError(oss.str().c_str(), __FILE__, __LINE__);
        }
      bool found_am = false;
      if (am.size() > 0) {
          for (unsigned int li=0; amtypes[li] != '\0'; li++) {
              if (amtypes[li] == am[0] || AMTYPES[li] == am[0]) {
                l[i] = li;
                found_am = true;
                break;
              }
            }
        }
      if (not found_am) {
          ExEnv::err0() << indent
               << scprintf("GaussianShell: bad angular momentum: \"%s\"\n",
                           am.c_str());
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
            }
        }
      for (int j=0; j<nprim; j++) {
        coef_blk[i * nprim + j] = keyval->doublevalue("coef",i,j);
        if (keyval->error() != KeyVal::OK) {
          std::ostringstream oss;
          oss << scprintf("GaussianShell: error reading coef:%d:%d: %s\n",
                          i,j,keyval->errormsg());
          ExEnv::err0() << indent << oss.str().c_str();
          keyval->errortrace(ExEnv::err0());
          throw InputError(oss.str().c_str(), __FILE__, __LINE__);
        }
      }
    }

  if (normalized) return Normalized;
  else return Unnormalized;
}

void
GaussianShell::chomp(std::vector<double>& exp, std::vector<double>& coef_blk,
                     unsigned int ncontr, double epsilon) {
  const size_t nprim = exp.size();
  std::vector<bool> chomp_prim(nprim, false);
  bool have_prim_to_chomp = false;

  // for each primitive, make sure that there is at least 1 nonzero coefficient
  for(size_t p=0; p<nprim; ++p) {
    MPQC_ASSERT(exp[p] >= 0.0);

    bool have_nonzero_coef = false;
    for(unsigned int c=0; c<ncontr; ++c) {
      if (fabs(coef_blk[c*nprim + p]) >= epsilon) {
        have_nonzero_coef = true;
        break;
      }
    }

    if (not have_nonzero_coef) {
      chomp_prim[p] = true;
      have_prim_to_chomp = true;
    }
  }

  if (have_prim_to_chomp) {
    size_t nprim_new = 0;
    for(size_t p=0; p<nprim; ++p)
      if (chomp_prim[p] == false)
        ++nprim_new;

    std::vector<double> exp_new;
    std::vector<double> coef_blk_new(ncontr * nprim_new);
    for(size_t p=0; p<nprim; ++p) {
      if (chomp_prim[p] == false) {
        exp_new.push_back(exp[p]);
        const size_t p_new = exp_new.size() - 1;
        for(size_t c=0; c<ncontr; ++c) {
          coef_blk_new[c * nprim_new + p_new] = coef_blk[c * nprim + p];
        }
      }
    }

    exp.swap(exp_new);
    coef_blk.swap(coef_blk_new);
  }
}

void
GaussianShell::init_computed_data()
{
  chomp(exp, coef_blk, ncontraction(), epsilon());

  coef.resize(ncontraction());
  for(int c=0; c<ncontraction(); ++c) {
    coef[c] = &coef_blk[c * nprimitive()];
  }

  unsigned int max = 0;
  unsigned int min = UINT_MAX;
  unsigned int nc = 0;
  unsigned int nf = 0;
  has_pure_ = 0;
  has_cartesian_ = 0;
  contr_to_func_.resize(ncontraction());
  for (int i=0; i<ncontraction(); i++) {
      const unsigned am = l[i];
      if (max < am) max = am;
      if (min > am) min = am;

      nc += ncartesian(i);

      contr_to_func_[i] = nf;
      nf += nfunction(i);

      if (is_pure(i)) has_pure_ = 1;
      else has_cartesian_ = 1;
    }
  max_am_ = max;
  min_am_ = min;
  ncart_ = nc;
  nfunc = nf;
  
  // map functions back to contractions
  func_to_contr_.resize(nfunc);
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

unsigned int GaussianShell::ncartesian_with_aminc(int aminc) const
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
  // Convert the contraction coefficients from coefficients over
  // normalized primitives to coefficients over unnormalized primitives
  for (int gc=0; gc<ncontraction(); gc++) {
      for (int p=0; p<nprimitive(); p++) {
	  const double c = 0.25/exp[p];
	  const double ss = pow(M_PI/(exp[p]+exp[p]),1.5);
	  coef[gc][p]
	    *= 1.0/sqrt(::norm(l[gc],l[gc],c,ss));
	}
  }
}

double GaussianShell::coefficient_norm(int con,int prim) const
{
  const double c = 0.25/exp[prim];
  const double ss = pow(M_PI/(exp[prim]+exp[prim]),1.5);
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
  double result = 0.0;
  for (int i=0; i<nprimitive(); i++) {
    for (int j=0; j<nprimitive(); j++) {
      const double c = 0.50/(exp[i] + exp[j]);
      const double ss = pow(M_PI/(exp[i]+exp[j]),1.5);
      result += coef[gc][i] * coef[gc][j] *
               ::norm(l[gc],l[gc],c,ss);
      }
    }

  return 1.0/sqrt(result);
}
 
void GaussianShell::normalize_shell()
{
  for (int gc=0; gc<ncontraction(); gc++) {
    // Normalize the contraction coefficients
    const double normalization = shell_normalization(gc);
    for (int p=0; p<nprimitive(); p++) {
      coef[gc][p] *= normalization;
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
  os << indent << "GaussianShell:\n" << incindent
     << indent << "ncontraction = " << ncontraction() << endl
     << indent << "nprimitive = " << nprimitive() << endl << indent
     << "exponents:";

  for (int i=0; i<nprimitive(); i++)
      os << scprintf(" %f",exp[i]);

  os << endl << indent << "l:";
  for (int i=0; i<ncontraction(); i++)
      os << scprintf(" %d", l[i]);

  os << endl << indent << "type:";
  for (int i=0; i<ncontraction(); i++)
      os << scprintf(" %s", puream[i]?"pure":"cart");
  os << endl;

  for (int i=0; i<ncontraction(); i++) {
      os << indent << scprintf("coef[%d]:",i);
      for (int j=0; j<nprimitive(); j++)
          os << scprintf(" %f",coef[i][j]);
      os << endl;
    }

  os << decindent;
}

GaussianShell::~GaussianShell()
{
}

unsigned int
GaussianShell::nfunction(int con) const
{
  return puream[con]?
           ((l[con]<<1)+1):
           (((l[con]+2)*(l[con]+1))>>1);
}

bool
GaussianShell::equiv(const GaussianShell& s) const
{
  if (nprimitive() != s.nprimitive()) return 0;
  for (int i=0; i<nprimitive(); i++) {
      if (fabs((exp[i] - s.exp[i])/exp[i]) > epsilon()) return 0;
    }
  if (ncontraction() != s.ncontraction()) return 0;
  for (int i=0; i<ncontraction(); i++) {
      if (l[i] != s.l[i]) return 0;
      if (puream[i] != s.puream[i]) return 0;
      for (int j=0; j<nprimitive(); j++) {
        if (coef[i][j] != 0.0) { 
          if (fabs((coef[i][j] - s.coef[i][j])/coef[i][j]) > epsilon()) return 0;
        }
        else {
          if (fabs((coef[i][j] - s.coef[i][j])) > epsilon()) return 0;
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

#ifdef MPQC_NEW_FEATURES
boost::property_tree::ptree&
GaussianShell::write_xml(
    boost::property_tree::ptree& parent,
    const XMLWriter& writer
)
{
  using boost::property_tree::ptree;
  ptree& child = get_my_ptree(parent);
  child.put("nprimitive", nprimitive());
  child.put("ncontraction", ncontraction());
  child.put("nfunction", nfunction());
  child.put("has_pure", has_pure());
  {
    ptree& prims_tree = child.add_child("primitives", ptree());
    for(int iprim = 0; iprim < nprimitive(); ++iprim){
      ptree& prim_tree = prims_tree.add_child("primitive", ptree());
      prim_tree.put("<xmlattr>.index", iprim);
      prim_tree.put("exponent", exponent(iprim));
    }
  }
  {
    ptree& cons_tree = child.add_child("contractions", ptree());
    for(int icon = 0; icon < ncontraction(); ++icon){
      ptree& con_tree = cons_tree.add_child("contraction", ptree());
      con_tree.put("<xmlattr>.index", icon);
      con_tree.put("am", am(icon));
      con_tree.put("amchar", amchar(icon));
      con_tree.put("nfunction", nfunction(icon));
      con_tree.put("ncartesian", ncartesian(icon));
      con_tree.put("pure", is_pure(icon));
      con_tree.put("function_offset", contraction_to_function(icon));
      ptree& coeffs_tree = con_tree.add_child("coefficients", ptree());
      for(int iprim = 0; iprim < nprimitive(); ++iprim){
        ptree& coef_tree = coeffs_tree.add_child("coefficient", ptree());
        coef_tree.put("<xmlattr>.index", iprim);
        coef_tree.put("unnormalized", coefficient_unnorm(icon, iprim));
        coef_tree.put("normalized", coefficient_norm(icon, iprim));
      }
    }
  }
  return child;
}
#endif // MPQC_NEW_FEATURES

void
sc::ToStateOut(const GaussianShell &s, StateOut &so, int &count) {
  so.put(sc::class_desc<GaussianShell>()->version());
  ToStateOut(s.l, so, count);
  ToStateOut(s.puream, so, count);
  ToStateOut(s.exp,so, count);
  ToStateOut(s.coef_blk,so, count);
}

void
sc::FromStateIn(GaussianShell &s, StateIn &si, int &count){
  int version; si.get(version);
  if (version < 3) {
    throw FileOperationFailed("cannot restore from old GaussianShell archives",
                              __FILE__, __LINE__, 0,
                              FileOperationFailed::Corrupt,
                              s.class_desc());
  }

  FromStateIn(s.l, si, count);
  FromStateIn(s.puream, si, count);
  FromStateIn(s.exp, si, count);
  FromStateIn(s.coef_blk, si, count);

  s.init_computed_data();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
