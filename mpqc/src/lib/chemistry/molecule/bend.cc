
/* bend.cc -- implementation of the bending simple internal coordinate class
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
 *      February, 1993
 */

#include <string.h>
#include <math.h>

#include <chemistry/molecule/simple.h>
#include <chemistry/molecule/localdef.h>
#include <chemistry/molecule/chemelem.h>

#define CLASSNAME BendSimpleCo
#define PARENTS public SimpleCo
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
BendSimpleCo::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SimpleCo::_castdown(cd);
  return do_castdowns(casts,cd);
}
SimpleCo_IMPL(BendSimpleCo)

BendSimpleCo::BendSimpleCo() : SimpleCo(3) {}

BendSimpleCo::BendSimpleCo(const BendSimpleCo& s)
  : SimpleCo(3)
{
  *this=s;
  }

BendSimpleCo::BendSimpleCo(const char *refr, int a1, int a2, int a3)
  : SimpleCo(3,refr)
{
  atoms[0]=a1; atoms[1]=a2; atoms[2]=a3;
  }

BendSimpleCo::BendSimpleCo(const RefKeyVal &kv)
  : SimpleCo(kv,3)
{
}

BendSimpleCo::~BendSimpleCo()
{
}

BendSimpleCo& BendSimpleCo::operator=(const BendSimpleCo& s)
{
  if(label_) delete[] label_;
  label_=new char[strlen(s.label_)+1]; strcpy(label_,s.label_);
  atoms[0]=s.atoms[0]; atoms[1]=s.atoms[1]; atoms[2]=s.atoms[2];
  return *this;
  }

double BendSimpleCo::calc_intco(Molecule& m, double *bmat, double coeff)
{
  Point u1(3),u2(3);
  int a=atoms[0]-1; int b=atoms[1]-1; int c=atoms[2]-1;

  norm(u1,m[a].point(),m[b].point());
  norm(u2,m[c].point(),m[b].point());

  double co=scalar(u1,u2);

  value_=acos(co);

  if(bmat) {
    double uu,ww,vv;
    double si=s2(co);
    double r1i = 1.0/(si*dist(m[a].point(),m[b].point()));
    double r2i = 1.0/(si*dist(m[c].point(),m[b].point()));
#if OLD_BMAT
    r1i /= bohr;
    r2i /= bohr;
#endif    
    for (int j=0; j < 3; j++) {
      uu = (co*u1[j]-u2[j])*r1i;
      ww = (co*u2[j]-u1[j])*r2i;
      vv = -uu-ww;
      bmat[a*3+j] += coeff*uu;
      bmat[b*3+j] += coeff*vv;
      bmat[c*3+j] += coeff*ww;
    }
  }

  return value_;
}

double BendSimpleCo::calc_force_con(Molecule& m)
{
  int a=atoms[1]-1; int b=atoms[0]-1; int c=atoms[2]-1;

  double rad_ab =   m[a].element().atomic_radius()
                  + m[b].element().atomic_radius();

  double rad_ac =   m[a].element().atomic_radius()
                  + m[c].element().atomic_radius();

  double r_ab = dist(m[a].point(),m[b].point());
  double r_ac = dist(m[a].point(),m[c].point());

  double k = 0.089 + 0.11/pow((rad_ab*rad_ac),-0.42) *
                           exp(-0.44*(r_ab+r_ac-rad_ab-rad_ac));

#if OLD_BMAT
  // return force constant in mdyn*ang/rad^2
  return k*4.359813653;
#else  
  return k;
#endif  
}

const char *
BendSimpleCo::ctype() const
{
  return "BEND";
}

double
BendSimpleCo::radians() const
{
  return value_;
}

double
BendSimpleCo::degrees() const
{
  return value_*rtd;
}

double
BendSimpleCo::preferred_value() const
{
  return value_*rtd;
}
