
/* lin.cc -- implementation of the linear bending internal coordinate classes
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

#define CLASSNAME LinIPSimpleCo
#define PARENTS public SimpleCo
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
LinIPSimpleCo::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SimpleCo::_castdown(cd);
  return do_castdowns(casts,cd);
}
SimpleCo_IMPL(LinIPSimpleCo)

LinIPSimpleCo::LinIPSimpleCo() : SimpleCo(4) {}

LinIPSimpleCo::LinIPSimpleCo(const LinIPSimpleCo& s)
  : SimpleCo(4)
{
  *this=s;
  }

LinIPSimpleCo::LinIPSimpleCo(const char *refr, int a1, int a2, int a3, int a4)
  : SimpleCo(4,refr)
{
  atoms[0]=a1; atoms[1]=a2; atoms[2]=a3; atoms[3]=a4;
  }

LinIPSimpleCo::~LinIPSimpleCo()
{
}

LinIPSimpleCo::LinIPSimpleCo(const RefKeyVal &kv) :
  SimpleCo(kv,4)
{
}

// LinIPSimpleCo::LinIPSimpleCo(KeyVal *kv, const char *lab, int n) :
//   SimpleCo(4)
// {
//   label_=kv->pcharvalue(lab,n,1);
//   atoms[0]=kv->intvalue(lab,n,2);
//   atoms[1]=kv->intvalue(lab,n,3);
//   atoms[2]=kv->intvalue(lab,n,4);
//   atoms[3]=kv->intvalue(lab,n,5);
//   }

LinIPSimpleCo& LinIPSimpleCo::operator=(const LinIPSimpleCo& s)
{
  if(label_) delete[] label_;
  label_=new char[strlen(s.label_)+1];
  strcpy(label_,s.label_);
  atoms[0]=s.atoms[0]; atoms[1]=s.atoms[1]; atoms[2]=s.atoms[2];
  atoms[3]=s.atoms[3];
  return *this;
  }

double LinIPSimpleCo::calc_intco(Molecule& m, double *bmat, double coeff)
{
  int a=atoms[0]-1; int b=atoms[1]-1; int c=atoms[2]-1; int d=atoms[3]-1;
  Point u1(3),u2(3),u3(3);

  norm(u1,m[a].point(),m[b].point());
  norm(u2,m[d].point(),m[b].point());
  norm(u3,m[c].point(),m[b].point());

  double co=scalar(u1,u2);
  double co2=scalar(u3,u2);

  value_ = pi-acos(co)-acos(co2);

  if (bmat) {
    double uu,ww,vv;
    Point z1(3),z2(3);
    normal(u2,u1,z2);
    normal(u1,z2,z1);
    normal(u3,u2,z2);
    normal(z2,u3,u1);
    double r1 = dist(m[a].point(),m[b].point());
    double r2 = dist(m[c].point(),m[b].point());
#if OLD_BMAT
    r1 *= bohr;
    r2 *= bohr;
#endif    
    for (int j=0; j < 3; j++) {
      uu=z1[j]/r1;
      ww=u1[j]/r2;
      vv = -uu-ww;
      bmat[a*3+j] += coeff*uu;
      bmat[b*3+j] += coeff*vv;
      bmat[c*3+j] += coeff*ww;
    }
  }

  return value_;
}

double
LinIPSimpleCo::calc_force_con(Molecule&m)
{
  return 1.0;
}

const char * LinIPSimpleCo::ctype() const
{
  return "LINIP";
}

double
LinIPSimpleCo::radians() const
{
  return value_;
}

double
LinIPSimpleCo::degrees() const
{
  return value_*rtd;
}

double
LinIPSimpleCo::preferred_value() const
{
  return value_*rtd;
}
