
/* tors.cc -- implementation of the torsion internal coordinate class
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

#include "simple.h"
#include "localdef.h"
#include "chemelem.h"


#define CLASSNAME TorsSimpleCo
#define PARENTS public SimpleCo
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
TorsSimpleCo::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SimpleCo::_castdown(cd) };
  return do_castdowns(casts,cd);
}
SimpleCo_IMPL(TorsSimpleCo)


TorsSimpleCo::TorsSimpleCo() : SimpleCo(4) {}

TorsSimpleCo::TorsSimpleCo(const TorsSimpleCo& s)
  : SimpleCo(4)
{
  *this=s;
  }

TorsSimpleCo::TorsSimpleCo(const char *refr, int a1, int a2, int a3, int a4)
  : SimpleCo(4,refr)
{
  atoms[0]=a1; atoms[1]=a2; atoms[2]=a3; atoms[3]=a4;
  }

TorsSimpleCo::~TorsSimpleCo()
{
}

TorsSimpleCo::TorsSimpleCo(KeyVal &kv):
  SimpleCo(4)
{
  ref=kv.pcharvalue(0);
  atoms[0]=kv.intvalue(1);
  atoms[1]=kv.intvalue(2);
  atoms[2]=kv.intvalue(3);
  atoms[3]=kv.intvalue(4);
}

TorsSimpleCo::TorsSimpleCo(KeyVal *kv, const char *lab, int n) :
  SimpleCo(4)
{
  ref=kv->pcharvalue(lab,n,1);
  atoms[0]=kv->intvalue(lab,n,2);
  atoms[1]=kv->intvalue(lab,n,3);
  atoms[2]=kv->intvalue(lab,n,4);
  atoms[3]=kv->intvalue(lab,n,5);
  }

TorsSimpleCo& TorsSimpleCo::operator=(const TorsSimpleCo& s)
{
  if(ref) delete[] ref;
  ref=new char[strlen(s.ref)+1];
  strcpy(ref,s.ref);
  atoms[0]=s.atoms[0]; atoms[1]=s.atoms[1]; atoms[2]=s.atoms[2];
  atoms[3]=s.atoms[3];
  return *this;
  }

void TorsSimpleCo::print(ostream& os, const char *pad) const
{
  os << pad << "Torsion:\n";
  if(ref) os << pad << "  ref   = " << ref << endl;
  if(atoms) {
    os << pad << "  atoms = " << atoms[0] << " " << atoms[1];
    os << " " << atoms[2] << " " << atoms[3] << endl;
    }
  os << pad << "  omega = " << value() << endl;
  os.flush();
  }

void TorsSimpleCo::print(FILE *of, const char *pad) const
{
  fprintf(of,"%sTorsion:\n",pad);
  if(ref) fprintf(of,"%s  ref   = %s\n",pad,ref);
  if(atoms) fprintf(of,"%s  atoms = %d %d %d %d\n",pad,
     atoms[0],atoms[1],atoms[2],atoms[3]);
  fprintf(of,"%s  omega = %lf\n",pad,value());
  fflush(of);
  }

double TorsSimpleCo::calc_intco(Molecule& m, double *bmat, double coeff)
{
  int a=atoms[0]-1; int b=atoms[1]-1; int c=atoms[2]-1; int d=atoms[3]-1;
  Point u1(3),u2(3),u3(3),z1(3),z2(3);

  norm(u1,m[a].point(),m[b].point());
  norm(u2,m[c].point(),m[b].point());
  norm(u3,m[c].point(),m[d].point());

  normal(u1,u2,z1);
  normal(u3,u2,z2);

  double co=scalar(z1,z2);
  u1[0]=z1[1]*z2[2]-z1[2]*z2[1];
  u1[1]=z1[2]*z2[0]-z1[0]*z2[2];
  u1[2]=z1[0]*z2[1]-z1[1]*z2[0];
  double co2=scalar(u1,u2);

  if ((co+1.0) < 1.0e-8) co= -1.0;

  // save the old value of the torsion so we can make sure the discontinuity
  // at -pi/2 doesn't bite us

  double oldval = -value_;  

  value_=(co2<0) ? -acos(-co) : acos(-co);

  // ok, we want omega between 3*pi/2 and -pi/2, so if omega is > pi/2
  // (omega is eventually -omega), then knock 2pi off of it
  if(value_ > pih) value_ -= tpi;

  // the following tests to see if the new coordinate has crossed the
  // 3pi/2 <--> -pi/2 boundary...if so, then we add or subtract 2pi as
  // needed to prevent the transformation from internals to cartesians
  // from blowing up
  if(fabs(oldval) > 0.1) {
    while(oldval-value_ > 6.0) value_ += tpi;
    while(oldval-value_ < -6.0) value_ -= tpi;
    }

  value_ = -value_;

  if (bmat) {
    double uu,vv,ww,zz;
    norm(u1,m[a].point(),m[b].point());
    norm(u2,m[c].point(),m[b].point());
    norm(u3,m[c].point(),m[d].point());
    normal(u1,u2,z1);
    normal(u3,u2,z2);
    co=scalar(u1,u2); double si=s2(co);
    co2=scalar(u2,u3); double si2=s2(co2);
    double r1 = bohr*dist(m[a].point(),m[b].point());
    double r2 = bohr*dist(m[c].point(),m[b].point());
    double r3 = bohr*dist(m[c].point(),m[d].point());
    for (int j=0; j < 3; j++) {
      uu = z1[j]/(r1*si);
      zz = z2[j]/(r3*si2);
      vv = (r1*co/r2-1.0)*uu-zz*r3*co2/r2;
      ww = -uu-vv-zz;
      bmat[a*3+j] += coeff*uu;
      bmat[b*3+j] += coeff*vv;
      bmat[c*3+j] += coeff*ww;
      bmat[d*3+j] += coeff*zz;
    }
  }

  return value_;
}

double TorsSimpleCo::calc_force_con(Molecule& m)
{
  int a=atoms[1]-1; int b=atoms[2]-1;

  double rad_ab = (m[a].element().atomic_radius() +
                   m[b].element().atomic_radius()) / 0.52917706;

  double r_ab = dist(m[a].point(),m[b].point());

  double k = 0.0015 + 14.0*pow(1.0,0.57)/pow((rad_ab*r_ab),4.0) *
                           exp(-2.85*(r_ab-rad_ab));

  // return force constant in mdyn*ang/rad^2
  return k*4.359813653;
  }
