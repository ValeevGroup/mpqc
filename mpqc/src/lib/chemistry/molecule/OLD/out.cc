
/* out.cc -- implementation of the out-of-plane internal coordinate class
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


DescribedClass_IMPL(OutSimpleCo,1,"SimpleCo","")
SavableState_IMPL(OutSimpleCo)
SimpleCo_IMPL(OutSimpleCo)
void * OutSimpleCo::_castdown(const ClassDesc *cd)
{
  if(&class_desc_ == cd) return this;
  return SimpleCo::_castdown(cd);
  }

OutSimpleCo::OutSimpleCo() : SimpleCo(4) {}

OutSimpleCo::OutSimpleCo(const OutSimpleCo& s)
  : SimpleCo(4)
{
  *this=s;
  }

OutSimpleCo::OutSimpleCo(const char *refr, int a1, int a2, int a3, int a4)
  : SimpleCo(4,refr)
{
  atoms[0]=a1; atoms[1]=a2; atoms[2]=a3; atoms[3]=a4;
  }

OutSimpleCo::OutSimpleCo(KeyVal *kv, const char *lab, int n) :
  SimpleCo(4)
{
  ref=kv->pcharvalue(lab,n,1);
  atoms[0]=kv->intvalue(lab,n,2);
  atoms[1]=kv->intvalue(lab,n,3);
  atoms[2]=kv->intvalue(lab,n,4);
  atoms[3]=kv->intvalue(lab,n,5);
  }

OutSimpleCo::~OutSimpleCo()
{
}

OutSimpleCo& OutSimpleCo::operator=(const OutSimpleCo& s)
{
  if(ref) delete[] ref;
  ref=new char[strlen(s.ref)+1];
  strcpy(ref,s.ref);
  atoms[0]=s.atoms[0]; atoms[1]=s.atoms[1]; atoms[2]=s.atoms[2];
  atoms[3]=s.atoms[3];
  return *this;
  }

void OutSimpleCo::print(ostream& os, const char *pad) const
{
  os << pad << "Out-of-plane:\n";
  if(ref) os << pad << "  ref   = " << ref << endl;
  if(atoms) {
    os << pad << "  atoms = " << atoms[0] << " " << atoms[1];
    os << " " << atoms[2] << " " << atoms[3] << endl;
    }
  os << pad << "  theta = " << value() << endl;
  os.flush();
  }

void OutSimpleCo::print(FILE *of, const char *pad) const
{
  fprintf(of,"%sOut-of-plane:\n",pad);
  if(ref) fprintf(of,"%s  ref   = %s\n",pad,ref);
  if(atoms) fprintf(of,"%s  atoms = %d %d %d %d\n",pad,
        atoms[0],atoms[1],atoms[2],atoms[3]);
  fprintf(of,"%s  theta = %lf\n",pad,value());
  fflush(of);
  }

double OutSimpleCo::calc_intco(Molecule& m, double *bmat, double coeff)
{
  int a=atoms[0]-1; int b=atoms[1]-1; int c=atoms[2]-1; int d=atoms[3]-1;
  Point u1(3),u2(3),u3(3),z1(3);

  norm(u1,m[a].point(),m[b].point());
  norm(u2,m[c].point(),m[b].point());
  norm(u3,m[d].point(),m[b].point());

  normal(u2,u3,z1);
  double st=scalar(u1,z1);
  double ct=s2(st);

  value_ = (st<0) ? -acos(ct) : acos(ct);

  if (bmat) {
    double uu,vv;
    Point ww(3),xx(3),zz(3);
    double cphi1 = scalar(u2,u3);
    double sphi1 = s2(cphi1);
    double cphi2 = scalar(u3,u1);
    double cphi3 = scalar(u2,u1);
    double den = ct * sphi1*sphi1;
    double sthta2 = (cphi1*cphi2-cphi3)/
              (den*bohr*dist(m[c].point(),m[b].point()));
    double sthta3 = (cphi1*cphi3-cphi2)/
              (den*bohr*dist(m[d].point(),m[b].point()));
    for(int j=0; j < 3; j++) {
      ww[j] = z1[j]*sthta2;
      zz[j] = z1[j]*sthta3;
    }
    normal(z1,u1,xx);
    normal(u1,xx,z1);
    double r1i = 1.0/(bohr*dist(m[a].point(),m[b].point()));
    for(j=0; j < 3; j++) {
      uu = z1[j]*r1i;
      vv = -uu-ww[j]-zz[j];
      bmat[a*3+j] += coeff*uu;
      bmat[b*3+j] += coeff*vv;
      bmat[c*3+j] += coeff*ww[j];
      bmat[d*3+j] += coeff*zz[j];
    }
  }

  return value_;
}


double OutSimpleCo::calc_force_con(Molecule& m)
{
  int x=atoms[0]-1;
  int a=atoms[1]-1; int b=atoms[2]-1; int c=atoms[3]-1;

  double rad_ab = (m[a].element().atomic_radius()
                +  m[b].element().atomic_radius()) / 0.52917706;

  double rad_ac = (m[a].element().atomic_radius()
                +  m[c].element().atomic_radius()) / 0.52917706;

  double rad_ax = (m[a].element().atomic_radius()
                +  m[x].element().atomic_radius()) / 0.52917706;

  double r_ax = dist(m[a].point(),m[x].point());

  calc_intco(m);

  double k = 0.0025 + 0.0061*pow((rad_ab*rad_ac),0.80)*pow(cos(value()),4.0) *
                           exp(-3.0*(r_ax-rad_ax));

  // return force constant in mdyn*ang/rad^2
  return k*4.359813653;
  }
