
/* stre.cc -- implementation of the stretch internal coordinate class
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


#define CLASSNAME StreSimpleCo
#define PARENTS public SimpleCo
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
StreSimpleCo::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SimpleCo::_castdown(cd) };
  return do_castdowns(casts,cd);
}
SimpleCo_IMPL(StreSimpleCo)

StreSimpleCo::StreSimpleCo() : SimpleCo(2) {}

StreSimpleCo::StreSimpleCo(const StreSimpleCo& s)
  : SimpleCo(2)
{
  *this=s;
  }

StreSimpleCo::StreSimpleCo(const char *re, int a1, int a2)
  : SimpleCo(2,re)
{
  atoms[0]=a1; atoms[1]=a2;
  }

StreSimpleCo::StreSimpleCo(KeyVal &kv)
  : SimpleCo(2)
{
  label_=kv.pcharvalue(0);
  atoms[0]=kv.intvalue(1);
  atoms[1]=kv.intvalue(2);
}

StreSimpleCo::StreSimpleCo(KeyVal *kv, const char *lab, int n)
  : SimpleCo(2)
{
  label_=kv->pcharvalue(lab,n,1);
  atoms[0]=kv->intvalue(lab,n,2);
  atoms[1]=kv->intvalue(lab,n,3);
  }

StreSimpleCo::~StreSimpleCo()
{
}

StreSimpleCo& StreSimpleCo::operator=(const StreSimpleCo& s)
{
  if(label_) delete[] label_;
  label_=new char[strlen(s.label_)+1]; strcpy(label_,s.label_);
  atoms[0]=s.atoms[0]; atoms[1]=s.atoms[1];

  return *this;
  }

double StreSimpleCo::calc_force_con(Molecule& m)
{
  int a=atoms[0]-1; int b=atoms[1]-1;
  double rad_ab = (m[a].element().atomic_radius()
                 + m[b].element().atomic_radius()) / 0.52917706;

  calc_intco(m);

  double k = 0.3601 * exp(-1.944*(value()-rad_ab));

  // return force constant in mdyn/ang
  return k*4.359813653/(0.52917706*0.52917706);
  }

double StreSimpleCo::calc_intco(Molecule& m, double *bmat, double coeff)
{
  int a=atoms[0]-1; int b=atoms[1]-1;
  value_ = dist(m[a].point(),m[b].point());
  if(bmat) {
    Point uu(3);
    norm(uu,m[a].point(),m[b].point());
    bmat[a*3] += coeff*uu[0]; bmat[b*3] -= coeff*uu[0];
    bmat[a*3+1] += coeff*uu[1]; bmat[b*3+1] -= coeff*uu[1];
    bmat[a*3+2] += coeff*uu[2]; bmat[b*3+2] -= coeff*uu[2];
  }

  return angstrom();
}

void StreSimpleCo::print(ostream& os, const char *pad) const
{
  os << pad << "Stretch:\n";
  if(label_) os << pad << "  ref   = " << label() << endl;
  if(atoms) os << pad << "  atoms = " << atoms[0] << " " << atoms[1] << endl;
  os << pad << "  len   = " << value() << endl;
  os.flush();
  }

void StreSimpleCo::print(FILE *of, const char *pad) const
{
  fprintf(of,"%sStretch:\n",pad);
  if(label_) fprintf(of,"%s  ref   = %s\n",pad,label());
  if(atoms) fprintf(of,"%s  atoms = %d %d\n",pad,atoms[0],atoms[1]);
  fprintf(of,"%s  len   = %lf\n",pad,value());
  fflush(of);
  }

