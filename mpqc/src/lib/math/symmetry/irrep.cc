
/* irrep.cc -- implementation of the point group classes
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
 *      Internet: seidl@janed.com
 *      June, 1993
 */


#include <stdlib.h>

#include <util/misc/newstring.h>
#include <math/symmetry/pointgrp.h>
#include <util/misc/formio.h>

/////////////////////////////////////////////////////////////////////////

IrreducibleRepresentation::IrreducibleRepresentation() :
  g(0), degen(0), nrot_(0), ntrans_(0), complex_(0), symb(0), rep(0)
{
}

IrreducibleRepresentation::IrreducibleRepresentation(
  int order, int d, const char *lab) :
  g(0), degen(0), nrot_(0), ntrans_(0), complex_(0), symb(0), rep(0)
{
  init(order,d,lab);
}


IrreducibleRepresentation::IrreducibleRepresentation(
  const IrreducibleRepresentation& ir) :
  g(0), degen(0), nrot_(0), ntrans_(0), complex_(0), symb(0), rep(0)
{
  *this = ir;
}

IrreducibleRepresentation::~IrreducibleRepresentation()
{
  init();
}

IrreducibleRepresentation&
IrreducibleRepresentation::operator=(const IrreducibleRepresentation& ir)
{
  init(ir.g,ir.degen,ir.symb);

  nrot_ = ir.nrot_;
  ntrans_ = ir.ntrans_;
  complex_ = ir.complex_;
  
  for (int i=0; i < g; i++)
    rep[i]=ir.rep[i];
  
  return *this;
}

void
IrreducibleRepresentation::init(int order, int d, const char *lab)
{
  g=order;
  degen=d;
  ntrans_=nrot_=complex_=0;

  if (symb)
    delete[] symb;

  symb = new_string(lab);

  if (rep) {
    delete[] rep;
    rep=0;
  }

  if (g) {
    rep = new SymRep[g];
    for (int i=0; i < g; i++)
      rep[i].set_dim(d);
  }
}

void
IrreducibleRepresentation::print(ostream& os) const
{
  if (!g)
    return;

  int i,d;
  
  os << node0 << indent << scprintf("%-5s",symb);

  for (i=0; i < g; i++)
    os << node0 << scprintf(" %6.3f",character(i));
  os << node0 << " | " << ntrans_ << " t, " << nrot_ << " R\n";

  for (d=0; d < nproj(); d++) {
    os << node0 << indent << "     ";
    for (i=0; i < g; i++)
      os << node0 << scprintf(" %6.3f",p(d,i));
    os << node0 << endl;
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
// End:
