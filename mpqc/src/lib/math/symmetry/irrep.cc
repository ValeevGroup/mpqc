
#include <stdlib.h>

#include <util/misc/newstring.h>
#include <math/symmetry/pointgrp.h>
#include <util/misc/formio.h>
#include <iomanip.h>

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
  
  ios::fmtflags oldf = os.setf(ios::left);
  
  os << indent << setw(5) << symb;

  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  
  for (i=0; i < g; i++)
    os << " " << setw(6) << setprecision(3) << character(i);
  os << " | " << ntrans_ << " t, " << nrot_ << " R\n";

  for (d=0; d < nproj(); d++) {
    os << indent << "     ";
    for (i=0; i < g; i++)
      os << " " << setw(6) << setprecision(3) << p(d,i);
    os << endl;
  }

  os.setf(oldf);
}
