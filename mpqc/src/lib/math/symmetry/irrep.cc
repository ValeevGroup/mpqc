
#include <stdlib.h>

#include <util/misc/newstring.h>
#include <math/symmetry/pointgrp.h>

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
IrreducibleRepresentation::print(FILE *fp, const char *off)
{
  if (!g)
    return;

  int i,d;
  
  fprintf(fp,"%s%-5s",off,symb);
  for (i=0; i < g; i++)
    fprintf(fp," %6.3f",character(i));
  fprintf(fp," | %d t, %d R\n",ntrans_,nrot_);

  for (d=0; d < nproj(); d++) {
    fprintf(fp,"%s     ",off);
    for (i=0; i < g; i++)
      fprintf(fp," %6.3f",p(d,i));
    fprintf(fp,"\n");
  }
}
