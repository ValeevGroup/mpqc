
/* pointgrp.cc -- implementation of the point group classes
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
 *      June, 1993
 */

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <util/unix/cct_cprot.h>

#include <math/symmetry/pointgrp.h>

// this is for use by CharacterTable after making the gamma array...don't
// call free for symb or rep, just make them zero

void IrreducibleRepresentation::init()
{
  g=degen=nrot_=ntrans_=0;
  symb=0; rep=0;
}

IrreducibleRepresentation::IrreducibleRepresentation() :
  g(0), degen(0), nrot_(0), ntrans_(0), symb(0), rep(0)
{
}

IrreducibleRepresentation::IrreducibleRepresentation(
                                             int g_, int d, const char *lab)
  : g(g_), degen(d), nrot_(0), ntrans_(0), symb(0), rep(0)
{
  if(lab) { symb = new char[strlen(lab)+1]; strcpy(symb,lab); }
  if(g) { rep = new double[g]; }
}


IrreducibleRepresentation::IrreducibleRepresentation(
                                       const IrreducibleRepresentation& ir)
  : g(0), degen(0), nrot_(0), ntrans_(0), symb(0), rep(0)
{
  *this = ir;
}

IrreducibleRepresentation::~IrreducibleRepresentation()
{
  if (symb) delete[] symb;
  if (rep) delete[] rep;
  init();
}

IrreducibleRepresentation&
IrreducibleRepresentation::operator=(const IrreducibleRepresentation& ir)
{
  g = ir.g; degen = ir.degen; nrot_ = ir.nrot_; ntrans_ = ir.ntrans_;

  if (symb) delete[] symb; symb=0;
  if (ir.symb) { symb = new char[strlen(ir.symb)+1]; strcpy(symb,ir.symb); }

  if (rep) delete[] rep;
  rep = new double[g];
  for(int i=0; i < g; i++) rep[i] = ir.rep[i];

  return *this;
}

void IrreducibleRepresentation::print(FILE *fp, const char *off)
{
  if(!g) return;

  fprintf(fp,"%s%-5s",off,symb);
  for (int i=0; i < g; i++) fprintf(fp," %6.3f",rep[i]);
  fprintf(fp," | %d t, %d R\n",ntrans_,nrot_);
}

////////////////////////////////////////////////////////////////////////

SymmetryOperation::SymmetryOperation()
{
  memset(d,'\0',sizeof(double)*9);
}

SymmetryOperation::~SymmetryOperation()
{
}

void
SymmetryOperation::print(FILE* outfile) const
{
  fprintf(outfile,"        1          2          3\n");
  fprintf(outfile,"  1  %10.7f %10.7f %10.7f\n",d[0][0],d[0][1],d[0][2]);
  fprintf(outfile,"  2  %10.7f %10.7f %10.7f\n",d[1][0],d[1][1],d[1][2]);
  fprintf(outfile,"  3  %10.7f %10.7f %10.7f\n",d[2][0],d[2][1],d[2][2]);
}

////////////////////////////////////////////////////////////////////////

CharacterTable::CharacterTable()
  : g(0), nt(0), pg(C1), nirrep_(0), gamma_(0), symop(0), symb(0)
{
}

CharacterTable::CharacterTable(const CharacterTable& ct)
  : g(0), nt(0), pg(C1), nirrep_(0), gamma_(0), symop(0), symb(0)
{
  *this = ct;
}

CharacterTable::~CharacterTable()
{
  if (symb) delete[] symb; symb=0;
  if (gamma_) delete[] gamma_; gamma_=0;
  if (symop) delete[] symop; symop=0;
  g=nt=nirrep_=0;
}

CharacterTable& CharacterTable::operator=(const CharacterTable& ct)
{
  g=ct.g; nt=ct.nt; pg=ct.pg; nirrep_=ct.nirrep_;
  
  if (symb) delete[] symb; symb=0;
  if (ct.symb) { symb = new char[strlen(ct.symb)+1]; strcpy(symb,ct.symb); }

  if (gamma_) delete[] gamma_; gamma_=0;
  if (ct.gamma_) {
    gamma_ = new IrreducibleRepresentation[nirrep_];
    for (int i=0; i < nirrep_; i++) {
      gamma_[i].init();
      gamma_[i] = ct.gamma_[i];
    }
  }

  if (symop) delete[] symop; symop=0;
  if (ct.symop) {
    symop = new SymmetryOperation[g];
    for (int i=0; i < g; i++) {
      symop[i] = ct.symop[i];
    }
  }

  return *this;
}

void CharacterTable::print(FILE *fp, const char *off)
{
  if (!g || !nirrep_) return;

  fprintf(fp,"%spoint group %s\n",off,symb);

  char * myoff = new char[strlen(off)+3];
  sprintf(myoff,"%s  ",off);

  for (int i=0; i < nirrep_; i++)
    gamma_[i].print(fp,myoff);

  for(i=0; i < g; i++) symop[i].print();
}

// given matrix R in so and E in frame, perform similarity transform
// R' = ~E * R * E

static void
sim_transform(SymmetryOperation& so, const SymmetryOperation& frame)
{
  SymmetryOperation foo;

 // foo = ~E*R 
  for (int i=0; i < 3; i++) {
    for (int j=0; j < 3; j++) {
      double t=0;
      for (int k=0; k < 3; k++) t += frame(k,i)*so(k,j);
      foo(i,j) = t;
    }
  }

 // R' = ~E*R*E = foo*E
  for (i=0; i < 3; i++) {
    for (int j=0; j < 3; j++) {
      double t=0;
      for (int k=0; k < 3; k++) t += foo(i,k)*frame(k,j);
      so(i,j) = t;
    }
  }
}

CharacterTable::CharacterTable(const char *cpg, const SymmetryOperation& frame)
  : g(0), nt(0), pg(C1), nirrep_(0), gamma_(0), symop(0), symb(0)
{
  // first parse the point group symbol, this will give us the order of the
  // point group(g), the type of point group (pg), the order of the principle
  // rotation axis (nt), and the number of irreps (nirrep_)

  if (!cpg)
    err_quit("CharacterTable::CharacterTable: null point group");

  symb = new char[strlen(cpg)+1];
  for (int i=0; i < strlen(cpg); i++) symb[i] = tolower(cpg[i]);
  symb[i] = '\0';

  if (parse_symbol() < 0)
    err_quit("CharacterTable::CharacterTable: invalid point group %s",cpg);

  if (make_table() < 0)
    err_quit("CharacterTable::CharacterTable: could not make table");

  for (i=0; i < g; i++) sim_transform(symop[i],frame);
}

CharacterTable::CharacterTable(const char *cpg)
  : g(0), nt(0), pg(C1), nirrep_(0), gamma_(0), symop(0), symb(0)
{
  // first parse the point group symbol, this will give us the order of the
  // point group(g), the type of point group (pg), the order of the principle
  // rotation axis (nt), and the number of irreps (nirrep_)

  if (!cpg)
    err_quit("CharacterTable::CharacterTable: null point group");

  symb = new char[strlen(cpg)+1];
  for (int i=0; i < strlen(cpg); i++) symb[i] = tolower(cpg[i]);
  symb[i] = '\0';

  if (parse_symbol() < 0)
    err_quit("CharacterTable::CharacterTable: invalid point group %s",cpg);

  if (make_table() < 0)
    err_quit("CharacterTable::CharacterTable: could not make table");
}

int CharacterTable::parse_symbol()
{
  // default to C1 symmetry
  g=1; pg=C1; nt=1; nirrep_=1;

  if (!symb) return 0;

  if (!strcmp(symb,"c1")) return 0;

  if (!strcmp(symb,"ci")) {
    g = 2; pg = CI; nirrep_ = 2; nt = 2;
    return 0;
    }

  if(!strcmp(symb,"cs")) {
    g = 2; pg = CS; nirrep_ = 2; nt = 0;
    return 0;
    }

  if (symb[0] == 'c') {
    int nab,ne;

    if (symb[1] == '\0') return -1;

    nt = atoi(&symb[1]);
    ne = (nt%2) ? nt/2 : nt/2-1;
    nab = (nt%2) ? 1 : 2;

    if (symb[2] != '\0') {
      if(symb[2] == 'v') {
        g  = 2*nt; pg = CNV; nirrep_ = 2*nab + ne;
	}
      else if (symb[2] == 'h') {
        g  = 2*nt; pg = CNH; nirrep_ = 2*(nab+ne);
	}
      else {
        return -1;
	}
      }
    else {
      g = nt; pg = CN; nirrep_ = nab+ne;
      }

    return 0;
    }

  if (symb[0] == 'd') {
    int nab,ne;

    if(symb[1] == '\0') return -1;

    nt = atoi(&symb[1]);
    ne = (nt%2) ? nt/2 : nt/2-1;
    nab = (nt%2) ? 1 : 2;

    if (symb[2] != '\0') {
      if (symb[2] == 'd') {
        g = 4*nt; pg = DND; nirrep_ = nt+3;
        }
      else if (symb[2] == 'h') {
        g = 4*nt; pg = DNH; nirrep_ = 4*nab + 2*ne;
        }
      else {
        return -1;
        }
      }
    else {
      g = 2*nt; pg = DN; nirrep_ = 2*nab + ne;
      }

    return 0;
    }

  if (symb[0] == 's') {
    if (symb[1] == '\0') return -1;

    nt = atoi(&symb[1]);

   // only S2n groups make sense */
    if(nt%2) return -1;

    g = nt; pg = SN; nirrep_ = nt/2+1;

    return 0;
    }

  if (symb[0] == 't') {
    if (symb[1] != '\0') {
      if (symb[1] == 'd') {
        g = 24; pg = TD; nirrep_ = 5;
        }
      else if(symb[1] == 'h') {
        g = 24; pg = TH; nirrep_ = 6;
        }
      else {
        return -1;
        }
      }
    else {
      g = 12; pg = T; nirrep_ = 3;
      }

    return 0;
    }

  if (symb[0] == 'o') {
    if (symb[1] != '\0') {
      if (symb[1] == 'h') {
        pg = OH; g = 48; nirrep_ = 10;
        }
      else {
        return -1;
        }
      }
    else {
      g = 24; pg = O; nirrep_ = 5;
      }

    return 0;
    }

  if (symb[0] == 'i') {
    if (symb[1] != '\0') {
      if (symb[1] == 'h') {
        g = 120; pg = IH; nirrep_ = 10;
        }
      else {
        return -1;
        }
      }
    else {
      g = 60; pg = I; nirrep_ = 5;
      }
    }

  return -1;
}

////////////////////////////////////////////////////////////////////////

#define CLASSNAME PointGroup
#define PARENTS virtual public SavableState
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
PointGroup::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

PointGroup::PointGroup()
  : symb(0)
{
  set_symbol("c1");
  frame(0,0) = frame(1,1) = frame(2,2) = 1;
  origin_[0] = origin_[1] = origin_[2] =0;
}

PointGroup::PointGroup(const char *s)
  : symb(0)
{
  set_symbol(s);
  frame(0,0) = frame(1,1) = frame(2,2) = 1;
  origin_[0] = origin_[1] = origin_[2] =0;
}

PointGroup::PointGroup(const char *s, SymmetryOperation& so)
  : symb(0)
{
  set_symbol(s);
  frame = so;
  origin_[0] = origin_[1] = origin_[2] =0;
}

PointGroup::PointGroup(const char *s, SymmetryOperation& so, Point& or)
  : symb(0)
{
  set_symbol(s);
  frame = so;
  origin_ = or;
}

PointGroup::PointGroup(const RefKeyVal& kv)
  : symb(0)
{
  if (kv->exists("symmetry"))
    symb = kv->pcharvalue("symmetry");
  else
    set_symbol("c1");

  if (kv->exists("symmetry_frame")) {
    for (int i=0; i < 3; i++)
      for (int j=0; j < 3; j++) 
        frame(i,j) = kv->doublevalue("symmetry_frame",i,j);
  } else {
    frame(0,0) = frame(1,1) = frame(2,2) = 1;
  }

  if (kv->exists("origin")) {
    for (int i=0; i < 3; i++)
      origin_[i] = kv->doublevalue("origin",i);
  } else {
    origin_[0] = origin_[1] = origin_[2] =0;
  }
}

PointGroup::PointGroup(StateIn& si) :
  symb(0),
  SavableState(si,PointGroup::class_desc_),
  origin_(si)
{
  si.getstring(symb);
  for (int i=0; i < 3; i++)
    for (int j=0; j < 3; j++)
      si.get(frame(i,j));
}

PointGroup::PointGroup(PointGroup& pg)
  : symb(0)
{
  *this = pg;
}

PointGroup::~PointGroup()
{
  if (symb) { delete[] symb; symb=0; }
}

PointGroup& PointGroup::operator=(PointGroup& pg)
{
  set_symbol(pg.symb);
  frame = pg.frame;
  origin_ = pg.origin_;
  return *this;
}

void PointGroup::set_symbol(const char *sym)
{
  if (sym) {
    if (symb) delete[] symb;
    symb = new char[strlen(sym)+1];
    strcpy(symb,sym);
  } else {
    set_symbol("c1");
  }
}

void PointGroup::save_data_state(StateOut& so)
{
  origin_.save_object_state(so);
  so.putstring(symb);

  for (int i=0; i < 3; i++)
    for (int j=0; j < 3; j++)
      so.put(frame(i,j));
}

void PointGroup::restore_data_state(int version, StateIn& si)
{
  Point p(si); 
  origin_=p;

  si.getstring(symb);
  for (int i=0; i < 3; i++)
    for (int j=0; j < 3; j++)
      si.get(frame(i,j));
}

CharacterTable PointGroup::char_table() const
{
  CharacterTable ret(symb,frame);
  return ret;
}
