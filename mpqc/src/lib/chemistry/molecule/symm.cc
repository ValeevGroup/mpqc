
/* symm.cc -- implementation of the symmetrized internal coordinate class
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

#include "symm.h"
#include "localdef.h"
#include "simple.h"
#include "simpleQCList.h"
#include "symmQCList.h"

#include <string.h>
#include <math.h>

///////////////////////////////////////////////////////////////////////

static SymmCoList * make_symm_from_simple(SimpleCoList*, int =0);

///////////////////////////////////////////////////////////////////////

DescribedClass_REF_def(SymmCo);

#define CLASSNAME SymmCo
#define PARENTS virtual public SavableState
#define HAVE_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
SymmCo::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

SymmCo::SymmCo() :
  label_(0), n_simple(0), simple(0), coeff_(0), val(0)
{
}

SymmCo::SymmCo(const char *lab, int n) : 
  n_simple(n), val(0), label_(0), coeff_(0)
{
  if(lab) { label_ = new char[strlen(lab)+1]; strcpy(label_,lab); }
  coeff_ = new double[n_simple];
  simple = new RefSimpleCo[n_simple];
  for (int i=0; i < n_simple; i++) coeff_[i]=1.0;
  }

SymmCo::SymmCo(const char *lab, RefSimpleCo& sim) : 
  n_simple(1), val(0), label_(0), coeff_(0)
{
  if(lab) { label_ = new char[strlen(lab)+1]; strcpy(label_,lab); }
  coeff_ = new double[1];
  simple = new RefSimpleCo[1]; simple[0]=sim;
  coeff_[0]=1.0;
  }

static RefSimpleCo&
find_simple(SimpleCoList *list, const char *val)
{
  for(SimpleCoListIter p=list; p; p++) {
    if (!strcmp(p->reference(),val)) {
      return p.this_object();
      }
    }

  err_msg("cannot find simple coordinate %s\n",val);
  return 0;
  }

// sample format for SymmCo input
//  simp: (
//    st1<StreSimpleCo> = [ "st1" 1 2 ]
//    st2<StreSimpleCo> = [ "st2" 1 3 ]
//    )
//  symmco<SymmCo> = [
//       label = "asymmetric stretch"
//       simp = [ $:simp:st1 $:simp:st2 ]
//       coef = [ 1 -1 ]
//     ]
SymmCo::SymmCo(KeyVal &kv):
  n_simple(0), val(0), label_(0), coeff_(0)
{
  n_simple=kv.count("simp");

  coeff_ = new double[n_simple];
  simple = new RefSimpleCo[n_simple];
  for(int i=0; i < n_simple; i++) coeff_[i]=1.0;

  for(i=0; i < n_simple; i++) {
    simple[i] = kv.describedclassvalue("simp",i);
    if (kv.exists("coef",i)) coeff_[i]=kv.doublevalue("coef",i);
    }

  label_ = kv.pcharvalue("label");
  }

SymmCo::SymmCo(SimpleCoList *simple_list, KeyVal *kv, const char *lab, int n) :
  n_simple(0), val(0), label_(0), coeff_(0)
{
  int type=kv->count(lab,n);
  n_simple=kv->count(lab,n,1);

  coeff_ = new double[n_simple];
  simple = new RefSimpleCo[n_simple];
  for(int i=0; i < n_simple; i++) coeff_[i]=1.0;

  for(i=0; i < n_simple; i++) {
    char *v=kv->Va_pcharvalue(lab,3,n,1,i);
    if(type==3) coeff_[i]=kv->Va_doublevalue(lab,3,n,2,i);
    if(v[0]=='-') { coeff_[i]*=-1.0; v++; }

    simple[i] = find_simple(simple_list,v);
    delete[] v;
    }

  label_ = kv->pcharvalue(lab,n,0);
  }

SymmCo::SymmCo(const SymmCo& s) :
  label_(0), simple(0), coeff_(0)
{
  *this=s;
}

SymmCo::~SymmCo()
{
  init();
  }

void SymmCo::init()
{
  if(label_) delete[] label_; label_=0;
  if(coeff_) delete[] coeff_; coeff_=0;
  if(simple) delete[] simple; simple=0;
  n_simple=0;
  }

SymmCo& SymmCo::operator=(const SymmCo& s)
{
  n_simple=s.n_simple; val=s.val;

  if(label_) delete[] label_;
  label_=new char[strlen(s.label_)+1]; strcpy(label_,s.label_);

  if(coeff_) delete[] coeff_; coeff_=new double[n_simple]; 
  if(simple) delete[] simple; simple=new RefSimpleCo[n_simple];

  for (int i=0; i < n_simple; i++) {
    coeff_[i] = s.coeff_[i];
    simple[i] = s.simple[i];
    }

  return *this;
  }

int SymmCo::operator==(SymmCo& sc)
{
  if(label_ && !sc.label_ || !label_ && sc.label_) return 0;
  if(label_ && strcmp(label_,sc.label_)) return 0;

  if(n_simple != sc.n_simple) return 0;

  if(coeff_ && !sc.coeff_ || !coeff_ && sc.coeff_) return 0;
  if(coeff_)
    for(int i=0; i < n_simple; i++) if(coeff_[i]-sc.coeff_[i]>1.0e-15) return 0;

  if(simple && !sc.simple || !simple && sc.simple) return 0;
  SimpleCoListIter p1,p2;
  if(simple)
    for (int i=0; i < n_simple; i++) if (simple[i] != sc.simple[i]) return 0;

  return 1;
  }

void SymmCo::save_data_state(StateOut& so)
{
  so.putstring(label_);
  so.put(n_simple);
  so.put(val);
  so.put(coeff_,n_simple);
  for (int i=0; i < n_simple; i++) simple[i]->save_state(so);
}

SymmCo::SymmCo(StateIn& si):
  SavableState(si,class_desc_)
{
  si.getstring(label_);
  si.get(n_simple);
  si.get(val);
  si.get(coeff_);
  simple = new RefSimpleCo[n_simple];
  for (int i=0; i < n_simple; i++) {
      simple[i] = SimpleCo::restore_state(si);
    }
}

void SymmCo::print(ostream& os, const char *off) const
{
  os << off << "SymmCo:\n";

  os.setf(ios::fixed,ios::floatfield); os.precision(10);

  char *myoff = new char[strlen(off)+3]; sprintf(myoff,"%s  ",off);
  os << myoff << "label = " << label_ << endl;
  os << myoff << "val   = " << val << endl;
  os << myoff << "coeff = [ ";
  for(int i=0; i < n_simple; i++)
    os << coeff_[i] << " ";
  os << "]" << endl;

  for(i=0; i < n_simple; i++) simple[i]->print(os,myoff);
  delete[] myoff;
  os.flush();
  }

void SymmCo::print(FILE *of, const char *off) const
{
  fprintf(of,"%sSymmCo:\n",off);
  char *myoff = new char[strlen(off)+3]; sprintf(myoff,"%s  ",off);

  fprintf(of,"%slabel = %s\n",myoff,label_);
  fprintf(of,"%sval   = %20.10f\n",myoff,val);
  fprintf(of,"%scoeff = [ ",myoff);
  for(int i=0; i < n_simple; i++)
    fprintf(of,"%14.10f ",coeff_[i]);
  fprintf(of,"]\n");

  for (i=0; i < n_simple; i++) simple[i]->print(of,myoff);
  delete[] myoff;
  fflush(of);
  }

/////////////////////////////////////////////////////////////////////

/*
 * normalize the internal coordinate
 */

void SymmCo::normalize()
{
  if(!coeff_ || !n_simple) return;

  double c1=0;
  for(int i=0; i < n_simple; i++) c1 += coeff_[i]*coeff_[i];
  if(c1 < 1.0e-12) return;
  c1 = 1.0/sqrt(c1);
  for(i=0; i < n_simple; i++) coeff_[i] *= c1;
  }
  

/*
 * calculate the value of the internal coordinate.  if bmat is non-null,
 * then also calculate the row of the internal to cartesian coordinate
 * matrix for this internal coordinate.  One day I'll read Wilson, Decius
 * and Cross and explain this.
 */

void SymmCo::calc_intco(Molecule& m, double *bmat)
{
  const double bohr=0.52917706;
  int i,j,k,l;
  int a,b,c,d;
  double u[3],v[3],w[3],x[3],z[3];

  double c1=0;
  val=0;

  if(bmat) bzero(bmat,sizeof(double)*m.natom()*3);

  for(i=0; i < n_simple; i++) {
    c1 += coeff_[i]*coeff_[i];  // this ensures normalization

    val += coeff_[i] * simple[i]->calc_intco(m,bmat,coeff_[i]);
  }

  c1 = 1.0/sqrt(c1);
  val *= c1;
  if(bmat)
    for(i=0; i < m.natom()*3; i++) bmat[i] *= c1;
}

////////////////////////////////////////////////////////////////////////

ostream& operator<<(ostream& os,SymmCo& sc)
{
  int i;

  if(sc.label_==0 || sc.coeff_==0 || sc.simple==0) return os;

  os.setf(ios::fixed,ios::floatfield); os.precision(10);

  os << " ";
  os.setf(ios::left,ios::adjustfield); os.width(10); os << sc.label() << endl;

  for(i=0; i < sc.n_simple; i++) {
    os << " " << *(sc.simple[i]);
    os.width(16); os << sc.coeff_[i] << endl;
    }

  return os;
  }

///////////////////////////////////////////////////////////////////////

void Geom_print_pretty(SymmCoList *sc) { Geom_print_pretty(cout,sc); }

void
Geom_print_pretty(ostream& os, SymmCoList *sc)
{
  for(SymmCoListIter sp=sc; sp; sp++)
    os << *sp.this_object();
  }

/////////////////////////////////////////////////////////////////////////

/*
 * appends symmetrized coordinates from file pointed to by kv, and
 * referenced by the array "stype", to "slist"
 *
 * if no simple coordinates are provided, they are read from the
 * input
 */

void Geom_add_symm(KeyVal *kv, const char *stype, SymmCoList *slist,
                                                         SimpleCoList *list)
{
  int read_simp=0;

  if(list==0) { list = Geom_read_simples(kv); read_simp=1; }

  int nsymm = kv->count(stype);

  for(int i=0; i < nsymm; i++)
    slist->add(new SymmCo(list,kv,stype,i));

  if(read_simp) delete list;
  }

/*
 * read symmco array "stype" from input
 */

SymmCoList * 
Geom_read_symm(KeyVal *kv, const char *stype, SimpleCoList *list)
{
  SymmCoList *slist = new SymmCoList;
  Geom_add_symm(kv,stype,slist,list);
  return slist;
  }

/*
 * calculate the values of all the simple internal coordinates used
 * to make up the symmetrized coordinates
 */

void
Geom_calc_simples(SymmCoList *slist, Molecule& m)
{
  for(SymmCoListIter sp=slist; sp; sp++) {
    for (int i=0; i < sp->nsimple(); i++)
      sp->get_simple(i)->calc_intco(m);
  }
}

/*
 * normalize each symmco in the list
 */

void
Geom_normalize(SymmCoList *slist)
{
  for(SymmCoListIter sp=slist; sp; sp++)
    sp->normalize();
  }

/*
 * given a list of symmetrized coordinates and a geometry,
 * calculate the cartesian coordinate to internal coordinate
 * transformation matrix
 */

DMatrix
Geom_make_bmat(SymmCoList *slist, Molecule& m)
{
  int i,nsym=0;
  for(SymmCoListIter sp=slist; sp; sp++,nsym++) ;

  DMatrix bmat(nsym,m.natom()*3);

  for(i=0,sp=slist; sp; sp++,i++) {
    sp->calc_intco(m,bmat[i]);
  }

  return bmat;
  }

//////////////////////////////////////////////////////////////////////

/*
 * functions for automagically creating symmetrized coordinates, guess
 * hessians, etc...
 */

/*
 * given a list of simples, generate a list of symmetrized coordinates
 * with each symmco containing only one simpleco
 *
 * if sort==true then sort the simples by type
 */

SymmCoList *
Geom_symm_from_simple(SimpleCoList *list, int sort)
{
  if(!list) return 0;

  return make_symm_from_simple(list,sort);
  }

static SymmCoList *
make_symm_from_simple(SimpleCoList *list, int sort)
{
  SimpleCoListIter p;
  SymmCoList *slist = new SymmCoList;

  for(p=list; p; p++)
    slist->add(new SymmCo("a",p.this_object()));

  return slist;
  }

//////////////////////////////////////////////////////////////////

/*
 * given a redundant list of coordinates, form a set of unique, non-redundant
 * coordinates.  the redundant list can be either simples or symmcos
 * if just_a1 is true, then only form the totally symmetric coordinates
 */

SymmCoList *
Geom_form_symm(Molecule& m, SimpleCoList *list, int just_a1)
{
  int free_list = 0;

  if(!list) {
    if((list=Geom_form_simples(m))==0) {
      err_msg("Geom_form_symm: could not form simples list");
      return 0;
      }
    free_list=1;
    }

  SymmCoList *slist = make_symm_from_simple(list);

  if(free_list) delete list;

  SymmCoList *ret = Geom_form_symm(m,slist,just_a1);
  delete slist;

  return ret;
  }

SymmCoList *
Geom_form_symm(Molecule& m, SymmCoList *list, int just_a1)
{
  int i,j;
  SymmCoList *slist;

  if(list==0) {
    err_msg("Geom_form_symm: list is null");
    return 0;
    }

 // normalize each symmco in the list
  for(SymmCoListIter p=list; p; p++)
    p->normalize();

 // form the transformation matrix from redundant coords to non-redund.
  DMatrix K = Geom_form_K(m,list,just_a1);

  int nsym = K.ncol();
  int nred = K.nrow();

  if(nsym==0) {
    err_msg("Geom_form_symm: trouble making K matrix");
    return 0;
    }

  slist = new SymmCoList;

 // ok, now we have the K matrix, the columns of which give us the
 // contribution from each red. coord to the ith non-red. coord.
 // this gets a little hairy since the red coords can themselves be
 // linear combinations of simple coords
  for(i=0; i < nsym; i++) {
    int nonzero=0;

   // nonzero is the number of simples which will be in this coordinate
    for(j=0,p=0; j < nred; j++,p++)
      if(fabs(K[j][i]) > 1.0e-12) nonzero += p->nsimple();

    char label[80];
    sprintf(label,"gencoord%10d",i+1);
    SymmCo *sp = new SymmCo(label,nonzero);

    int sim=0;
    for(j=0,p=0; j < nred; j++,p++) {
      if(fabs(K[j][i]) > 1.0e-12) {
        for(int k=0; k < p->nsimple(); k++) {
          sp->simple[sim] = p.this_object()->simple[k];
	  sp->coeff_[sim] = K[j][i]*p.this_object()->coeff_[k];
	  sim++;
	  }
	}
      }

    slist->add(sp);
    }

  return slist;
  }

/////////////////////////////////////////////////////////////////////////
/*
 * given a list of redundant coordinates, form the B matrix for that set,
 * then form the matrix B*B~ and diagonalize it.  use the non-zero
 * eigenvectors to form the matrix K.  if just_a1 is true, then also
 * remove from K non-totally symmetric coordinates
 */

DMatrix
Geom_form_K(Molecule& m, SimpleCoList *list, int just_a1)
{
  int free_list=0;

  if(!list) {
    if((list=Geom_form_simples(m))==0) {
      err_msg("Geom_form_K: could not form simples list");
      DMatrix foo(0,0);
      return foo;
      }
    free_list=1;
    }

  SymmCoList *slist = make_symm_from_simple(list);

  if(free_list) delete list;

  DMatrix ret = Geom_form_K(m,slist,just_a1);
  delete slist;

  return ret;
  }
  

DMatrix
Geom_form_K(Molecule& m, SymmCoList *list, int just_a1)
{
  int i,j;

  if(!list) {
    err_msg("Geom_form_K: list is null");
    DMatrix foo(0,0);
    return foo;
    }
  
 // use list to form bmat
  DMatrix bmat = Geom_make_bmat(list,m);
  int nred = bmat.nrow();  // nred is the number of redundant coordinates

 // and form b*b~
  DMatrix bmbt = bmat * bmat.transpose();

 // now diagonalize bmbt, this should give you the 3n-6(5) symmetrized
 // internal coordinates
  DMatrix vecs(nred,nred);
  DVector vals(nred);

  //for(i=0; i < bmbt.nrow(); i++) bmbt(i,i) += 100.0;
  bmbt.diagonalize(vals,vecs);
  //vals.print("evals");
  //for(i=0; i < vals.dim(); i++) vals[i] -= 100.0;

  DVector coords;

  if(just_a1) {
   // ok, hopefully multiplying bmat*cart_coords will tell me which
   // coordinates are totally symmetric...we'll see
    DMatrix bm = vecs.transpose() * bmat;
    DVector geom(bmat.ncol());
    for(i=0; i < geom.dim()/3; i++) {
      geom[3*i] = m[i].point()[0];
      geom[3*i+1] = m[i].point()[1];
      geom[3*i+2] = m[i].point()[2];
      }
    coords = bm * geom;
    }

  int nonzero=0;
  if(just_a1) {
    for(i=0; i < coords.dim(); i++) if(fabs(coords[i])>1.0e-8) nonzero++;
    }
  else {
    for(i=0; i < vals.dim(); i++) if(vals[i] > 1.0e-8) nonzero++;
    }
    
  DMatrix ret(nred,nonzero);

  int coordno=0;
  for(i=0; i < nred; i++) {

   // nonzero eigenvalues are the non-redundant coordinates
    if(vals[i] > 1.0e-8) {

      int nonzero=0;
      for(j=0; j < nred; j++) if(fabs(vecs[j][i]) > 1.0e-8) nonzero++;

      if(!nonzero) {
	err_msg("Geom_form_K: no nonzero coordinates");
	DMatrix foo(0,0);
	return foo;
	}

     // if we only want the totally symmetric coords, eliminate coords
     // for which bmat[i]*cart_coords is zero

      if(!just_a1 || (fabs(coords[i]) > 1.0e-8)) {
	for(int ii=0; ii < nred; ii++)
	  ret(ii,coordno) = vecs(ii,i);
	coordno++;
        }
      }
    }

  return ret;
  }

/////////////////////////////////////////////////////////////////


// given a list of simple internal coordinates, or nothing, create an
// approximate hessian ala Fischer and Alml:of

DMatrix
Geom_form_hessian(Molecule& m, SimpleCoList *list)
{
  int free_list=0;

  if(!list) {
    if((list=Geom_form_simples(m))==0) {
      err_msg("Geom_form_hessian: could not form simples list");
      DMatrix foo(0,0);
      return foo;
      }
    free_list=1;
    }

  int count=0;
  for(SimpleCoListIter p=list; p; p++) count++;

  if(!count) {
    DMatrix hess(0,0);
    return hess;
    }

  DMatrix hess(count,count);
  hess.zero();

  int i;
  for(i=0,p=0; p; p++,i++)
    hess(i,i) = p->calc_force_con(m);

  if(free_list) delete list;

  return hess;
  }

static DMatrix
gen_inverse(DMatrix& mat)
{
  DMatrix vecs(mat.nrow(),mat.nrow());
  DVector vals(mat.nrow());

  mat.diagonalize(vals,vecs);

  DMatrix lam(mat.nrow(),mat.nrow());
  lam.zero();

  for(int i=0; i < lam.nrow(); i++)
    if(vals[i] > 1.0e-8) lam(i,i) = 1.0/vals[i];

  lam = vecs * lam * vecs.transpose();

  return lam;
  }


DMatrix
Geom_form_hessian(Molecule& m, SimpleCoList *simples, SymmCoList *symm)
{
  int i;

  DMatrix chessian;
  {
   // form the hessian in simple internal coordinates
    DMatrix ihessian = Geom_form_hessian(m,simples);

   // now form the cartesian hessian from ihessian 
    SymmCoList *rsymm = make_symm_from_simple(simples);
    DMatrix rbmat = Geom_make_bmat(rsymm,m);
    delete rsymm;

    chessian = rbmat.transpose() * ihessian * rbmat;
    }

 // and transform chessian to symm coords

  DMatrix bmat = Geom_make_bmat(symm,m);
  DMatrix bmbt = gen_inverse(bmat*bmat.transpose());

  DMatrix hessian = bmbt * bmat * chessian * bmat.transpose() * bmbt;

  return hessian;
  }
