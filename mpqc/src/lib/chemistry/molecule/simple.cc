
/* simple.cc -- implementation of the simple internal coordinate classes
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

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "simple.h"
#include "simpleQCList.h"
#include "chemelem.h"
#include "localdef.h"

#include <math/topology/bitarray.h>

//////////////////////////////////////////////////////////////////////

static void add_bonds(SimpleCoList*, BitArray&, Molecule&);
static void add_bends(SimpleCoList*, BitArray&, Molecule&);
static void add_tors(SimpleCoList*, BitArray&, Molecule&);
static void add_out(SimpleCoList*, BitArray&, Molecule&);

static int linear(Molecule&,int,int,int);
static int hterminal(Molecule&, BitArray&, int);

//////////////////////////////////////////////////////////////////////

DescribedClass_REF_def(SimpleCo);

#define CLASSNAME SimpleCo
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
SimpleCo::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

SimpleCo::SimpleCo()
 : natoms_(0), ref(0), value_(0), atoms(0) {}

SimpleCo::SimpleCo(int na, const char *re) :
  natoms_(na), atoms(0), value_(0), ref(0)
{
  atoms=new int[na]; bzero(atoms,sizeof(int)*na);
  if(re) { ref=new char[strlen(re)+1]; strcpy(ref,re); }
  }

SimpleCo::~SimpleCo()
{
  init();
  }

void SimpleCo::init()
{
  if(atoms) delete[] atoms; atoms=0;
  if(ref) delete[] ref; ref=0;
  natoms_=0; value_=0;
  }

void SimpleCo::save_data_state(StateOut& so)
{
  so.put(value_);
  so.put(natoms_);
  so.put(atoms,natoms_);
  so.putstring(ref);
}

SimpleCo::SimpleCo(StateIn& si):
  SavableState(si,class_desc_)
{
  si.get(value_);
  si.get(natoms_);
  si.get(atoms);
  si.getstring(ref);
}

int SimpleCo::operator==(SimpleCo& sc)
{
  if(ref && !sc.ref || !ref && sc.ref) return 0;
  if(ref && strcmp(ref,sc.ref)) return 0;

  if(atoms && !sc.atoms || !atoms && sc.atoms) return 0;
  if(atoms)
    for(int i=0; i < natoms_; i++) if (atoms[i]!=sc.atoms[i]) return 0;

  return 1;
  }

////////////////////////////////////////////////////////////////////

SimpleCoList::SimpleCoList(KeyVal &keyval)
{
  int nsimp=keyval.count();

  for(int i=0; i < nsimp; i++) {
      add(keyval.describedclassvalue(i));
    }

}

SimpleCoList *
Geom_read_simples(KeyVal *keyval)
{
  int nsimp=keyval->count("simp");

  SimpleCoList *list = new SimpleCoList;

  for(int i=0; i < nsimp; i++) {
    char *val = keyval->pcharvalue("simp",i,0);

    if (!strcmp("stre",val))
      list->add(new Stre(keyval,"simp",i));
    else if (!strcmp("bend",val))
      list->add(new Bend(keyval,"simp",i));
    else if (!strcmp("tors",val))
      list->add(new Tors(keyval,"simp",i));
    else if (!strcmp("out",val))
      list->add(new Out(keyval,"simp",i));
    else if (!strcmp("linip",val))
      list->add(new LinIP(keyval,"simp",i));
    else if (!strcmp("linop",val))
      list->add(new LinOP(keyval,"simp",i));
    else
      err_msg("Geom_read_simples: invalid coordinate type %s",val);

    delete[] val;
    }

  return list;
  }

void Geom_print_pretty(SimpleCoList *list) { Geom_print_pretty(cout,list); }

void
Geom_print_pretty(ostream& os, SimpleCoList *list, const double *coeff)
{
  int i;
  SimpleCoListIter p;

  for(i=0,p=list; p; p++,i++) {
    os << *(p.this_object());
    if(coeff) { os.width(16); os << coeff[i]; }
    os << endl;
    }
  }

ostream& operator<<(ostream& os, SimpleCo& sc)
{
  if (sc.reference()==0 || sc.natoms()==0) return os;

  os.setf(ios::fixed,ios::floatfield); os.precision(10);

  os << "  ";
  os.setf(ios::left,ios::adjustfield); os.width(10); os << sc.reference();
  os.setf(ios::right,ios::adjustfield); os.width(6); os << sc.ctype();
  for (int i=0; i < sc.natoms(); i++) { os.width(6); os << sc[i]; }
  os.width(40-6*sc.natoms());  os << sc.preferred_value();

  return os;
  }

void
Geom_calc_simples(SimpleCoList *list, Molecule &m)
{
  for(SimpleCoListIter p=list; p; p++)
    p->calc_intco(m);
  }

/////////////////////////////////////////////////////////////////////

SimpleCoList * Geom_form_simples(Molecule& m)
{
  int i,j,ij;
  SimpleCoList *ret=0;

 // let's go through the geometry and find all the close contacts
 // bonds is a lower triangle matrix of 1's and 0's indicating whether
 // there is a bond between atoms i and j

  BitArray bonds(m.natom(),m.natom());

  for(i=0; i < m.natom(); i++) {
    double at_rad_i = m[i].element().atomic_radius();

    for(j=0; j < i; j++) {
      double at_rad_j = m[j].element().atomic_radius();

      if (bohr*dist(m[i].point(),m[j].point()) < 1.1*(at_rad_i+at_rad_j))
        bonds.set(i,j);
      }
    }

  ret = new SimpleCoList;

  add_bonds(ret,bonds,m);
  add_bends(ret,bonds,m);
  add_tors(ret,bonds,m);
  add_out(ret,bonds,m);

  return ret;
  }

/*
 * the following are translations of functions written by Gregory Humphreys
 * at the NIH
 */

/*
 * for each bonded pair, add an entry to the simple coord list
 */

static void
add_bonds(SimpleCoList *list, BitArray& bonds, Molecule& m)
{
  int i,j,ij;
  int refc=0;
  char ref[80];

  for(i=ij=0; i < m.natom(); i++) {
    for(j=0; j <= i; j++,ij++) {
      if(bonds[ij]) {
        refc++;
        sprintf(ref,"s%d",refc);
        list->add(new Stre(ref,j+1,i+1));
        }
      }
    }
  }

/*
 * for each pair of bonds sharing an atom, add a bend to the simple
 * coord list
 *
 * check each bend to see if it is linear.  if so, then we'll have to add
 * in-plane and out-of-plane linear bends as well
 *
 * let's do this later...I think I only want to do this for symmetric
 * tops, but I'm not sure...anyway, for now, let's just eliminate linear
 * bends so things like co2 won't blow up
 */

static int
linear(Molecule& m, int i, int j, int k)
{
  double dij = dist(m[i].point(),m[j].point());
  double dik = dist(m[i].point(),m[k].point());
  double djk = dist(m[j].point(),m[k].point());

 // if dij+djk==dik, then this bug is linear
  if((dij+djk - dik) < 1.0e-5) return 1;
  
  return 0;
  }

static void
add_bends(SimpleCoList *list, BitArray& bonds, Molecule& m)
{
  int i,j,k;
  int refc=0;
  char ref[80];

  int n = m.natom();

  for(i=0; i < n; i++) {
    for(j=0; j < n; j++) {
      if(bonds(i,j)) {
        for(k=0; k < i; k++) {
          if(bonds(j,k)) {
            if(linear(m,i,j,k)) continue;

	    refc++;
	    sprintf(ref,"b%d",refc);
	    list->add(new Bend(ref,k+1,j+1,i+1));
	    }
          }
	}
      }
    }
  }

/*
 * for each pair of bends which share a common bond, add a torsion
 */

/*
 * just look at the heavy-atom skeleton. return true if i is a terminal
 * atom.
 */

static int
hterminal(Molecule& m, BitArray& bonds, int i)
{
  int nh=0;
  for (int j=0; j < m.natom(); j++)
    if (bonds(i,j) && m[j].element().mass() > 1.1) nh++;
  return (nh==1);
}

static void
add_tors(SimpleCoList *list, BitArray& bonds, Molecule& m)
{
  int i,j,k,l;
  int refc=0;
  char ref[80];

  int n = m.natom();

  for(j=0; j < n; j++) {
    for(k=0; k < j; k++) {
      if(bonds(j,k)) {
        for(i=0; i < n; i++) {
          if(k==i) continue;

         // no hydrogen torsions, ok?
	  if (m[i].element().mass() < 1.1 && !hterminal(m,bonds,j)) continue;

          if (bonds(j,i)) {
	    if (linear(m,i,j,k)) continue;

            for (l=0; l < n; l++) {
              if (l==j || l==i) continue;

             // no hydrogen torsions, ok?
	      if (m[l].element().mass() < 1.1 && !hterminal(m,bonds,k))
                continue;

              if (bonds(k,l)) {
		if(linear(m,j,k,l)) continue;

		refc++;
		sprintf(ref,"t%d",refc);
		list->add(new Tors(ref,l+1,k+1,j+1,i+1));
		}
	      }
	    }
          }
	}
      }
    }
  }

static void
add_out(SimpleCoList *list, BitArray& bonds, Molecule& m)
{
  int i,j,k,l;
  int refc=0;
  char ref[80];

  int n = m.natom();

 // first find all tri-coordinate atoms
  for(i=0; i < n; i++) {
    if(bonds.degree(i)!=3) continue;

   // then look for terminal atoms connected to i
    for(j=0; j < n; j++) {
      if(bonds(i,j) && bonds.degree(j)==1) {

        for(k=0; k < n; k++) {
          if(k!=j && bonds(i,k)) {
            for(l=0; l < k; l++) {
              if(l!=j && bonds(i,l)) {
		refc++;
		sprintf(ref,"o%d",refc);
		list->add(new Out(ref,j+1,i+1,l+1,k+1));
		}
	      }
	    }
          }
	}
      }
    }
  }
